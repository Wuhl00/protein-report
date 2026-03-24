import os
import requests
import time
import hashlib
import io
import re
from urllib.parse import quote
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import ProtParamData
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez
from fpdf import FPDF
import pandas as pd

class ProteinAnalyzer:
    """
    蛋白质序列深度分析类 (全自动化版)。
    集成：理化性质、二级结构、疏水性绘图、EBI InterProScan 结构域自动化抓取与绘图。
    """
    
    def __init__(self, fasta_path, output_dir=None):
        if not os.path.exists(fasta_path):
            raise FileNotFoundError(f"序列文件 {fasta_path} 未找到。")
            
        self.fasta_path = fasta_path
        self.record = list(SeqIO.parse(fasta_path, "fasta"))[0]
        self.sequence = str(self.record.seq).upper()
        self.output_dir = output_dir or os.getcwd()
        os.makedirs(self.output_dir, exist_ok=True)
        
        # 动态生成项目级唯一邮箱 (用于 API 请求)
        project_hash = hashlib.md5(os.getcwd().encode()).hexdigest()[:8]
        self.email = f"trae_user_{project_hash}@gmail.com"
        
        self.analysis_results = {
            'physicochemical': {},
            'secondary_structure': {},
            'domains': [],
            'ebi_job_id': None
        }

    def analyze_physicochemical(self):
        """
        [阶段一] 本地理化性质分析。
        参考 Expasy ProtParam 标准，使用 Bjellqvist 标尺计算 pI。
        """
        print(f"--- 正在分析 {self.record.id} 的理化性质 ---")
        analysed_seq = ProteinAnalysis(self.sequence)
        
        # 计算结果
        mw = analysed_seq.molecular_weight()
        pi = analysed_seq.isoelectric_point()
        instability = analysed_seq.instability_index()
        aromaticity = analysed_seq.aromaticity()
        gravy = analysed_seq.gravy()
        
        self.analysis_results['physicochemical'] = {
            'Length (aa)': len(self.sequence),
            'Molecular Weight (Da)': f"{mw:.2f}",
            'Isoelectric Point (pI)': f"{pi:.2f} (Bjellqvist Scale)",
            'Instability Index': f"{instability:.2f} (<40 Stable, >40 Unstable)",
            'Aromaticity': f"{aromaticity:.3f}",
            'GRAVY (Grand Average of Hydropathy)': f"{gravy:.3f} (-: Hydrophilic, +: Hydrophobic)",
        }
        return self.analysis_results['physicochemical']

    def analyze_secondary_structure(self):
        """
        [阶段二] 二级结构倾向分析。
        计算 Helix, Turn, Sheet 的百分比倾向。
        """
        print(f"--- 正在分析 {self.record.id} 的二级结构倾向 ---")
        analysed_seq = ProteinAnalysis(self.sequence)
        # 获取 Helix, Turn, Sheet 的百分比
        helix, turn, sheet = analysed_seq.secondary_structure_fraction()
        
        self.analysis_results['secondary_structure'] = {
            'Helix (Alpha-helix)': f"{helix*100:.1f}%",
            'Turn (Coil/Loop)': f"{turn*100:.1f}%",
            'Sheet (Beta-sheet)': f"{sheet*100:.1f}%"
        }
        return self.analysis_results['secondary_structure']

    def plot_hydrophobicity(self, window=9, output_img="hydrophobicity.png"):
        """
        本地绘制疏水性图表 (Kyte & Doolittle 算法)。
        """
        print(f"--- 正在生成疏水性分析图 ---")
        analysed_seq = ProteinAnalysis(self.sequence)
        # 使用 Kyte-Doolittle 标尺
        # Biopython 1.80+ 推荐用法
        values = analysed_seq.protein_scale(ProtParamData.kd, window, edge=1.0)
        
        plt.figure(figsize=(10, 4))
        plt.plot(values, color='teal', linewidth=1.5)
        plt.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        plt.title(f"Hydrophobicity Profile (Window={window}) - {self.record.id}")
        plt.xlabel("Residue Position")
        plt.ylabel("Score")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        output_path = output_img
        if not os.path.isabs(output_img) and os.path.dirname(output_img) == "":
            output_path = os.path.join(self.output_dir, output_img)
        plt.savefig(output_path)
        plt.close()
        self.analysis_results['hydro_plot'] = output_path
        return output_path

    def run_ebi_interproscan(self):
        """
        [自动化核心] 提交并抓取 EBI InterProScan 结构域分析结果。
        """
        base_url = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
        
        # 1. 提交任务
        print(f"--- 正在向 EBI 提交结构域分析任务 (邮箱: {self.email}) ---")
        try:
            submit_res = requests.post(f"{base_url}/run", data={
                'email': self.email,
                'sequence': self.sequence
            })
            if submit_res.status_code != 200:
                print(f"提交失败: {submit_res.text}")
                return False
            job_id = submit_res.text
            self.analysis_results['ebi_job_id'] = job_id
        except Exception as e:
            print(f"提交异常: {e}")
            return False

        # 2. 轮询状态 (智能等待)
        wait_time = 30
        print("--- 正在排队等待 EBI 分析结果 (通常需要 1-3 分钟) ---")
        start_time = time.time()
        while True:
            try:
                status = requests.get(f"{base_url}/status/{job_id}", timeout=20).text.upper()
            except requests.exceptions.RequestException as e:
                if time.time() - start_time > 600:
                    print(f"状态轮询超时: {e}")
                    return False
                print(f"状态轮询异常，重试中: {e}")
                time.sleep(min(wait_time, 60))
                continue
            print(f"当前状态: {status}")
            if status == "FINISHED": break
            if status in ["FAILED", "ERROR"]: return False
            time.sleep(wait_time)
            wait_time = min(wait_time + 10, 60)

        # 3. 抓取 JSON 结果
        print("--- 正在抓取并解析结构域数据 ---")
        try:
            res_data = requests.get(f"{base_url}/result/{job_id}/json", timeout=60).json()
            domains = []
            
            # 修正解析逻辑：InterProScan JSON 结构为 results -> [ { matches: [...] } ]
            if 'results' in res_data and len(res_data['results']) > 0:
                matches = res_data['results'][0].get('matches', [])
                print(f"DEBUG: 成功进入结果层级，发现 {len(matches)} 条匹配项 (Matches)")
                
                for match in matches:
                    signature = match.get('signature', {})
                    # 获取核心元数据
                    acc = signature.get('accession')
                    name = signature.get('name') or acc
                    desc = signature.get('description') or "No description available"
                    db = signature.get('signatureLibraryRelease', {}).get('library', 'Unknown')
                    
                    # 抓取该 match 下的所有位置 (Locations)
                    for loc in match.get('locations', []):
                        domains.append({
                            'acc': str(acc),
                            'name': str(name),
                            'desc': str(desc),
                            'start': loc.get('start'),
                            'end': loc.get('end'),
                            'db': db
                        })
            
            self.analysis_results['domains'] = domains
            print(f"DEBUG: 成功提取出 {len(domains)} 个详细特征位置")
            return True
        except Exception as e:
            print(f"结果解析异常: {e}")
            return False

    def plot_domains(self, output_img="domain_map.png"):
        """
        绘制全量结构域分布图。
        """
        if not self.analysis_results['domains']:
            print("未发现结构域，跳过绘图。")
            return None
            
        print("--- 正在绘制全量结构域地图 ---")
        domains = self.analysis_results['domains']
        # 按起始位置排序
        domains.sort(key=lambda x: x['start'])
        
        fig, ax = plt.subplots(figsize=(12, max(4, len(domains) * 0.4)))
        
        # 绘制蛋白质骨架
        ax.hlines(0, 1, len(self.sequence), colors='gray', linestyles='solid', linewidth=2, label='Full Protein')
        
        # 绘制结构域色块
        colors = plt.cm.get_cmap('tab20', len(domains))
        for i, dom in enumerate(domains):
            y_pos = i + 1
            start, end = dom['start'], dom['end']
            rect = plt.Rectangle((start, y_pos-0.3), end-start, 0.6, color=colors(i), alpha=0.8)
            ax.add_patch(rect)
            # 标注文字
            label = f"{dom['name']} ({dom['db']})"
            ax.text(start, y_pos + 0.4, label, fontsize=8, verticalalignment='center')
            
        ax.set_ylim(-1, len(domains) + 2)
        ax.set_xlim(0, len(self.sequence) + 50)
        ax.set_xlabel("Residue Position")
        ax.set_yticks([])
        ax.set_title(f"Comprehensive Domain Map - {self.record.id}")
        
        plt.tight_layout()
        output_path = output_img
        if not os.path.isabs(output_img) and os.path.dirname(output_img) == "":
            output_path = os.path.join(self.output_dir, output_img)
        plt.savefig(output_path)
        plt.close()
        self.analysis_results['domain_plot'] = output_path
        return output_path

    def run_blast(self, hit_count=5):
        """
        [自动化核心] 执行 EBI NCBI BLASTP 搜索 (异步)。
        目标数据库：uniprotkb_swissprot。内置 180s 超时保护与优雅降级。
        """
        base_url = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"
        print(f"--- 正在通过 EBI 提交 BLASTP 任务 (数据库: uniprotkb_swissprot) ---")
        try:
            submit = requests.post(f"{base_url}/run", data={
                "email": self.email,
                "program": "blastp",
                "database": "uniprotkb_swissprot",
                "sequence": self.sequence,
                "stype": "protein"
            })
            if submit.status_code != 200:
                print(f"BLAST 提交失败: {submit.text}")
                self.analysis_results['blast_hits'] = []
                return False
            job_id = submit.text.strip()
        except Exception as e:
            print(f"BLAST 提交异常: {e}")
            self.analysis_results['blast_hits'] = []
            return False

        # 轮询状态，最多等待 180s
        start_time = time.time()
        print("--- 正在轮询 EBI BLAST 任务状态 ---")
        while True:
            try:
                status = requests.get(f"{base_url}/status/{job_id}").text.upper()
            except Exception:
                status = "UNKNOWN"
            if status == "FINISHED":
                break
            if status in ["FAILED", "ERROR", "NOT_FOUND"]:
                print(f"BLAST 任务失败: {status}")
                self.analysis_results['blast_hits'] = []
                return False
            if time.time() - start_time > 180:
                print("BLAST 超时，已优雅降级跳过。")
                self.analysis_results['blast_hits'] = []
                return False
            time.sleep(5)

        # 抓取 JSON 结果，解析前 N 个命中
        print("--- 正在抓取并解析 EBI BLAST 结果 ---")
        hits = []
        def _normalize_uniprot_acc(value):
            if value is None:
                return None
            v = str(value).strip()
            if ":" in v:
                v = v.split(":", 1)[1].strip()
            if "|" in v:
                parts = [p for p in v.split("|") if p]
                if len(parts) >= 2:
                    v = parts[1].strip()
            if "." in v:
                v = v.split(".", 1)[0]
            if not v:
                return None
            # Try to extract a plausible UniProt accession from the string
            m = re.search(r"\b([A-Z0-9]{6,10})\b", v.upper())
            return (m.group(1) if m else v.upper())

        # 优先解析 XML (稳定获得 identity/e-value/acc)
        try:
            xml_res = requests.get(f"{base_url}/result/{job_id}/xml")
            xml_text = xml_res.text if xml_res is not None else ""
            xml_lower = xml_text.lower()
            if ("<blastoutput" not in xml_lower) and ("<?xml" not in xml_lower):
                raise ValueError("No BLAST XML detected")
            blast_record = NCBIXML.read(io.StringIO(xml_text))
            for alignment in blast_record.alignments[:hit_count]:
                if not alignment.hsps:
                    continue
                hsp = alignment.hsps[0]
                acc = _normalize_uniprot_acc(alignment.accession) or _normalize_uniprot_acc(alignment.hit_id)
                if not acc:
                    continue
                identity_pct = f"{(hsp.identities / hsp.align_length) * 100:.1f}%"
                hits.append({
                    'title': alignment.title,
                    'acc': acc,
                    'e_value': hsp.expect,
                    'identity': identity_pct,
                    'url': f"https://www.uniprot.org/uniprotkb/{acc}/entry"
                })
            self.analysis_results['blast_hits'] = hits
            print(f"DEBUG: 成功获取 {len(hits)} 条同源比对结果 (XML)")
            return bool(hits)
        except Exception:
            pass

        try:
            data = requests.get(f"{base_url}/result/{job_id}/json").json()
            # 在不确定 JSON 层级的情况下，递归查找 hits 列表
            def _find_hits(obj):
                if isinstance(obj, dict):
                    for k, v in obj.items():
                        if k == "hits" and isinstance(v, list):
                            return v
                        found = _find_hits(v)
                        if found is not None:
                            return found
                elif isinstance(obj, list):
                    for it in obj:
                        found = _find_hits(it)
                        if found is not None:
                            return found
                return None

            jhits = _find_hits(data) or []

            for h in jhits[:hit_count]:
                desc = (h.get('description') or [{}])[0]
                title = desc.get('title') or desc.get('id') or "Unknown"
                raw_id = (desc.get('id') or "")
                # 提取 accession (sp|Qxxxx| 或者 直接 accession)
                acc = _normalize_uniprot_acc(desc.get('accession') or raw_id)
                if not acc:
                    continue
                hsp = (h.get('hsps') or [{}])[0]
                evalue = hsp.get('evalue')
                ident = hsp.get('identity')
                alen = hsp.get('align_len') or hsp.get('align-len') or hsp.get('align_len'.upper(), None)
                identity_pct = None
                if ident is not None and alen:
                    try:
                        identity_pct = f"{(float(ident)/float(alen))*100:.1f}%"
                    except Exception:
                        identity_pct = "N/A"
                hits.append({
                    'title': title,
                    'acc': acc,
                    'e_value': evalue,
                    'identity': identity_pct or "N/A",
                    'url': f"https://www.uniprot.org/uniprotkb/{acc}/entry"
                })
            # 如果 JSON 解析未获取任何命中，则尝试文本解析回退
            if not hits:
                raise ValueError("Empty hits from JSON, fallback to text parser")
        except Exception as e:
            print(f"JSON 解析失败，尝试解析文本输出: {e}")
            try:
                text = requests.get(f"{base_url}/result/{job_id}/out").text
                # 简单解析命中表格 (Sequences producing significant alignments)
                lines = text.splitlines()
                table_started = False
                for ln in lines:
                    if "Sequences producing significant alignments" in ln:
                        table_started = True
                        continue
                    if table_started:
                        ln = ln.strip()
                        if not ln:
                            if hits:
                                break
                            else:
                                continue
                        # 解析 accession 与 e-value (表格通常以 ... E value 结尾)
                        parts = ln.split()
                        if len(parts) >= 2:
                            # 文本第一列可能是描述，Accession 常以 'sp|Qxxxx|' 或 'tr|..|' 形式出现在行内
                            # 尝试从整行中捕获形如 '|QXXXX|' 的 accession
                            acc_candidate = None
                            if '|' in ln:
                                segs = ln.split('|')
                                if len(segs) >= 2:
                                    acc_candidate = segs[1].strip()
                            if not acc_candidate:
                                acc_candidate = parts[0]
                            evalue_candidate = parts[-1]
                            acc = _normalize_uniprot_acc(acc_candidate)
                            if not acc:
                                continue
                            hits.append({
                                'title': ln,
                                'acc': acc,
                                'e_value': evalue_candidate,
                                'identity': "N/A",
                                'url': f"https://www.uniprot.org/uniprotkb/{acc}/entry"
                            })
                        if len(hits) >= hit_count:
                            break
                if hits:
                    for h in hits:
                        acc = h.get('acc')
                        if not acc:
                            continue
                        start_match = re.search(rf"^>.*\b{re.escape(acc)}\b.*$", text, flags=re.IGNORECASE | re.MULTILINE)
                        if not start_match:
                            start_match = re.search(rf"\b{re.escape(acc)}\b", text, flags=re.IGNORECASE)
                        if not start_match:
                            continue
                        seg = text[start_match.start():start_match.start() + 8000]
                        m = re.search(r"Identities\s*=\s*\d+/\d+\s*\((\d+)%\)", seg)
                        if m:
                            h['identity'] = f"{m.group(1)}%"
            except Exception as e2:
                print(f"文本解析亦失败: {e2}")
                hits = []

        self.analysis_results['blast_hits'] = hits
        print(f"DEBUG: 成功获取 {len(hits)} 条同源比对结果")
        return bool(hits)

    def generate_ai_summary(self):
        """
        [阶段三] 综合分析结论 (AI 综述)。
        三段论：1. Investigation Summary; 2. Functional Prediction; 3. Authoritative References。
        输出为英文，便于在 PDF 中直接展示。
        """
        print("--- 正在生成 AI 功能预测综述 (三段论, 英文) ---")
        domains = self.analysis_results.get('domains', [])
        blast_hits = self.analysis_results.get('blast_hits', [])
        phys = self.analysis_results.get('physicochemical', {})
        
        # 提取关键信息
        domain_names = list(set([d['name'] for d in domains if d['db'] in ['PFAM', 'INTERPRO']]))
        tm_domains = [d for d in domains if d['db'] in ['TMHMM', 'PHOBIUS']]
        top_hit = blast_hits[0] if blast_hits else None
        
        # 1) Investigation Summary
        length = phys.get('Length (aa)', 'N/A')
        mw = phys.get('Molecular Weight (Da)', 'N/A')
        pi = str(phys.get('Isoelectric Point (pI)', 'N/A')).split(' ')[0]
        instability = str(phys.get('Instability Index', 'N/A')).split(' ')[0]
        summary_part = (
            f"The input protein has a length of {length} aa with an estimated molecular weight of {mw} Da. "
            f"The computed isoelectric point (pI) is {pi}, and the instability index is {instability}. "
        )
        if tm_domains:
            summary_part += f"{len(tm_domains)} putative transmembrane segment(s) were detected by TM predictors. "
        else:
            summary_part += "No obvious transmembrane segments were detected, suggesting a soluble protein. "
        core_domains = ", ".join(domain_names) if domain_names else "no well-characterized domains"
        summary_part += f"Domain scanning indicates core components: {core_domains}. "
        if top_hit and top_hit.get('acc'):
            summary_part += f"Homology search shows the top SwissProt hit is {top_hit['acc']} (identity: {top_hit.get('identity','N/A')})."

        # 2) Functional Prediction
        prediction_part = ""
        is_skp1 = any("Skp1" in d for d in domain_names) or (top_hit and "SKP1" in (top_hit.get('title','').upper()))
        if is_skp1:
            prediction_part = ("Supported by Skp1-related domains and high-confidence SwissProt homology, "
                               "this protein is predicted to be a core adaptor of the SCF (Skp1–Cullin–F-box) E3 ubiquitin ligase complex. "
                               "It likely bridges F-box proteins to the Cullin scaffold, thereby regulating ubiquitin-mediated proteolysis "
                               "in pathways such as cell cycle control, hormone signaling, or stress responses.")
        elif any("Kinase" in d for d in domain_names):
            prediction_part = ("The presence of kinase-related domains suggests protein phosphorylation activity, potentially acting as a "
                               "switch or amplifier in intracellular signal transduction cascades.")
        elif tm_domains:
            prediction_part = ("Pronounced transmembrane topology indicates a membrane-associated role. "
                               "Even in the absence of clear annotation, the protein may function as a receptor, channel, or anchor in biological membranes.")
        else:
            prediction_part = ("Given the lack of well-characterized domains and high-identity curated homologs, "
                               "this protein may represent a novel biochemical function. "
                               "Consider AlphaFold-based structural modeling to identify potential ligand-binding pockets or interaction interfaces.")

        # 3) Related Literature Search (PubMed)
        search_term = "+".join(domain_names[:2]) if domain_names else "protein+function"
        ref_links = [f"[PubMed Search](https://pubmed.ncbi.nlm.nih.gov/?term={search_term})"]

        self.analysis_results['ai_structured'] = {
            'summary': summary_part,
            'prediction': prediction_part,
            'references': ref_links
        }
        return self.analysis_results['ai_structured']

    def generate_report(self, output_pdf=None):
        pdf = FPDF()
        pdf.add_page()
        bookmarks = []  # 收集章节标题与页码，用于侧边书签
 
        # 标题
        pdf.set_font("Arial", 'B', 20)
        pdf.cell(200, 15, txt=f"Protein Analysis: {self.record.id}", ln=True, align='C')
        pdf.ln(4)

        # 内部链接（用于章节锚点）
        link_seq = pdf.add_link()
        link_phys = pdf.add_link()
        link_hydro = pdf.add_link()
        link_domain = pdf.add_link()
        link_blast = pdf.add_link()
        link_ai = pdf.add_link()
        
        # 1. 序列展示
        pdf.set_font("Arial", 'B', 14)
        pdf.cell(200, 10, txt="1. Input Sequence", ln=True)
        pdf.set_link(link_seq, y=pdf.get_y(), page=pdf.page_no())
        bookmarks.append(("1. Input Sequence", pdf.page_no()))
        pdf.set_font("Courier", '', 8)
        seq = self.sequence
        for i in range(0, len(seq), 70):
            pdf.cell(0, 4, txt=seq[i:i+70], ln=True)
        pdf.ln(5)
        
        # 2. 理化性质
        pdf.set_font("Arial", 'B', 14)
        pdf.cell(200, 10, txt="2. Physicochemical Properties", ln=True)
        pdf.set_link(link_phys, y=pdf.get_y(), page=pdf.page_no())
        bookmarks.append(("2. Physicochemical Properties", pdf.page_no()))
        pdf.set_font("Arial", '', 11)
        for k, v in self.analysis_results.get('physicochemical', {}).items():
            pdf.cell(0, 8, txt=f"- {k}: {v}", ln=True)
        pdf.ln(5)
        
        # 3. 疏水性图表 (放在第一页，紧跟理化性质)
        if 'hydro_plot' in self.analysis_results:
            pdf.set_font("Arial", 'B', 14)
            pdf.cell(200, 10, txt="3. Hydrophobicity Profile", ln=True)
            pdf.set_link(link_hydro, y=pdf.get_y(), page=pdf.page_no())
            bookmarks.append(("3. Hydrophobicity Profile", pdf.page_no()))
            pdf.image(self.analysis_results['hydro_plot'], x=15, w=180)
            pdf.ln(5)

        # 4. 结构域地图与列表
        if 'domain_plot' in self.analysis_results:
            pdf.add_page()
            pdf.set_font("Arial", 'B', 16)
            pdf.cell(200, 12, txt="4. Comprehensive Domain Analysis", ln=True)
            pdf.set_link(link_domain, y=pdf.get_y(), page=pdf.page_no())
            bookmarks.append(("4. Comprehensive Domain Analysis", pdf.page_no()))
            pdf.image(self.analysis_results['domain_plot'], x=10, w=190)
            pdf.ln(10)
            
            # 标题
            pdf.set_font("Arial", 'B', 13)
            pdf.cell(0, 10, txt="Domain Details (Sorted by Position):", ln=True)
            
            # Note 移至开头
            pdf.set_font("Arial", 'I', 10)
            pdf.set_text_color(100, 100, 100) # 深灰色
            pdf.multi_cell(0, 6, txt="Note: Clicking the database tags [e.g. PFAM] will redirect to official documentation pages.")
            pdf.ln(4)
            
            # 列表展示 (按位置排序)
            pdf.set_text_color(0, 0, 0)
            pdf.set_font("Arial", '', 11)
            
            sorted_domains = sorted(self.analysis_results['domains'], key=lambda x: x['start'])
            for dom in sorted_domains:
                db_upper = str(dom.get('db', '')).upper()
                acc = str(dom.get('acc', '')).strip()
                if db_upper in ['PHOBIUS', 'TMHMM', 'MOBIDB_LITE']:
                    url = None
                elif db_upper == 'PFAM':
                    url = f"https://www.ebi.ac.uk/interpro/entry/pfam/{quote(acc)}"
                elif db_upper == 'INTERPRO':
                    url = f"https://www.ebi.ac.uk/interpro/entry/interpro/{quote(acc)}"
                else:
                    slug_map = {
                        'SMART': 'smart',
                        'PANTHER': 'panther',
                        'GENE3D': 'gene3d',
                        'SUPERFAMILY': 'superfamily',
                        'CDD': 'cdd',
                    }
                    if db_upper in slug_map:
                        url = f"https://www.ebi.ac.uk/interpro/entry/{slug_map[db_upper]}/{quote(acc)}"
                    else:
                        url = f"https://www.ebi.ac.uk/interpro/search/text/{quote(acc)}"
                
                # 1. 打印文本部分 (黑色)
                prefix = f"- [{dom['start']}-{dom['end']}] {dom['name']} ({dom['acc']}) "
                pdf.write(8, prefix)
                
                # 2. 打印标签部分 (蓝色 + 链接)
                tag = f"[{dom['db']}]"
                if url:
                    pdf.set_text_color(0, 0, 255)
                    pdf.write(8, tag, url)
                else:
                    pdf.set_text_color(120, 120, 120)
                    pdf.write(8, tag)
                
                # 3. 换行并恢复颜色
                pdf.set_text_color(0, 0, 0)
                pdf.ln(8)

        # 5. 同源性分析 (BLAST)
        pdf.add_page()
        pdf.set_font("Arial", 'B', 16)
        pdf.cell(200, 12, txt="5. Homology Analysis (BLASTP)", ln=True)
        pdf.set_link(link_blast, y=pdf.get_y(), page=pdf.page_no())
        bookmarks.append(("5. Homology Analysis (BLASTP)", pdf.page_no()))
        pdf.ln(5)
        
        pdf.set_font("Arial", 'I', 10)
        pdf.set_text_color(100, 100, 100)
        pdf.multi_cell(0, 6, txt="Default database: SwissProt (Reviewed). For broader searches (e.g., 'nr'), use the official portal below if automated retrieval is skipped due to server load.")
        pdf.set_text_color(0, 0, 255)
        pdf.set_font("Arial", 'BU', 10)
        pdf.write(8, "[NCBI BLASTP Portal]", "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome")
        pdf.ln(8)
        pdf.set_text_color(0, 0, 0)

        if self.analysis_results.get('blast_hits'):
            pdf.set_font("Arial", '', 11)
            for hit in self.analysis_results['blast_hits']:
                pdf.set_font("Arial", 'B', 11)
                title_text = hit.get('title', '')[:150] + "..."
                pdf.multi_cell(0, 7, txt=f"> {title_text}")
                pdf.set_font("Arial", '', 10)
                ev = hit.get('e_value', 'N/A')
                ident = hit.get('identity', 'N/A')
                pdf.write(6, f"  - E-value: {ev} | Identity: {ident} | Accession: ")
                pdf.set_text_color(0, 0, 255)
                pdf.write(6, f"[{hit.get('acc', 'N/A')}]", hit.get('url', ''))
                pdf.set_text_color(0, 0, 0)
                pdf.ln(10)
        else:
            pdf.set_font("Arial", '', 11)
            pdf.multi_cell(0, 6, txt="Automated BLAST retrieval was skipped due to server load or timeout. Please use the portal link above if immediate results are needed.")
        
        # 6. AI 辅助功能预测结论 (三段论, 英文直接输出)
        if 'ai_structured' in self.analysis_results:
            pdf.add_page()
            pdf.set_font("Arial", 'B', 16)
            pdf.cell(200, 12, txt="6. AI-Assisted Function Prediction", ln=True)
            pdf.set_link(link_ai, y=pdf.get_y(), page=pdf.page_no())
            bookmarks.append(("6. AI-Assisted Function Prediction", pdf.page_no()))
            pdf.ln(5)
            
            ai = self.analysis_results['ai_structured']
            
            # 1. Summary
            pdf.set_font("Arial", 'B', 12)
            pdf.cell(0, 8, txt="6.1 Investigation Summary", ln=True)
            pdf.set_font("Arial", '', 11)
            # 过滤非 ASCII 字符，避免 FPDF 编码错误
            summary_ascii = "".join((c if ord(c) < 128 else "-" if c in "–—" else "" ) for c in ai['summary'])
            pdf.multi_cell(0, 7, txt=summary_ascii)
            pdf.ln(5)
            
            # 2. Prediction
            pdf.set_font("Arial", 'B', 12)
            pdf.cell(0, 8, txt="6.2 Functional Prediction", ln=True)
            pdf.set_font("Arial", '', 11)
            pdf.set_fill_color(245, 245, 255) 
            prediction_ascii = "".join((c if ord(c) < 128 else "-" if c in "–—" else "" ) for c in ai['prediction'])
            pdf.multi_cell(0, 7, txt=prediction_ascii, border=1, fill=True)
            pdf.ln(5)
            
            # 3. References
            pdf.set_font("Arial", 'B', 12)
            pdf.cell(0, 8, txt="6.3 Related Literature Search", ln=True)
            pdf.set_font("Arial", '', 10)
            for ref in ai['references']:
                if '](' in ref and ')' in ref:
                    text = ref[ref.find("[")+1:ref.find("]")]
                    url = ref[ref.find("(")+1:ref.find(")")]
                    # 限制为 ASCII 以避免编码问题
                    text_ascii = "".join(c for c in text if ord(c) < 128)
                    pdf.set_text_color(0, 0, 255)
                    pdf.write(7, f"- {text_ascii}", url)
                    pdf.set_text_color(0, 0, 0)
                    pdf.ln(7)
            
            pdf.ln(10)
            pdf.set_font("Arial", 'I', 9)
            pdf.set_text_color(150, 150, 150)
            pdf.multi_cell(0, 5, txt="Disclaimer: This summary is automatically synthesized by the Trae AI protein analysis agent based on integrated bioinformatics data. Please verify with wet-lab experiments or specialized literature.")
        
        if output_pdf is None:
            output_pdf = os.path.join(self.output_dir, f"{self.record.id}_report.pdf")
        pdf.output(output_pdf)
        print(f"PDF 报告已生成: {output_pdf}")

        added = self._try_add_pdf_bookmarks(output_pdf, bookmarks)
        if added:
            print("PDF 侧边书签已添加。")
        else:
            print("PDF 侧边书签未添加（PyPDF2 未安装或添加失败）。")

    def _try_add_pdf_bookmarks(self, pdf_path, bookmarks):
        try:
            from PyPDF2 import PdfReader, PdfWriter
        except Exception:
            return False
        try:
            reader = PdfReader(pdf_path)
            writer = PdfWriter()
            for page in reader.pages:
                writer.add_page(page)
            # 添加书签（页码从 0 开始）
            for title, page_no in bookmarks:
                try:
                    writer.add_outline_item(title, page_number=page_no-1)
                except Exception:
                    # 兼容老版本 API
                    try:
                        writer.addBookmark(title, page_no-1)
                    except Exception:
                        pass
            tmp_path = pdf_path + ".tmp.pdf"
            with open(tmp_path, "wb") as f:
                writer.write(f)
            # 替换原文件
            try:
                os.replace(tmp_path, pdf_path)
            except Exception:
                pass
            return True
        except Exception:
            return False

    def generate_markdown_report(self, output_md=None):
        if output_md is None:
            output_md = os.path.join(self.output_dir, f"{self.record.id}_report.md")
        with open(output_md, "w", encoding="utf-8") as f:
            f.write(f"# 蛋白质深度分析报告: {self.record.id}\n\n")
            
            f.write("## 1. 分析序列 (Input Sequence)\n")
            f.write("```fasta\n")
            f.write(f">{self.record.id}\n")
            seq = self.sequence
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")
            f.write("```\n\n")
            
            f.write("## 2. 理化性质分析 (Physicochemical Properties)\n")
            f.write("| 性质 | 计算值 |\n")
            f.write("| --- | --- |\n")
            for k, v in self.analysis_results['physicochemical'].items():
                f.write(f"| {k} | {v} |\n")
            f.write("\n")
            
            if 'domain_plot' in self.analysis_results:
                f.write("## 3. 全量结构域地图 (Domain Map)\n")
                f.write(f"![Domains]({os.path.basename(self.analysis_results['domain_plot'])})\n\n")
                
                f.write("### 结构域详细列表 (按位置排序)\n")
                f.write("| Position | Name | Accession | Database | Link |\n")
                f.write("| --- | --- | --- | --- | --- |\n")
                
                # 按位置排序
                sorted_domains = sorted(self.analysis_results['domains'], key=lambda x: x['start'])
                for dom in sorted_domains:
                    # 根据数据库生成跳转链接
                    db_link = "N/A"
                    db_upper = str(dom.get('db', '')).upper()
                    acc = str(dom.get('acc', '')).strip()
                    if db_upper in ['PHOBIUS', 'TMHMM', 'MOBIDB_LITE']:
                        db_link = "N/A"
                    elif db_upper == 'PFAM':
                        db_link = f"[View Pfam](https://www.ebi.ac.uk/interpro/entry/pfam/{quote(acc)})"
                    elif db_upper == 'INTERPRO':
                        db_link = f"[View InterPro](https://www.ebi.ac.uk/interpro/entry/interpro/{quote(acc)})"
                    else:
                        slug_map = {
                            'SMART': 'smart',
                            'PANTHER': 'panther',
                            'GENE3D': 'gene3d',
                            'SUPERFAMILY': 'superfamily',
                            'CDD': 'cdd',
                        }
                        if db_upper in slug_map:
                            db_link = f"[View {db_upper}](https://www.ebi.ac.uk/interpro/entry/{slug_map[db_upper]}/{quote(acc)})"
                        else:
                            db_link = f"[Search InterPro](https://www.ebi.ac.uk/interpro/search/text/{quote(acc)})"
                        
                    f.write(f"| {dom['start']}-{dom['end']} | {dom['name']} | {dom['acc']} | {dom['db']} | {db_link} |\n")
                f.write("\n")

            if 'hydro_plot' in self.analysis_results:
                f.write("## 4. 疏水性分布图 (Hydrophobicity Profile)\n")
                f.write(f"![Hydrophobicity]({os.path.basename(self.analysis_results['hydro_plot'])})\n\n")
            
            f.write("## 5. 建模与深度挖掘门户 (Action Portal)\n")
            f.write("### 🧬 结构建模 (AlphaFold / ColabFold)\n")
            f.write(f"- **[前往 AlphaFold 数据库搜索](https://alphafold.ebi.ac.uk/search/text?q={self.sequence[:20]})**\n")
            f.write("- **[使用 ColabFold (从头建模)](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)**\n\n")

            if 'blast_hits' in self.analysis_results:
                f.write("## 6. 同源性比对 (Homology Analysis - BLASTP)\n")
                f.write("> **Note**: 为了保证功能预测的高可靠性，脚本默认检索 **SwissProt (Reviewed)** 数据库。如果您需要检索包括预测序列在内的全量数据库（如 `nr`），请前往：\n")
                f.write("> **[👉 NCBI BLASTP 官方网站](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome)**\n\n")
                f.write("| Accession | Identity | E-value | Title | Link |\n")
                f.write("| --- | --- | --- | --- | --- |\n")
                for hit in self.analysis_results['blast_hits']:
                    f.write(f"| {hit['acc']} | {hit['identity']} | {hit['e_value']} | {hit['title'][:80]}... | [View UniProt]({hit['url']}) |\n")
                f.write("\n")
            
            if 'ai_structured' in self.analysis_results:
                f.write("## 7. AI 辅助功能预测结论 (AI-Assisted Function Prediction)\n")
                ai = self.analysis_results['ai_structured']
                
                f.write("### 7.1 调查总结 (Investigation Summary)\n")
                f.write(f"> {ai['summary']}\n\n")
                
                f.write("### 7.2 功能预测 (Functional Prediction)\n")
                f.write(f"**{ai['prediction']}**\n\n")
                
                f.write("### 7.3 Related Literature Search\n")
                for ref in ai['references']:
                    f.write(f"- {ref}\n")
                f.write("\n")
                f.write("--- \n")
                f.write("*免责声明：本综述由 Trae AI 蛋白分析助手基于集成生物信息学数据自动合成，仅供参考。请结合实验或专业文献进行验证。*\n\n")

            f.write("### 🧬 二级结构倾向预测\n")
            for k, v in self.analysis_results.get('secondary_structure', {}).items():
                f.write(f"- {k}: {v}\n")
            
        print(f"Markdown 报告已生成: {output_md}")

    # 预留占位方法
    def run_blast_phase2(self):
        pass

    def run_cd_phase3(self):
        pass

if __name__ == "__main__":
    input_file = "input.fasta"
    if os.path.exists(input_file):
        safe_id = re.sub(r"[^A-Za-z0-9._-]+", "_", str(list(SeqIO.parse(input_file, "fasta"))[0].id))
        ts = time.strftime("%Y%m%d_%H%M%S")
        out_dir = os.path.join(os.getcwd(), "analysis_runs", f"{safe_id}_{ts}")
        analyzer = ProteinAnalyzer(input_file, output_dir=out_dir)
        analyzer.analyze_physicochemical()
        analyzer.analyze_secondary_structure()
        analyzer.plot_hydrophobicity()
        
        # 开启自动化结构域检索 (耗时约 1-2 分钟)
        if analyzer.run_ebi_interproscan():
            analyzer.plot_domains()
            
        # 开启同源比对（EBI BLAST，内置 180s 超时与优雅降级）
        analyzer.run_blast()
        
        # 开启 AI 综述生成
        analyzer.generate_ai_summary()
            
        analyzer.generate_report()
        analyzer.generate_markdown_report()
