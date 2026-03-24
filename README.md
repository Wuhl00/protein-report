protein-report — Quick Start

Overview
  protein-report is a Python tool for end-to-end protein sequence analysis.
  It takes a protein FASTA sequence and generates a bookmarked PDF report plus a Markdown report.

What it does
  - Physicochemical properties (local, Biopython ProtParam)
  - Hydropathy plot (Kyte–Doolittle)
  - Domain analysis (EBI InterProScan; async submit/poll/fetch)
  - Homology search (EBI BLASTP against SwissProt/Reviewed; 180s timeout with graceful degradation)
  - Structured AI summary (English): Investigation Summary / Functional Prediction / Related Literature Search (PubMed link)
  - PDF sidebar bookmarks (PyPDF2)

Minimal files needed (for non-Trae users)
  - protein_analyzer.py
  - requirements.txt
  - input.fasta (your input; can be replaced each run)
  - README.txt (this file)

Installation
  1) Create and activate a Python environment (Python ≥3.8, Recommended: 3.10+).
  2) Install dependencies:
       pip install -r requirements.txt

Input
  Put your protein sequence in FASTA format into input.fasta.
  You can replace input.fasta content for each analysis run.

Example:
  >MyProtein
  MAA....(amino acids)....KKK

Run
  From the project root:
    python protein_analyzer.py

Output
  Each run creates a separate output folder to avoid overwriting:
    analysis_runs\<FASTA_ID>_YYYYMMDD_HHMMSS\

  Inside the folder you will find:
    - <FASTA_ID>_report.pdf
    - <FASTA_ID>_report.md
    - hydrophobicity.png
    - domain_map.png

Notes on network dependency
  - InterProScan and BLAST rely on EBI web services.
  - Transient network errors are retried.
  - BLAST has a hard timeout (180 seconds). If it times out, the report is still generated.

Troubleshooting
  - If PDF bookmarks do not appear, confirm PyPDF2 is installed and rerun.
  - If BLAST hits are missing, it may have timed out; rerun later or use the NCBI BLAST portal link in the report.
