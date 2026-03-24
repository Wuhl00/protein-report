---
name: protein-report
description: Protein sequence analysis Skill. Takes a protein FASTA sequence and automatically runs physicochemical analysis, hydropathy plotting, EBI InterProScan domain analysis, EBI BLAST homology search, and a structured AI summary. Outputs a bookmarked PDF and a Markdown report (one output folder per run).
---

# Protein Sequence Deep Analysis Skill (protein-report)

This Skill is designed for protein sequences and turns a multi-step bioinformatics workflow into a single input action. Provide a protein FASTA sequence and get a ready-to-read PDF report (with sidebar bookmarks) plus a Markdown report (easy to edit/share).

## Features

- **Physicochemical properties (local)**: Biopython ProtParam metrics including length, molecular weight, pI, instability index, aromaticity, and GRAVY.
- **Hydropathy plot**: Kyte-Doolittle hydropathy profile.
- **Domain analysis (EBI InterProScan)**: Asynchronous submit/poll/fetch; extracts domain locations and generates a domain map plus a sorted domain table.
- **Homology search (EBI BLASTP / SwissProt)**: Asynchronous BLASTP against SwissProt/Reviewed; generates clickable UniProt links; includes a 180s timeout and graceful degradation (never blocks the full report).
- **Structured AI summary (English)**: Investigation Summary / Functional Prediction / Related Literature Search (PubMed search link only).
- **Outputs**:
  - PDF: One-stop report, sidebar bookmarks + clickable external links.
  - Markdown: Easy to edit/copy/share (images referenced via relative paths).

## Tech Stack

- **Parsing & analysis**: `Biopython` (FASTA parsing, ProtParam metrics)
- **Plotting**: `Matplotlib` (hydropathy + domain map)
- **Domain search**: `EBI InterProScan REST API` (async submit/poll/fetch)
- **Homology search**: `EBI NCBI BLAST REST API` (async submit/poll/fetch; fallback to text parsing when needed)
- **PDF reporting**: `fpdf` (PDF generation) + `PyPDF2` (write sidebar bookmarks/outline)

## Usage

- **Input**: A protein sequence in standard FASTA format (single-letter amino acids).
- **How to run**:
  1. Put the sequence into `input.fasta` at the project root
  2. Run `python protein_analyzer.py`
- **Output location**: One output folder per run to avoid overwriting `.png/.pdf/.md`:
  - `analysis_runs/<FASTA_ID>_YYYYMMDD_HHMMSS/`
- **Network dependency**:
  - InterProScan and BLAST rely on EBI services. Transient network errors are retried; timeouts gracefully degrade without blocking the entire report.

## Example

### 1) End-to-end analysis
User input:
```fasta
>Sample_Protein
MAVSRSSRLRLGRALAAAAAATAVALPAVAVAGPPAVAAAAA
```
Outputs:
- `analysis_runs/Sample_Protein_YYYYMMDD_HHMMSS/Sample_Protein_report.pdf`
- `analysis_runs/Sample_Protein_YYYYMMDD_HHMMSS/Sample_Protein_report.md`
- The same folder also contains `hydrophobicity.png` and `domain_map.png`
