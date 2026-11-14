# NEOHELIX — Sequence Alignment Workbench

A small educational Flask web app built around the team's Needleman–Wunsch and Smith–Waterman code. It adds a basic progressive star multiple-sequence alignment workflow.

## Included features

- Needleman–Wunsch global pairwise alignment, with score matrix and traceback path
- Smith–Waterman local pairwise alignment, with score matrix and traceback path
- Progressive star MSA for 2–8 DNA/RNA sequences
- Consensus sequence and per-position conservation
- Pairwise identity matrix for the resulting MSA
- Basic DNA composition statistics in the API
- NCBI nucleotide accession retrieval using Biopython Entrez
- Download aligned results as FASTA
- Simple HTML/CSS/JavaScript interface; no frontend framework

## Important algorithm note

NW and SW are pairwise algorithms: each aligns two sequences. The MSA feature is a separate, educational **progressive star strategy**. It selects a centre sequence using summed pairwise NW scores, aligns every other sequence to that centre, and merges insertions into shared columns. This is a simplified teaching implementation, not a production-grade MSA package such as MAFFT or MUSCLE.

The scoring model uses a match score, mismatch score and one constant gap penalty. It does not implement affine gap-open/gap-extension penalties. Input is DNA/RNA; U is converted to T. Protein alignment and substitution matrices such as BLOSUM are not implemented.

## Run locally on macOS

```bash
cd neohelix_demo
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python app.py
```

Open `http://127.0.0.1:5000` in your browser for optional local testing. This is not the public deployment URL.

## Deploy as a public website

This project uses Flask/Python, so **GitHub Pages alone cannot run it**. GitHub stores the source code; a Python-capable host runs the backend.

1. Create a GitHub repository and upload the contents of this folder (not the outer ZIP).
2. On Render, choose **New → Web Service**, connect the GitHub repository, and select the project root.
3. Use build command: `pip install -r requirements.txt`
4. Use start command: `gunicorn app:app`
5. Deploy. Render will provide a public HTTPS URL after a successful deployment.

The included `render.yaml` can help configure the service. Hosting providers may change their free-tier availability and terms.

## NCBI retrieval

The NCBI section asks for a **Nucleotide (nuccore) accession** and a contact email. Example: `NM_000546.6` (human TP53 mRNA). Biopython Entrez uses the email to identify the requester. The FASTA parser uses Biopython's `fasta-blast` format so NCBI comment lines at the beginning of a response are ignored.

**SRA accession warning:** identifiers such as `SRR501130` refer to sequencing runs in NCBI's Sequence Read Archive, not ordinary `nuccore` nucleotide records. Full SRA run reads are retrieved using SRA-specific tools such as the SRA Toolkit (`prefetch` and `fasterq-dump`) and can be very large. This small hosted demo intentionally does not download entire SRA runs. Use a nucleotide-record accession such as `NM_000546.6` for the NCBI fetch demonstration. See the [NCBI SRA download guide](https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/) and [NCBI E-utilities documentation](https://www.ncbi.nlm.nih.gov/books/NBK25499/) for the distinction.

## Suggested viva explanation

- `needleman_wunsch`: initializes the first row/column with gap penalties, fills each cell using diagonal/up/left scores, then traces back from the bottom-right.
- `smith_waterman`: uses the same recurrence but includes zero as a candidate, remembers the highest-scoring cell, and traces back until a zero cell.
- `progressive_star_msa`: chooses a centre, aligns every sequence to it with NW, then merges insertions into shared slots.
- `alignment_analysis`: calculates a majority consensus, column conservation and pairwise identity from aligned rows.
- `app.py`: validates input and exposes the Python functions through Flask routes used by the browser.
- `static/app.js`: sends input to Flask and renders returned results. The alignment algorithms remain in Python.

## Project structure

```text
neohelix_demo/
├── app.py
├── backend.py
├── templates/
│   └── index.html
├── static/
│   ├── app.js
│   └── style.css
├── requirements.txt
├── Procfile
├── render.yaml
└── README.md
```
