# NEOHELIX — Sequence Alignment Workbench

A small educational Flask web app built around the team's Needleman–Wunsch and Smith–Waterman code. It adds a basic progressive star multiple-sequence alignment workflow with a colour-coded alignment grid.

## Included features

- Needleman–Wunsch global pairwise alignment, with score matrix and traceback path
- Smith–Waterman local pairwise alignment, with score matrix and traceback path
- Progressive star MSA for 2–8 DNA/RNA sequences
- **Colour-coded alignment grid** with three views (bases, conservation, differences)
- Consensus sequence and per-position conservation
- Pairwise identity matrix for the resulting MSA
- Basic DNA composition statistics in the API
- NCBI nucleotide accession retrieval using Biopython Entrez
- Download aligned results as FASTA
- Simple HTML/CSS/JavaScript interface; no frontend framework

## Important algorithm note

NW and SW are pairwise algorithms: each aligns two sequences. The MSA feature is a separate, educational **progressive star strategy**. It selects a centre sequence using summed pairwise NW scores, aligns every other sequence to that centre, and merges insertions into shared columns. This is a simplified teaching implementation, not a production-grade MSA package such as MAFFT or MUSCLE.

The scoring model uses a match score, mismatch score and one constant gap penalty. It does not implement affine gap-open/gap-extension penalties. Input is DNA/RNA; U is converted to T. Protein alignment and substitution matrices such as BLOSUM are not implemented.

## Colour-coded alignment grid

After a multiple alignment runs, the result is drawn as a colour-coded grid by `static/msa_grid.js`. The module has no dependencies and injects its own CSS, so it does not need changes to `style.css`.

**Views** (switch with the buttons above the grid):

- **Bases:** each cell is coloured by nucleotide: green for A, red for T/U, amber for G, blue for C. Gaps are dimmed.
- **Conservation:** cells are brighter where the column is more conserved across the sequences.
- **Differences:** cells that match the consensus show as dots; only bases that differ are coloured and outlined.

**Other elements:**

- A position ruler above the alignment, a conservation bar chart and the consensus row below it.
- Hovering a cell shows the sequence name, the alignment column, the residue number within the original sequence, the column conservation and the consensus base.
- A cell-size slider for zooming long alignments. Sequence names stay fixed while the grid scrolls sideways.

**How it connects to the rest of the app:** `static/app.js` receives `rows` and `analysis` from `/api/msa` and calls:

```js
MSAGrid.render(byId("msaGrid"), { rows: data.rows, analysis: data.analysis });
```

`index.html` loads `msa_grid.js` before `app.js`. If `analysis` is missing, the module recomputes the consensus and conservation using the same definition as `alignment_analysis` in `backend.py` (percentage of all sequences sharing the most common base in that column; gaps count in the denominator).

**Limitation:** the grid draws one element per cell, which is fine for the 8-sequence limit and a few thousand alignment columns. Very long alignments (for example, several full-length NCBI genes) would need a length cap or a canvas-based version.

## Run locally

Requires Python 3.10 or newer. Run these commands from the folder that contains `app.py` (if you cloned the repository, that is the repository folder).

### Windows (PowerShell)

```powershell
cd neohelix_demo
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
python app.py
```

If PowerShell refuses to activate the environment, run `Set-ExecutionPolicy -Scope CurrentUser RemoteSigned` once and activate again. In Command Prompt, activate with `.venv\Scripts\activate.bat` instead.

### macOS / Linux

```bash
cd neohelix_demo
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python app.py
```

On macOS, port 5000 may already be used by AirPlay Receiver. If so, run `PORT=5001 python app.py` and open port 5001. On Windows PowerShell, set it with `$env:PORT=5001; python app.py`.

### Open the app

Open `http://127.0.0.1:5000` in your browser for local testing. This is not the public deployment URL. After editing `app.js`, `msa_grid.js` or `style.css`, refresh the page with Ctrl+Shift+R (Cmd+Shift+R on macOS) so the browser does not use cached files. Changes to Python files need the server to be restarted.

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
- `static/msa_grid.js`: draws the returned alignment as a colour-coded grid (bases, conservation and difference views) in the browser. It only displays results; it does not align anything.

## Project structure

```text
neohelix_demo/
├── app.py
├── backend.py
├── templates/
│   └── index.html
├── static/
│   ├── app.js
│   ├── msa_grid.js
│   └── style.css
├── requirements.txt
├── Procfile
├── render.yaml
└── README.md
```
