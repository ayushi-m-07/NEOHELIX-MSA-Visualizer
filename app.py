"""Flask web application for the NEOHELIX educational MSA demo."""

import os
import re
from flask import Flask, jsonify, render_template, request
from Bio import Entrez, SeqIO
from io import StringIO

from backend import (
    alignment_analysis,
    needleman_wunsch,
    progressive_star_msa,
    sequence_stats,
    smith_waterman,
)

app = Flask(__name__)
DNA_LETTERS = set("ACGTRYSWKMBDHVN")

def clean_sequence(sequence):
    """Clean and validate a DNA/RNA sequence."""

    cleaned = re.sub(r"\s+", "", sequence).upper()

    if not cleaned:
        raise ValueError("Sequence cannot be empty.")

    invalid = set(cleaned) - DNA_LETTERS - {"U"}

    if invalid:
        raise ValueError(
            f"Invalid nucleotide characters: {', '.join(sorted(invalid))}"
        )

    # Normalize RNA uracil to DNA thymine
    return cleaned.replace("U", "T")

def parse_sequences(text):
    """Read FASTA text or one sequence per line."""
    records = []
    current_name = None
    current_sequence = []
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        raise ValueError("Paste at least two sequences.")

    if any(line.startswith(">") for line in lines):
        for line in lines:
            if line.startswith(">"):
                if current_name is not None:
                    records.append((current_name, clean_sequence("".join(current_sequence))))
                current_name = line[1:].strip() or "Sequence " + str(len(records) + 1)
                current_sequence = []
            else:
                if current_name is None:
                    raise ValueError("FASTA sequence data must follow a header beginning with >.")
                current_sequence.append(line)
        if current_name is not None:
            records.append((current_name, clean_sequence("".join(current_sequence))))
    else:
        for index, line in enumerate(lines, start=1):
            records.append(("Sequence " + str(index), clean_sequence(line)))

    if len(records) < 2:
        raise ValueError("MSA needs at least two sequences.")
    if len(records) > 8:
        raise ValueError("This demo supports up to 8 sequences.")
    return records


@app.get("/")
def home():
    return render_template("index.html")


@app.get("/health")
def health():
    return jsonify({"status": "ok"})


@app.post("/api/pairwise")
def pairwise():
    try:
        data = request.get_json(force=True)
        seq1 = clean_sequence(data.get("seq1", ""))
        seq2 = clean_sequence(data.get("seq2", ""))
        mode = data.get("mode", "global")
        match = int(data.get("match", 1))
        mismatch = int(data.get("mismatch", -1))
        gap = int(data.get("gap", -2))
        if mode == "local":
            result = smith_waterman(seq1, seq2, match, mismatch, gap)
        else:
            result = needleman_wunsch(seq1, seq2, match, mismatch, gap)
        return jsonify({"mode": mode, "result": result})
    except (ValueError, TypeError) as error:
        return jsonify({"error": str(error)}), 400


@app.post("/api/msa")
def msa():
    try:
        data = request.get_json(force=True)
        records = parse_sequences(data.get("sequences", ""))
        sequences = [sequence for _, sequence in records]
        match = int(data.get("match", 1))
        mismatch = int(data.get("mismatch", -1))
        gap = int(data.get("gap", -2))
        result = progressive_star_msa(sequences, match, mismatch, gap)
        analysis = alignment_analysis(result["aligned"])
        rows = [{"name": records[i][0], "sequence": result["aligned"][i],
                 "stats": sequence_stats(records[i][1])} for i in range(len(records))]
        return jsonify({"rows": rows, "center_index": result["center_index"],
                        "center_sequence": result["center_sequence"],
                        "pairwise_score_totals": result["pairwise_score_totals"],
                        "analysis": analysis})
    except (ValueError, TypeError) as error:
        return jsonify({"error": str(error)}), 400




@app.post("/api/ncbi")
def ncbi_fetch():
    try:
        data = request.get_json(force=True)

        user_input = data.get("accession", "").strip()
        email = data.get("email", "").strip()

        if not user_input:
            raise ValueError(
                "Enter an NCBI Gene ID, gene symbol, or nucleotide accession."
            )

        if not email or "@" not in email:
            raise ValueError("Enter a valid contact email for NCBI Entrez.")

        Entrez.email = email
        Entrez.tool = "NEOHELIX_MSA_Visualizer"

        # Reject SRA run identifiers
        if re.match(r"^(SRR|ERR|DRR)\d+$", user_input, re.IGNORECASE):
            raise ValueError(
                "SRA run IDs are not supported. Please enter a gene ID, "
                "gene symbol, or nucleotide accession instead."
            )

        # Reject protein accessions
        if re.match(r"^(NP_|XP_|YP_|WP_)", user_input, re.IGNORECASE):
            raise ValueError(
                "Protein accessions are not supported. "
                "NEOHELIX currently retrieves nucleotide sequences only."
            )

        gene_id = None
        nucleotide_id = None

        # CASE 1: Numeric Gene ID, e.g. 672 or 672 - BRCA1
        gene_match = re.match(
            r"^\s*(\d+)(?:\s*[-—]\s*.*)?$",
            user_input
        )

        if gene_match:
            gene_id = gene_match.group(1)

        # CASE 2: Nucleotide accession, e.g. NM_007294.4
        elif re.match(
            r"^(?:[A-Z]{1,6}_?\d{5,}|[A-Z]\d{5,})(?:\.\d+)?$",
            user_input,
            re.IGNORECASE
        ):
            nucleotide_id = user_input

        # CASE 3: Gene symbol, e.g. BRCA1
        else:
            with Entrez.esearch(
                db="gene",
                term=f'{user_input}[Gene Name] AND "Homo sapiens"[Organism]',
                retmax=3
            ) as handle:
                gene_results = Entrez.read(handle)

            gene_ids = gene_results["IdList"]

            if not gene_ids:
                raise ValueError(
                    f"No human gene found for '{user_input}'. "
                    "Check the gene symbol or Gene ID."
                )

            if len(gene_ids) > 1:
                raise ValueError(
                    f"Multiple human genes matched '{user_input}'. "
                    "Please enter the specific NCBI Gene ID."
                )

            gene_id = gene_ids[0]

        # Resolve Gene ID to a RefSeq mRNA record
        if gene_id is not None:
            search_term = (
                f"{gene_id}[GeneID] AND biomol_mrna[PROP] "
                f"AND srcdb_refseq[PROP] AND Homo sapiens[Organism]"
            )

            with Entrez.esearch(
                db="nuccore",
                term=search_term,
                retmax=20
            ) as handle:
                nucleotide_results = Entrez.read(handle)

            nucleotide_ids = nucleotide_results["IdList"]

            if not nucleotide_ids:
                raise ValueError(
                    f"No human RefSeq mRNA nucleotide records found "
                    f"for Gene ID {gene_id}."
                )

            # Use the first matching RefSeq mRNA record
            nucleotide_id = nucleotide_ids[0]

        # Retrieve GenBank XML rather than decoding FASTA text
        with Entrez.efetch(
            db="nuccore",
            id=nucleotide_id,
            rettype="gb",
            retmode="xml"
        ) as handle:
            records = Entrez.read(handle)

        if not records:
            raise ValueError("NCBI returned no nucleotide record.")

        record = records[0]

        # Extract sequence and metadata from GenBank XML
       
        raw_sequence = record["GBSeq_sequence"]
        sequence = clean_sequence(raw_sequence)

        resolved_accession = record["GBSeq_accession-version"]
        description = record["GBSeq_definition"]

        # Maintain the existing NEOHELIX sequence length limit
        if len(sequence) > 100000:
            raise ValueError(
                "This record is longer than the demo's "
                "100,000-base limit. Choose a shorter nucleotide record."
            )

        # Keep the same response structure expected by app.js
        return jsonify({
            "accession": resolved_accession,
            "description": description,
            "sequence": sequence,
            "length": len(sequence),
            "fasta": (
                ">" + resolved_accession + " " + description
                + "\n" + sequence + "\n"
            )
        })

    except Exception as error:
        return jsonify({
            "error": "NCBI retrieval failed: " + str(error)
        }), 400



if __name__ == "__main__":
    # Local development only. Public hosting is configured separately on Render.
    app.run(host="0.0.0.0", port=int(os.environ.get("PORT", 5000)), debug=False)
