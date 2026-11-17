"""Flask web application for the NEOHELIX educational MSA demo."""

import os
import re

from flask import Flask, jsonify, render_template, request
from Bio import Entrez

from backend import (
    alignment_analysis,
    needleman_wunsch,
    progressive_star_msa,
    sequence_stats,
    smith_waterman,
)

app = Flask(__name__)
DNA_LETTERS = set("ACGTRYSWKMBDHVN")
MAX_SEQUENCE_LENGTH = 100000


def clean_sequence(sequence):
    """Clean and validate a DNA/RNA sequence."""
    cleaned = re.sub(r"\s+", "", str(sequence)).upper()
    if not cleaned:
        raise ValueError("The retrieved record did not contain a nucleotide sequence.")
    invalid = set(cleaned) - DNA_LETTERS - {"U"}
    if invalid:
        raise ValueError(
            "The retrieved record is not a supported nucleotide sequence. "
            f"Unexpected characters: {', '.join(sorted(invalid))}"
        )
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


def make_result(accession, description, raw_sequence, source):
    sequence = clean_sequence(raw_sequence)
    if len(sequence) > MAX_SEQUENCE_LENGTH:
        raise ValueError(
            f"This sequence has {len(sequence):,} bases, exceeding NEOHELIX's "
            f"{MAX_SEQUENCE_LENGTH:,}-base retrieval limit. Choose a shorter transcript or record."
        )
    return {
        "source": source,
        "accession": accession,
        "description": description,
        "sequence": sequence,
        "length": len(sequence),
        "fasta": f">{accession} {description}\n{sequence}\n",
    }


def fetch_ncbi(user_input, email, species):
    """Retrieve a nucleotide record or resolve a gene symbol/NCBI Gene ID."""
    if not email or "@" not in email:
        raise ValueError("Enter a contact email for NCBI Entrez.")

    Entrez.email = email
    Entrez.tool = "NEOHELIX_MSA_Visualizer"
    organism = (species or "Homo sapiens").strip()

    if re.match(r"^(SRR|ERR|DRR)\d+$", user_input, re.IGNORECASE):
        raise ValueError("SRA run IDs are not supported by this nucleotide-record fetcher.")
    if re.match(r"^(NP_|XP_|YP_|WP_)", user_input, re.IGNORECASE):
        raise ValueError("Protein accessions are not supported. Choose a nucleotide sequence.")

    gene_id = None
    nucleotide_id = None

    gene_match = re.match(r"^(\d+)(?:\s*[-—]\s*.*)?$", user_input)
    accession_match = re.match(
        r"^(?:[A-Z]{1,6}_?\d{5,}|[A-Z]\d{5,})(?:\.\d+)?$",
        user_input,
        re.IGNORECASE
    )

    if gene_match:
        gene_id = gene_match.group(1)
    elif accession_match:
        nucleotide_id = user_input
    else:
        term = f'{user_input}[Gene Name] AND "{organism}"[Organism]'
        with Entrez.esearch(db="gene", term=term, retmax=10) as handle:
            gene_results = Entrez.read(handle)
        gene_ids = gene_results["IdList"]
        if not gene_ids:
            raise ValueError(f"No gene named '{user_input}' was found for {organism}.")
        if len(gene_ids) > 1:
            raise ValueError(
                f"'{user_input}' matched multiple genes for {organism}. "
                "Try a more specific gene symbol or NCBI Gene ID."
            )
        gene_id = gene_ids[0]

    if gene_id:
        term = (
            f"{gene_id}[GeneID] AND biomol_mrna[PROP] "
            f"AND srcdb_refseq[PROP] AND {organism}[Organism]"
        )
        with Entrez.esearch(db="nuccore", term=term, retmax=20) as handle:
            nucleotide_results = Entrez.read(handle)
        ids = nucleotide_results["IdList"]
        if not ids:
            raise ValueError(
                f"No RefSeq mRNA records were found for Gene ID {gene_id} "
                f"in {organism}."
            )
        nucleotide_id = ids[0]

    with Entrez.efetch(
        db="nuccore", id=nucleotide_id, rettype="gb", retmode="xml"
    ) as handle:
        records = Entrez.read(handle)
    if not records:
        raise ValueError("NCBI returned no nucleotide record.")

    record = records[0]
    return make_result(
        record.get("GBSeq_accession-version", nucleotide_id),
        record.get("GBSeq_definition", "NCBI nucleotide record"),
        record.get("GBSeq_sequence", ""),
        "NCBI / GenBank",
    )


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
        result = (
            smith_waterman(seq1, seq2, match, mismatch, gap)
            if mode == "local"
            else needleman_wunsch(seq1, seq2, match, mismatch, gap)
        )
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
        rows = [
            {"name": records[i][0], "sequence": result["aligned"][i],
             "stats": sequence_stats(records[i][1])}
            for i in range(len(records))
        ]
        return jsonify({
            "rows": rows,
            "center_index": result["center_index"],
            "center_sequence": result["center_sequence"],
            "pairwise_score_totals": result["pairwise_score_totals"],
            "analysis": analysis,
        })
    except (ValueError, TypeError) as error:
        return jsonify({"error": str(error)}), 400


@app.post("/api/fetch-sequence")
def fetch_sequence():
    try:
        data = request.get_json(force=True)
        user_input = str(data.get("identifier", "")).strip()
        email = str(data.get("email", "")).strip()
        species = str(data.get("species", "Homo sapiens")).strip()

        if not user_input:
            raise ValueError("Enter a gene symbol, database ID, or sequence accession.")

        result = fetch_ncbi(user_input, email, species)

        return jsonify(result)
    except Exception as error:
        return jsonify({"error": f"Sequence retrieval failed: {error}"}), 400



if __name__ == "__main__":
    app.run(
        host="0.0.0.0",
        port=int(os.environ.get("PORT", 5000)),
        debug=False
    )
