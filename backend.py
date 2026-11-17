"""Core sequence-alignment and analysis functions for NEOHELIX.

Educational implementation: simple scoring, constant gap penalty, and
progressive star alignment. It is not a replacement for production MSA tools.
"""

import numpy as np


def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2, include_visual_data=True):
    """Return a global alignment of two complete sequences."""
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    matrix = np.zeros((rows, cols), dtype=int)
    traceback = np.zeros((rows, cols), dtype=int)  # 0 diagonal, 1 up, 2 left

    # A leading gap of length k costs k * gap.
    for i in range(rows):
        matrix[i][0] = i * gap
        if i > 0:
            traceback[i][0] = 1
    for j in range(cols):
        matrix[0][j] = j * gap
        if j > 0:
            traceback[0][j] = 2

    # Fill the dynamic-programming matrix.
    for i in range(1, rows):
        for j in range(1, cols):
            substitution = match if seq1[i - 1] == seq2[j - 1] else mismatch
            diagonal = matrix[i - 1][j - 1] + substitution
            up = matrix[i - 1][j] + gap
            left = matrix[i][j - 1] + gap
            best = max(diagonal, up, left)
            matrix[i][j] = best
            traceback[i][j] = [diagonal, up, left].index(best)

    # Trace back from the bottom-right cell to recover both aligned strings.
    aligned1, aligned2 = "", ""
    i, j = rows - 1, cols - 1
    path = [(i, j)]
    while i > 0 or j > 0:
        if i > 0 and j > 0 and traceback[i][j] == 0:
            aligned1 = seq1[i - 1] + aligned1
            aligned2 = seq2[j - 1] + aligned2
            i -= 1
            j -= 1
        elif i > 0 and traceback[i][j] == 1:
            aligned1 = seq1[i - 1] + aligned1
            aligned2 = "-" + aligned2
            i -= 1
        else:
            aligned1 = "-" + aligned1
            aligned2 = seq2[j - 1] + aligned2
            j -= 1
        path.append((i, j))

        result = {
        "score": int(matrix[-1][-1]),
        "aligned_seq1": aligned1,
        "aligned_seq2": aligned2,
    }

    if include_visual_data:
        result["matrix"] = matrix.tolist()
        result["path"] = list(reversed(path))

    return result


def smith_waterman(seq1, seq2, match=1, mismatch=-1, gap=-2, include_visual_data=True):
    """Return the highest-scoring local alignment between two sequences."""
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    matrix = np.zeros((rows, cols), dtype=int)
    traceback = np.zeros((rows, cols), dtype=int)  # 0 diagonal, 1 up, 2 left, 3 stop
    best_score, best_i, best_j = 0, 0, 0

    for i in range(1, rows):
        for j in range(1, cols):
            substitution = match if seq1[i - 1] == seq2[j - 1] else mismatch
            diagonal = matrix[i - 1][j - 1] + substitution
            up = matrix[i - 1][j] + gap
            left = matrix[i][j - 1] + gap
            best = max(0, diagonal, up, left)
            matrix[i][j] = best
            if best == 0:
                traceback[i][j] = 3
            else:
                traceback[i][j] = [diagonal, up, left].index(best)
            if best > best_score:
                best_score, best_i, best_j = best, i, j

    aligned1, aligned2 = "", ""
    i, j = best_i, best_j
    end_i, end_j = i, j
    path = [(i, j)]
    while i > 0 and j > 0 and traceback[i][j] != 3:
        direction = traceback[i][j]
        if direction == 0:
            aligned1 = seq1[i - 1] + aligned1
            aligned2 = seq2[j - 1] + aligned2
            i -= 1
            j -= 1
        elif direction == 1:
            aligned1 = seq1[i - 1] + aligned1
            aligned2 = "-" + aligned2
            i -= 1
        else:
            aligned1 = "-" + aligned1
            aligned2 = seq2[j - 1] + aligned2
            j -= 1
        path.append((i, j))

        result = {
        "score": int(best_score),
        "aligned_seq1": aligned1,
        "aligned_seq2": aligned2,
        "start1": i, "end1": end_i,
        "start2": j, "end2": end_j,
    }

    if include_visual_data:
        result["matrix"] = matrix.tolist()
        result["path"] = list(reversed(path))

    return result


def _alignment_slots(aligned_center, aligned_other):
    """Split an aligned pair into insertions around each center base."""
    insertions = [""]
    aligned_bases = []
    center_index = 0
    for center_char, other_char in zip(aligned_center, aligned_other):
        if center_char == "-":
            insertions[center_index] += other_char
        else:
            aligned_bases.append(other_char)
            center_index += 1
            insertions.append("")
    return insertions, aligned_bases


def progressive_star_msa(sequences, match=1, mismatch=-1, gap=-2):
    """Build a simple star MSA by aligning every sequence to a chosen center.

    The center is the sequence with the highest sum of pairwise NW scores.
    Insertions are merged into shared slots around center-sequence positions.
    """
    if len(sequences) < 2:
        raise ValueError("Enter at least two sequences for alignment.")
    if len(sequences) > 8:
        raise ValueError("This educational demo supports up to 8 sequences.")

    count = len(sequences)
    totals = [0] * count
    for i in range(count):
        for j in range(i + 1, count):
            score = needleman_wunsch(sequences[i], sequences[j], match, mismatch, gap)["score"]
            totals[i] += score
            totals[j] += score
    center_index = totals.index(max(totals))
    center = sequences[center_index]

    pair_data = {}
    max_insertions = [0] * (len(center) + 1)
    for index, sequence in enumerate(sequences):
        if index == center_index:
            continue
        result = needleman_wunsch(center, sequence, match, mismatch, gap)
        insertions, aligned_bases = _alignment_slots(result["aligned_seq1"], result["aligned_seq2"])
        pair_data[index] = (insertions, aligned_bases)
        for slot, inserted in enumerate(insertions):
            max_insertions[slot] = max(max_insertions[slot], len(inserted))

    aligned_rows = []
    for index, sequence in enumerate(sequences):
        pieces = []
        if index == center_index:
            for position, base in enumerate(center):
                pieces.append("-" * max_insertions[position])
                pieces.append(base)
            pieces.append("-" * max_insertions[-1])
        else:
            insertions, aligned_bases = pair_data[index]
            for position, base in enumerate(aligned_bases):
                pieces.append(insertions[position] + "-" * (max_insertions[position] - len(insertions[position])))
                pieces.append(base)
            pieces.append(insertions[-1] + "-" * (max_insertions[-1] - len(insertions[-1])))
        aligned_rows.append("".join(pieces))

    return {
        "aligned": aligned_rows,
        "center_index": center_index,
        "center_sequence": center,
        "pairwise_score_totals": totals,
    }


def alignment_analysis(aligned_sequences):
    """Calculate consensus, per-column conservation and pairwise identity."""
    if not aligned_sequences:
        raise ValueError("No aligned sequences were supplied.")
    width = len(aligned_sequences[0])
    if any(len(row) != width for row in aligned_sequences):
        raise ValueError("All aligned sequences must have equal length.")

    consensus = ""
    conservation = []
    for column in range(width):
        residues = [row[column] for row in aligned_sequences if row[column] != "-"]
        if not residues:
            consensus += "-"
            conservation.append(0)
        else:
            counts = {base: residues.count(base) for base in set(residues)}
            most_common = max(counts, key=counts.get)
            consensus += most_common
            conservation.append(round(counts[most_common] / len(aligned_sequences) * 100, 1))

    identities = []
    for i in range(len(aligned_sequences)):
        row = []
        for j in range(len(aligned_sequences)):
            compared = list(zip(aligned_sequences[i], aligned_sequences[j]))
            # A gap opposite a nucleotide is treated as a non-match.
            matches = sum(a == b and a != "-" for a, b in compared)
            row.append(round(matches / len(compared) * 100, 2) if compared else 0)
        identities.append(row)

    return {"consensus": consensus, "conservation": conservation, "identity": identities}


def sequence_stats(sequence):
    """Return basic DNA composition statistics (excluding gaps)."""
    sequence = sequence.replace("-", "")
    length = len(sequence)
    counts = {base: sequence.count(base) for base in "ACGTN"}
    gc = counts["G"] + counts["C"]
    at = counts["A"] + counts["T"]
    return {
        "length": length,
        "counts": counts,
        "gc_percent": round(gc / length * 100, 2) if length else 0,
        "at_percent": round(at / length * 100, 2) if length else 0,
    }
