"""
alignment_cli_demo.py

Pure Python, terminal-based demo of Needleman-Wunsch (global) and
Smith-Waterman (local) pairwise alignment. No HTML, CSS, or JS —
just run this file directly.

Prints:
    - the aligned sequences, color-coded by nucleotide
    - the full DP score matrix
    - the traceback path highlighted within that matrix

Usage:
    python alignment_cli_demo.py
    (press Enter at the prompts to use the built-in example sequences)
"""

import numpy as np

# ---------------------------------------------------------------------
# Terminal colors (ANSI escape codes) — purely cosmetic, safe to ignore
# if your terminal doesn't render them; the algorithm output is the
# same either way.
# ---------------------------------------------------------------------
COLORS = {
    "A": "\033[92m",   # green
    "T": "\033[91m",   # red
    "U": "\033[91m",   # red (RNA)
    "G": "\033[93m",   # yellow
    "C": "\033[94m",   # blue
    "-": "\033[90m",   # grey
}
RESET = "\033[0m"
PATH_HIGHLIGHT = "\033[1;46;30m"  # bold, cyan background, black text


def colorize_base(ch: str) -> str:
    color = COLORS.get(ch.upper(), "")
    return f"{color}{ch}{RESET}" if color else ch


# ---------------------------------------------------------------------
# Alignment algorithms
# ---------------------------------------------------------------------
def needleman_wunsch(seq1: str, seq2: str, match: int = 1,
                      mismatch: int = -1, gap: int = -2) -> dict:
    """Global alignment. Forces the full length of both sequences in."""
    m, n = len(seq1) + 1, len(seq2) + 1
    matrix = np.zeros((m, n))
    traceback = np.zeros((m, n), dtype=int)  # 0=diag, 1=up, 2=left

    for i in range(m):
        matrix[i][0] = i * gap
    for j in range(n):
        matrix[0][j] = j * gap
    for i in range(1, m):
        traceback[i][0] = 1
    for j in range(1, n):
        traceback[0][j] = 2

    for i in range(1, m):
        for j in range(1, n):
            score = match if seq1[i - 1] == seq2[j - 1] else mismatch
            diagonal = matrix[i - 1][j - 1] + score
            up = matrix[i - 1][j] + gap
            left = matrix[i][j - 1] + gap
            best = max(diagonal, up, left)
            matrix[i][j] = best
            traceback[i][j] = [diagonal, up, left].index(best)

    aligned1, aligned2, path = "", "", []
    i, j = m - 1, n - 1
    path.append((i, j))
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

    return {
        "score": matrix[m - 1][n - 1],
        "aligned_seq1": aligned1,
        "aligned_seq2": aligned2,
        "matrix": matrix,
        "path": list(reversed(path)),
    }


def smith_waterman(seq1: str, seq2: str, match: int = 1,
                    mismatch: int = -1, gap: int = -2) -> dict:
    """Local alignment. Finds the best-scoring shared subregion."""
    m, n = len(seq1) + 1, len(seq2) + 1
    matrix = np.zeros((m, n))
    traceback = np.zeros((m, n), dtype=int)  # 0=diag, 1=up, 2=left, 3=stop

    best_score, bi, bj = 0, 0, 0

    for i in range(1, m):
        for j in range(1, n):
            score = match if seq1[i - 1] == seq2[j - 1] else mismatch
            diagonal = matrix[i - 1][j - 1] + score
            up = matrix[i - 1][j] + gap
            left = matrix[i][j - 1] + gap
            best = max(diagonal, up, left, 0)
            matrix[i][j] = best

            if best == 0:
                traceback[i][j] = 3
            else:
                traceback[i][j] = [diagonal, up, left].index(best)

            if best > best_score:
                best_score, bi, bj = best, i, j

    aligned1, aligned2, path = "", "", []
    i, j = bi, bj
    path.append((i, j))
    while i > 0 and j > 0 and traceback[i][j] != 3:
        if traceback[i][j] == 0:
            aligned1 = seq1[i - 1] + aligned1
            aligned2 = seq2[j - 1] + aligned2
            i -= 1
            j -= 1
        elif traceback[i][j] == 1:
            aligned1 = seq1[i - 1] + aligned1
            aligned2 = "-" + aligned2
            i -= 1
        else:
            aligned1 = "-" + aligned1
            aligned2 = seq2[j - 1] + aligned2
            j -= 1
        path.append((i, j))

    return {
        "score": best_score,
        "aligned_seq1": aligned1,
        "aligned_seq2": aligned2,
        "start1": i, "end1": bi,
        "start2": j, "end2": bj,
        "matrix": matrix,
        "path": list(reversed(path)),
    }


# ---------------------------------------------------------------------
# Display helpers
# ---------------------------------------------------------------------
def print_alignment(aligned1: str, aligned2: str, score: float) -> None:
    print(f"\nAlignment score: {score}")
    print("  " + "".join(colorize_base(c) for c in aligned1))
    print("  " + "".join(colorize_base(c) for c in aligned2))


def print_matrix(matrix: np.ndarray, seq1: str, seq2: str, path: list) -> None:
    print("\nDP matrix (highlighted cells = traceback path):")
    path_set = set(path)
    col_labels = [" "] + list(seq2)

    header = "      " + "".join(f"{c:>5}" for c in col_labels)
    print(header)

    row_labels = [" "] + list(seq1)
    for i, row in enumerate(matrix):
        label = row_labels[i]
        line = f"  {label} "
        for j, val in enumerate(row):
            cell = f"{int(val):>4}"
            if (i, j) in path_set:
                line += f" {PATH_HIGHLIGHT}{cell}{RESET}"
            else:
                line += f" {cell}"
        print(line)


def run_demo(seq1: str, seq2: str, mode: str, match: int,
             mismatch: int, gap: int) -> None:
    if mode == "global":
        result = needleman_wunsch(seq1, seq2, match, mismatch, gap)
        print("\n=== Needleman-Wunsch (global alignment) ===")
    else:
        result = smith_waterman(seq1, seq2, match, mismatch, gap)
        print("\n=== Smith-Waterman (local alignment) ===")
        print(f"Matched region — seq1[{result['start1']}:{result['end1']}], "
              f"seq2[{result['start2']}:{result['end2']}]")

    print_alignment(result["aligned_seq1"], result["aligned_seq2"], result["score"])
    print_matrix(result["matrix"], seq1, seq2, result["path"])


# ---------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------
if __name__ == "__main__":
    default_seq1, default_seq2 = "GTTACATTT", "TGTTACA"

    seq1 = input(f"Sequence 1 [{default_seq1}]: ").strip().upper() or default_seq1
    seq2 = input(f"Sequence 2 [{default_seq2}]: ").strip().upper() or default_seq2

    mode_input = input("Algorithm — 'global' (NW) or 'local' (SW) [global]: ").strip().lower()
    mode = "local" if mode_input.startswith("l") else "global"

    match = int(input("Match score [1]: ") or 1)
    mismatch = int(input("Mismatch score [-1]: ") or -1)
    gap = int(input("Gap score [-2]: ") or -2)

    run_demo(seq1, seq2, mode, match, mismatch, gap)
