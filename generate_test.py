
import random

bases = ["A", "T", "G", "C"]

seq1 = "".join(random.choices(bases, k=5000))

# Create a slightly different second sequence
seq2 = list(seq1)
seq2[250] = random.choice([b for b in bases if b != seq2[250]])
seq2 = "".join(seq2)

with open("test_sequences.fasta", "w") as f:
    f.write(f">Sequence_1\n{seq1}\n")
    f.write(f">Sequence_2\n{seq2}\n")

print("Test FASTA file created successfully!")
