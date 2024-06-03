import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from needleman_w import needleman_wunsch
from smith_w import smith_waterman
import time

# Test sequences
fasta_file1 = "/Users/johnny/NW-vs-SW/sequence1.fasta"
fasta_file2 = "/Users/johnny/NW-vs-SW/sequence2.fasta"
sequence1 = list(SeqIO.parse(fasta_file1, "fasta"))
sequence2 = list(SeqIO.parse(fasta_file2, "fasta"))

seq1 = str(sequence1[0].seq)
seq2 = str(sequence2[0].seq)

# Measure execution time for Needleman-Wunsch
start_time = time.time()
needleman_wunsch(seq1, seq2)
nw_time = time.time() - start_time

# Measure execution time for Smith-Waterman
start_time = time.time()
smith_waterman(seq1, seq2)
sw_time = time.time() - start_time

# Plot the expected time complexity curves
plt.figure(figsize=(10, 6))
max_length = max(len(seq1), len(seq2))
lengths = np.linspace(100, max_length, 100)
expected_times = (lengths ** 2) / 1e6  # Adjust the scaling factor as needed

plt.plot(lengths, expected_times, label="Expected Time Complexity (O(nm))", color='green', linestyle='--')

# Plot the measured execution times
plt.scatter(len(seq1), nw_time, marker='o', color='red', label="Measured NW Time")
plt.scatter(len(seq2), sw_time, marker='o', color='blue', label="Measured SW Time")
plt.ylim(0, 200)
plt.yticks(np.arange(0, 201, 10))
plt.xlabel("Sequence Length")
plt.ylabel("Execution Time (seconds)")
plt.title("Time Complexity Comparison: Needleman-Wunsch vs. Smith-Waterman")
plt.legend()
plt.grid(True)
plt.show()

print("Needleman-Wunsch execution time:", nw_time)
print("Smith-Waterman execution time:", sw_time)
