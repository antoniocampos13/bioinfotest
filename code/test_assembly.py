from dna_assembly import *
from Bio import pairwise2

genome = get_genome(genome_length=1000, seed=123)
reads = simulate_reads(genome, genome_length=1000, read_length=100, coverage=30, seed=123)

with open("input_teste.txt", "w") as f:
    for element in reads:
        f.write(element + "\n")

assemble("input_teste.txt","assembly.txt", k=55)

with open("assembly.txt", "r") as f:
    assembly = f.readline().strip()

pairwise2.align.globalxx(genome, assembly)