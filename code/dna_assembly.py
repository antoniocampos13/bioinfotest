import random
from collections import Counter, deque
from itertools import chain, combinations
from pathlib import Path

def get_kmers(sequences, k=None):
    """Generate list of k-mers found in a list of sequences

    Args:
        sequences (list): Sequencing reads
        k (int, optional): Size of k-mer. Defaults to None. If None k will equal the half of the read (rounded if read length is odd)

    Returns:
        dictionary: k-mer dictionary with counts of how many times it appeared on the reads
    """    

    if k:
        k = k
    else:
        k = len(min(sequences, key=len)) // 2

    res = []
    for seq in sequences:
        res.append(
            [
                seq[x:y]
                for x, y in combinations(range(len(seq) + 1), r=2)
                if len(seq[x:y]) == k
            ]
        )

    kmers = list(chain.from_iterable(res))

    return dict(Counter(kmers))


def get_edges(kmers):
    """Determine De Bruijn's edges from a dictionary of k-mers

    Args:
        kmers ([dict]): Dictionary of k-mers, output of get_kmers function

    Returns:
        [set]: A set containing tuples representing the edges of a De Bruijn's graph
    """    

    combs = combinations(kmers, 2)

    edges = set()

    for km1, km2 in combs:
        if km1[1:] == km2[:-1]:
            edges.add((km1[:-1], km2[:-1]))
        if km1[:-1] == km2[1:]:
            edges.add((km2[:-1], km1[:-1]))

    return edges


def get_superstring(edges):
    """Generates a contig based on edges of a De Bruijn's graph via Eulerian path walk

    Args:
        edges ([set]): A set containing tuples representing the edges of a De Bruijn's graph. Output of get_edges function

    Returns:
        str: The shortest superstring that contains all other strings represented in the edges input
    """
    deck = deque()

    node1, node2 = random.choice(tuple(edges))
    deck = deque([node1, node2])

    while len(deck) <= len(edges):

        for node1, node2 in edges:

            if node1 == deck[-1]:
                deck.append(node2)

            elif node2 == deck[0]:
                deck.appendleft(node1)

            else:
                continue

    path = list(deck)

    superstring = path[0] + "".join(map(lambda x: x[-1], path[1:]))
    return superstring


def assemble(input_file_path, output_file_path="output.txt", k=None):
    """Wrapper function for get_kmers, get_edges and get_superstring

    Args:
        file ([file]): Text file containing one sequencing read per line
        k (int, optional): Size of k-mer. Defaults to None. If None k will equal the half of the read (rounded if read length is odd)

    Returns:
        str: The shortest superstring that contains all other strings represented in the edges input
    """

    sequences = open_input(input_file_path)

    km = get_kmers(sequences, k=k)

    edges = get_edges(km)

    superstring = get_superstring(edges)

    export_output(superstring, output_file_path)

def get_genome(genome_length, seed=123):
    """Simulate a genome

    Args:
        genome_length (int, optional): Genome length.
        seed (int, optional): Seed for the random number generator. Defaults to 123.

    Returns:
        str: A simulated genome of desired length
    """
    random.seed(seed)

    simgenome = "".join(random.choices("ATCG", k=genome_length))

    return simgenome


def get_read(genome, read_length):
    """Simulate a read of desired length from a specified genomic sequence

    Args:
        genome (str): A simulated genome
        read_length (int): The length of each read (substring)

    Returns:
        list: A string containing a simulated sequencing read
    """
    idx = random.randrange(0, len(genome) - read_length + 1)

    return genome[idx : (idx + read_length)]


def simulate_reads(genome, genome_length, read_length, coverage, seed=123):
    """Simulate sequencing reads by taking substrings of a simulated genome

    Args:
        genome (str): A simulated genome
        genome_length (int): The length of the simulated genome
        read_length (int): The length of each read (substring).
        coverage (int, optional): The mean number of times a nucleotide was sequenced.
        seed (int, optional): Seed for the random number generator. Defaults to 123.

    Returns:
        list: A list containing simulated sequencing reads
    """
    random.seed(seed)

    n_reads =  ( coverage * genome_length ) // read_length

    simreads = [get_read(genome, read_length) for _ in range(n_reads)]

    return simreads

def open_input(file_path):
    """Function to open file with input (reads)

    Args:
        file_path (str): A string representing the path of the input file

    Returns:
        list: A list of reads to be handled by the other functions
    """
    wd = Path(file_path).resolve()

    with open(wd, "r") as file_in:
        sequences = [line.strip() for line in file_in]
    
    return sequences

def export_output(superstring, output_file_path):
    """Saves the assemble to disk

    Args:
        superstring (str): output of get_superstring function
        output_file_path (str): A string representing the path of the output file
    """
    wd = Path(output_file_path).resolve()

    with open(wd, "w") as file_out:
        file_out.write(superstring)