import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
import random
import re

from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

# testing
# filename = 'example_net.txt'
# graph = read_graph(filename)
# seq1 = gen_rand_strand(1,15)
# seq2 = gen_rand_strand(1,10)
# seqs = [seq1,seq2]

def gen_rand_strand(num_strands,len_strand):

    bases = ['A', 'T', 'C', 'G']

    strand_list=[]
    for i in range(num_strands):
        strand = ''.join(random.choice(bases) for _ in range(len_strand))
        strand_list.append(strand)

    print("Generated DNA sequence(s) are in FASTA format (3':__:5').")
    seq_heatmap(strand_list)

    return strand_list

def load_fasta_strand(fasta_path):
    
    strand_list=[]
    with open(str(fasta_path)) as fasta_file:
        for line in fasta_file:
            if "32630" in line:
                strand_list.append(str(next(fasta_file)).strip())
    
    print("Extracted DNA sequences from FASTA file (3':__:5').")
    seq_heatmap(strand_list)
    
    return strand_list

def seq_complement(single_strand,im_flag=True):
    
    if type(single_strand) != str:
        print("ERROR: Strand must be a objecttype:string, not objecttype:",str(type(single_strand)))
        return
    
    map_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    
    strand_complement = re.compile("|".join(map_dict.keys())).sub(lambda ele: map_dict[re.escape(ele.group(0))], single_strand.upper())
    
    for char in strand_complement:
        if (char not in map_dict) == True:
            print("ERROR: Strand contains non-DNA nucleotides; please input only DNA nucleotide sequences")
            return
        
    if im_flag == True:
        seq_heatmap(strand_complement,True)
        print("NOTE: Generated complement strand is in FASTA (5':__:3') format")
    
    return strand_complement

def seq_heatmap(seq, comp_flag=False):
    
    nucleotide_mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    
    if type(seq) == str or len(seq) == 1:
        
        if len(seq) == 1:
            seq = seq[0]
        
        numerical_seq = [nucleotide_mapping[base] for base in seq]
        dna_array = np.array(numerical_seq).reshape(1, -1)

        cmap = sns.color_palette(['#db5f57', '#d3db57', '#57db5f', '#57d3db'])

        plt.figure(figsize=(len(seq), 1))
        ax = sns.heatmap(dna_array, cmap=cmap, cbar=False, annot=False, linewidths=0.5)
        plt.xticks([])
        
        if comp_flag == True:
            plt.ylabel('Sequence #\n(5\':__:3\')')
        else:
            plt.ylabel('Sequence #\n(3\':__:5\')')
        
        for i in range(len(seq)):
            ax.text(i + 0.5, 0.5, seq[i], ha='center', va='center', color='Black')
            
        plt.show()
        
    else:
        
        max_len = len(max(seq, key=len))
        
        numerical_seq = []
        for i in range(len(seq)):
            if len(seq[i]) < max_len:
                numerical_seq_pad = [nucleotide_mapping[base] for base in seq[i]]
                numerical_seq_pad += [-1] * (max_len - len(numerical_seq_pad))
                numerical_seq.append(numerical_seq_pad)
            else:
                numerical_seq.append([nucleotide_mapping[base] for base in seq[i]])
        
        dna_array = np.array(numerical_seq)

        it = iter(seq)
        the_len = len(next(it))
        if not all(len(l) == the_len for l in it):
            cmap = sns.color_palette(['white','#db5f57', '#d3db57', '#57db5f', '#57d3db'])
        else:
            cmap = sns.color_palette(['#db5f57', '#d3db57', '#57db5f', '#57d3db'])

        plt.figure(figsize=(max_len, len(numerical_seq)))
        ax = sns.heatmap(dna_array, cmap=cmap, cbar=False, annot=False, linewidths=0.5)
        plt.xticks([])

        if comp_flag == True:
            plt.ylabel('Sequence #\n(5\':__:3\')')
        else:
            plt.ylabel('Sequence #\n(3\':__:5\')')
        
        for i in range(max_len):
            for j in range(len(numerical_seq)):
                if dna_array[j, i] != -1:
                    ax.text(i + 0.5, j + 0.5, seq[j][i] if j == 0 else seq[j][i], ha='center', va='center', color='Black')

        plt.show()

def seq_alignment(seq1, seq2):

    print("Assuming that BOTH sequence 1 and 2 are in FASTA (3':__:5') format.")
    sequence1 = Seq(seq1)
    sequence2 = Seq(seq2[::-1])

    print("Performing alignment...")
    alignments = pairwise2.align.globalxx(sequence1, sequence2) 

    best_alignment = alignments[0]
    print("\nBest Alignment:\n")
    print(format_alignment(*best_alignment))

    aligned_seq1 = best_alignment.seqA
    aligned_seq2 = best_alignment.seqB
    matches = sum(base1 == base2 for base1, base2 in zip(aligned_seq1, aligned_seq2))
    similarity_global = (matches / len(aligned_seq1)) * 100
    similarity_local =  (matches / max(len(seq1), len(seq2))) * 100
    
    print(f"Similarity (Global): {similarity_global:.2f}%")
    print(f"Similarity (Local): {similarity_local:.2f}%")
    print(f"Alignment Score: {best_alignment.score}")
    
    nucleotide_mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3,'-': 4}
    max_len = len(alignments[0][0])
    seq = [alignments[0][0],alignments[0][1]]
    
    numerical_seq = []
    for i in range(len(seq)):
        numerical_seq.append([nucleotide_mapping[base] for base in seq[i]])
    
    dna_array = np.array(numerical_seq)


    cmap = sns.color_palette(['#db5f57', '#d3db57', '#57db5f', '#57d3db','white'])


    plt.figure(figsize=(max_len, len(numerical_seq)))
    ax = sns.heatmap(dna_array, cmap=cmap, cbar=False, annot=False, linewidths=0.5)
    plt.xticks([])
    plt.yticks([0.5, 1.5], ['1       \n((3\':__:5\'))', '2       \n((5\':__:3\'))'])
    plt.ylabel('Sequence #')
    
    for i in range(max_len):
        for j in range(len(numerical_seq)):
            if dna_array[j, i] != -1:
                ax.text(i + 0.5, j + 0.5, seq[j][i] if j == 0 else seq[j][i], ha='center', va='center', color='Black')

    plt.show()

def read_graph(filename):

    G = nx.DiGraph()

    with open(filename, 'r') as f:
        for line in f:
            node1, node2 = line.strip().split()
            G.add_edge(node1, node2)

    nx.draw(G, with_labels=True, node_color='red', node_size=500, font_size=10, font_color='black', edge_color='black', arrowsize=10)
    plt.show()

    return G

def hamiltonian_check(graph):
    
    def is_safe(v, pos, path):

        if not graph.has_edge(path[pos - 1], v):
            return False
        if v in path:
            return False
        return True

    def hamiltonian_path_util(path, pos):

        if pos == len(graph.nodes):
            return True

        for v in graph.nodes:
            if is_safe(v, pos, path):
                path[pos] = v
                if hamiltonian_path_util(path, pos + 1):
                    return True
                path[pos] = -1

        return False

    path = [-1] * len(graph.nodes)
    path[0] = list(graph.nodes)[0]  # Start from the first node

    if not hamiltonian_path_util(path, 1):
        return False

    return True

def create_linked_list(nodes):
  linked_list = []
  for i in range(len(nodes) - 1):
    linked_list.append((nodes[i], nodes[i+1]))
  return linked_list

def adleman_sim(graph):

    # if hamiltonian_check(graph) == False:
    #     print("Error: directed graph does not have a hamiltonian path")

    edges = list(graph.edges())
    nodes = list(graph.nodes())
    nodes.sort()

    V = len(nodes)

    pairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    k = 10 # change at will

    node_strands = [''.join([ random.choice(list(pairs.keys())) for _ in range(k) ]) for _ in range(V)]
    print(f"Vertex Strands Generated: {node_strands}")

    edge_strands = [ node_strands[int(v1)][(k//2):] + node_strands[int(v2)][:(k//2)] for (v1, v2) in edges ]
    edge_strands = [ strand.replace(node_strands[0][-(k//2):], node_strands[0]).replace(node_strands[-1][:(k//2)], node_strands[-1]) for strand in edge_strands]

    print(f"\nEdge Strands Generated: {edge_strands}")

    complements = [ seq_complement(v,False) for v in node_strands ]
    print(f"\nComplements Generated: {complements}")

    all_edges = edge_strands * 500
    random.shuffle(all_edges)

    path_strands = []
    growing_strand = ""

    for idx, s in enumerate(all_edges):
        if len(growing_strand) == 0:
            growing_strand += s
        else:
            target_comp = seq_complement(growing_strand[-(k//2):],False) + seq_complement(s[:(k//2)],False)
            if target_comp in complements:
                growing_strand += s

        if growing_strand[-k:] == node_strands[-1]:
            path_strands.append(growing_strand)
            growing_strand = ""

    if len(path_strands) == 0:
        print("Warning: randomization of paths did not yield anything (potentially due to the number of edges).\nIncrease edge list multiplier, if possible.")
    else:
        print(f"\nStrands Created: {len(path_strands)} ")

    in_and_out_strands = [ path for path in path_strands if path[:k] == node_strands[0] and path[-k:] == node_strands[-1] ]
    print(f"\nStrands Starting at {nodes[0]} and Ending at {nodes[V-1]}: {len(in_and_out_strands)}")

    n_step_paths = [ path for path in in_and_out_strands if len(path) == V * 10 ]
    print(f"\nStrands with {V} Steps: {len(n_step_paths)}")
    
    included = n_step_paths
    for i in range(len(node_strands)):
        included = [ path for path in included if node_strands[i] in path ]
        if i == 0 or i == len(node_strands) - 1:
            print(f"Already Checked for {nodes[i]}")
        else:
            print(f"Eliminating Paths Not Including {nodes[i]}: {len(included)} Remaining")
    print(f"Strands Including All Vertices Once: {len(included)}")

    solution = included[0]
    path = []
    for i in range(0, len(solution), 10):
        vertex_strand = solution[i:i+10]
        vertex_num = node_strands.index(vertex_strand)
        decoded_letters = nodes[vertex_num]
        path.append(decoded_letters)

    print(f"Solution: {str(path)}")
    
    seq_heatmap(included[0],True)
    
    edge_list_to_color = create_linked_list(path)
    
    edge_colors = ['grey'] * len(graph.edges())
    for i, edge in enumerate(graph.edges()):
        if edge in edge_list_to_color:
            edge_colors[i] = 'red'  # Adjust color as needed

    nx.draw(graph, node_color='red', node_size=500, with_labels=True, edge_color=edge_colors, arrowsize=10)
    plt.show()