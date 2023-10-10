#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import textwrap
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path): # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file
    
    :raises ArgumentTypeError: If file doesn't exist
    
    :return: (str) Path 
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        defaut=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file (default contigs.fasta)")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as an image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    """Extract reads from fastq files.

    :param fastq_file: (str) Path to the fastq file.
    :return: A generator object that iterate the read sequences. 
    """
    with open(fastq_file, "rt") as reader :
        for line in reader :
            yield next(reader)
            next(reader)
            next(reader)


def cut_kmer(read, kmer_size):
    """Cut read into kmers of size kmer_size.
    
    :param read: (str) Sequence of a read.
    :return: A generator object that iterate the kmers of of size kmer_size.
    """
    for i in range(0, len(read)-kmer_size+1):
        yield read [i, i + kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    """Build a dictionnary object of all kmer occurrences in the fastq file

    :param fastq_file: (str) Path to the fastq file.
    :return: A dictionnary object that identify all kmer occurrences.
    """
    kmer_dict = {}  
    with open(fastq_file, 'r') as file:
        kmer_dict[kmer] = kmer_dict.get(kmer, 0) + 1  
    return kmer_dict
    
kmer_occurrences = build_kmer_dict(fastq_file, kmer_size)
print(kmer_occurrences)
    

def build_graph(kmer_dict):
    """Build the debruijn graph

    :param kmer_dict: A dictionnary object that identify all kmer occurrences.
    :return: A directed graph (nx) of all kmer substring and weight (occurrence).
    """
    debruijn_graph = nx.DiGraph()

    for kmer, count in kmer_dict.items():
        prefix = kmer[:-1]  
        suffix = kmer[1:]  
        debruijn_graph.add_edge(prefix, suffix, weight=count)
    return debruijn_graph
    
    kmer_dict = {'ATG': 2, 'TGC': 1, 'GCA': 3, 'CAT': 2}
    debruijn_graph = build_graph(kmer_dict)


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """Remove a list of path in a graph. A path is set of connected node in
    the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    modified_graph = graph.copy()
    for path in path_list:
        if delete_entry_node and path[0] in modified_graph:
            modified_graph.remove_node(path[0])
        if delete_sink_node and path[-1] in modified_graph:
            modified_graph.remove_node(path[-1])
    return modified_graph


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    """Select the best path between different paths

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param path_length_list: (list) A list of length of each path
    :param weight_avg_list: (list) A list of average weight of each path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    scores = [path_length_list[i] + weight_avg_list[i] for i in range(len(path_list))]
    best_path_index = scores.index(min(scores))
    best_path = path_list[best_path_index]  
    modified_graph = graph.copy()
    if delete_entry_node and best_path[0] in modified_graph:
        modified_graph.remove_node(best_path[0])
    if delete_sink_node and best_path[-1] in modified_graph:
        modified_graph.remove_node(best_path[-1])
    return modified_graph


def path_average_weight(graph, path):
    """Compute the weight of a path

    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])
    debruijn_graph = nx.DiGraph()
    path = ['A', 'B', 'C', 'D']
    average_weight = path_average_weight(debruijn_graph, path)
    print("Average Weight of the Path:", average_weight)


def solve_bubble(graph, ancestor_node, descendant_node):
    """Explore and solve bubble issue

    :param graph: (nx.DiGraph) A directed graph object
    :param ancestor_node: (str) An upstream node in the graph 
    :param descendant_node: (str) A downstream node in the graph
    :return: (nx.DiGraph) A directed graph object
    """
    paths = nx.all_simple_paths(graph, source=ancestor_node, target=descendant_node)
    if not paths:
        return graph

    shortest_path = min(paths, key=len)

    paths_to_remove = [path for path in paths if path != shortest_path]
    modified_graph = remove_paths(graph, paths_to_remove, delete_entry_node=True, delete_sink_node=True)
    return modified_graph


def simplify_bubbles(graph):
    """Detect and explode bubbles

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    simplified_graph = graph.copy()  
    for node in graph.nodes():
        successors = list(graph.successors(node))
        if len(successors) == 2:
            ancestor_node, descendant_node = successors[0], successors[1]
            if len(graph.successors(ancestor_node)) == 1 and len(graph.successors(descendant_node)) == 1:
                simplified_graph = solve_bubble(simplified_graph, ancestor_node, descendant_node)
    return simplified_graph


def solve_entry_tips(graph, starting_nodes):
    """Remove entry tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    modified_graph = graph.copy()  
    entry_tips = [node for node in graph.nodes() if modified_graph.in_degree(node) == 1 and node not in starting_nodes]

    for tip in entry_tips:
        predecessors = list(modified_graph.predecessors(tip))
        if tip in modified_graph and len(predecessors) == 1:
            modified_graph.remove_node(tip) 
    return modified_graph
    

def solve_out_tips(graph, ending_nodes):
    """Remove out tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    modified_graph = graph.copy()  
    out_tips = [node for node in graph.nodes() if modified_graph.out_degree(node) == 1 and node not in ending_nodes]

    for tip in out_tips:
        successors = list(modified_graph.successors(tip))
        if tip in modified_graph and len(successors) == 1:
            modified_graph.remove_node(tip)  
    return modified_graph


def get_starting_nodes(graph):
    """Get nodes without predecessors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without predecessors
    """
    starting_nodes = [node for node in nx.topological_sort(graph) if graph.in_degree(node) == 0]
    return starting_nodes


def get_sink_nodes(graph):
    """Get nodes without successors in a directed graph.

    :param graph: (nx.DiGraph) A directed graph object.
    :return: (list) A list of all nodes without successors.
    """
    sink_nodes = [node for node in graph.nodes() if graph.out_degree(node) == 0]
    return sink_nodes


def get_contigs(graph, starting_nodes, ending_nodes):
    """Extract the contigs from the graph.

    :param graph: (nx.DiGraph) A directed graph object.
    :param starting_nodes: (list) A list of nodes without predecessors.
    :param ending_nodes: (list) A list of nodes without successors.
    :return: (list) List of [contiguous sequence and their length].
    """
    contigs = []
    for start_node in starting_nodes:
        for end_node in ending_nodes:
            paths = list(nx.all_simple_paths(graph, source=start_node, target=end_node))
            for path in paths:
                contig_sequence = ''.join(graph.nodes[node]['sequence'] for node in path)
                contig_length = sum(graph.nodes[node]['length'] for node in path)
                contigs.append([contig_sequence, contig_length])
    return contigs


def save_contigs(contigs_list, output_file):
    """Write all contigs in fasta format.

    :param contigs_list: (list) List of [contiguous sequence and their length].
    :param output_file: (str) Path to the output file.
    """
    with open(output_file, 'w') as f:
        for idx, (sequence, length) in enumerate(contigs_list, start=1):
            f.write(f'>Contig_{idx}_Length_{length}\n{sequence}\n')
            

def draw_graph(graph, graphimg_file): 
    """Draw the graph and save it as an image.

    :param graph: (nx.DiGraph) A directed graph object.
    :param graphimg_file: (str) Path to the output image file.
    """
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmer_dict)
    
    # éliminant les "tips"
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    graph = simplify_bubbles(graph)
    graph = solve_entry_tips(graph, starting_nodes)
    graph = solve_out_tips(graph, ending_nodes)
    
    contigs = get_contigs(graph, starting_nodes, ending_nodes)
    save_contigs(contigs, args.output_file)
 
    if args.graphimg_file:
        draw_graph(graph, args.graphimg_file)


if __name__ == '__main__': # pragma: no cover
    main()