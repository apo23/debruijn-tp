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

import random
from random import randint
import statistics
from operator import itemgetter
import argparse
import os
import sys
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt
random.seed(9001)

__author__ = "Apollinaire Roubert"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Apollinaire Roubert"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Apollinaire Roubert"
__email__ = "apollinaire.roubert@etu.u-paris.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file, "r") as filin:
        for line in filin:
            yield next(filin).strip()
            next(filin)
            next(filin)


def cut_kmer(sequence, kmer_size):
    for i in range(len(sequence)-(kmer_size-1)):
        yield sequence[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    sequences = read_fastq(fastq_file)
    kmer_dict = {}
    for seq in sequences:
        kmers = cut_kmer(seq, kmer_size)
        for k in kmers:
            if k in seq:
                if k not in kmer_dict:
                    kmer_dict[k] = 0
                kmer_dict[k] += 1
    return kmer_dict


def build_graph(kmer_dict):
    G = nx.DiGraph()
    for k in kmer_dict:
        prefix = k[:-1]
        suffix = k[1:]
        G.add_node(prefix)
        G.add_node(suffix)
        G.add_edge(prefix, suffix, weight = kmer_dict[k])
    #nx.draw(G, with_labels=True)
    #plt.show()
    return G


def get_starting_nodes(graph):
    nodes_in = []
    for node in graph.nodes:
        if not graph.in_edges(node):
            nodes_in.append(node)
    return nodes_in


def get_sink_nodes(graph):
    nodes_out = []
    for node in graph.nodes:
        if not graph.out_edges(node):
            nodes_out.append(node)
    return nodes_out


def get_contigs(graph, nodes_in, nodes_out):
    contigs = []
    paths_list = []
    for node_in in nodes_in:
        for node_out in nodes_out:
            paths_list.append(nx.all_simple_paths(graph, node_in, node_out))
    for paths in paths_list:
        for path in paths:
            contig = path[0]
            for i in range(1,len(path)):
                contig += path[i][-1]
            contigs.append( (contig, len(contig)) )
    return contigs


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(tuple_list, filout_name):
    with open(filout_name + ".fasta", "w") as filout:
        contig_counter = 0
        for contig_tuple in tuple_list:
            filout.write(">contig_{} len={}\n".format(contig_counter,
                                                      contig_tuple[1]))
            filout.write(fill(contig_tuple[0]) + "\n")
            contig_counter += 1


def std(values):
    pass


def path_average_weight(graph, path):
    pass


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    pass


def select_best_path(graph, path_list, path_lengths, path_mean_weight,
                     delete_entry_node=False, delete_sink_node=False):
    pass


def solve_bubble(graph, node_old, node_new):
    pass


def simplify_bubbles(graph):
    pass


def solve_entry_tips(graph, nodes_in):
    pass


def solve_out_tips(graph, nodes_out):
    pass


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    gen = read_fastq(args.fastq_file)
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmer_dict)
    nodes_in = get_starting_nodes(graph)
    nodes_out = get_sink_nodes(graph)
    contigs = get_contigs(graph, nodes_in, nodes_out)
    save_contigs(contigs, "test")


if __name__ == '__main__':
    main()
