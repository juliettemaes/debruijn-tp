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
from operator import itemgetter
import os
from pathlib import Path
import random
from random import randint
import statistics
import sys
import textwrap
from networkx import (
    DiGraph,
    draw_networkx_edges,
    draw_networkx_nodes,
    all_simple_paths,
    lowest_common_ancestor,
    has_path,
    draw,
    random_layout,
)
import matplotlib as plt

random.seed(9001)



from typing import Iterator, Dict, List

plt.use("Agg")

__author__ = "Juliette Maes"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Juliette Maes"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Juliette Maes"
__email__ = "juliette_maes@icloud.com"
__status__ = "Developpement"


def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments():  # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(
        description=__doc__, usage="{0} -h".format(sys.argv[0])
    )
    parser.add_argument(
        "-i", dest="fastq_file", type=isfile, required=True, help="Fastq file"
    )
    parser.add_argument(
        "-k", dest="kmer_size", type=int, default=22, help="k-mer size (default 22)"
    )
    parser.add_argument(
        "-o",
        dest="output_file",
        type=Path,
        default=Path(os.curdir + os.sep + "contigs.fasta"),
        help="Output contigs in fasta file (default contigs.fasta)",
    )
    parser.add_argument(
        "-f", dest="graphimg_file", type=Path, help="Save graph as an image (png)"
    )
    return parser.parse_args()



def read_fastq(fastq_file: Path) -> Iterator[str]:
    """Extract reads from fastq files.

    :param fastq_file: (Path) Path to the fastq file.
    :return: A generator object that iterate the read sequences.
    """
    with open(fastq_file, "r", encoding="utf-8") as f:
        while True:
            try:
                # Skip the identifier line (starting with '@')
                _ = next(f).strip()
                # Read the sequence line
                sequence = next(f).strip()
                # Skip the '+' line and the quality score line
                _ = next(f).strip()
                _ = next(f).strip()

                # Yield the sequence
                yield sequence
            except StopIteration:
                # Stop iteration when the file ends
                break



def cut_kmer(read: str, kmer_size: int) -> Iterator[str]:
    """Cut read into kmers of size kmer_size.

    :param read: (str) Sequence of a read.
    :return: A generator object that provides the kmers (str) of size kmer_size.
    """
    for i in range(len(read) - kmer_size + 1):
        yield read[i : i + kmer_size]


def build_kmer_dict(fastq_file: Path, kmer_size: int) -> Dict[str, int]:
    """Build a dictionnary object of all kmer occurrences in the fastq file

    :param fastq_file: (str) Path to the fastq file.
    :return: A dictionnary object that identify all kmer occurrences.
    """
    kmer_dict = {}
    for read in read_fastq(fastq_file):
        for kmer in cut_kmer(read, kmer_size):
            if kmer in kmer_dict:
                kmer_dict[kmer] += 1
            else:
                kmer_dict[kmer] = 1
    return kmer_dict


def build_graph(kmer_dict: Dict[str, int]) -> DiGraph:
    """Build the debruijn graph

    :param kmer_dict: A dictionnary object that identify all kmer occurrences.
    :return: A directed graph (nx) of all kmer substring and weight (occurrence).
    """
    graph = DiGraph()

    for kmer, weight in kmer_dict.items():
        # Get the prefix and suffix of the kmer
        prefix = kmer[:-1]
        suffix = kmer[1:]

        # Add an edge with the weight as the count of occurrences
        if graph.has_edge(prefix, suffix):
            # If the edge already exists, update its weight
            graph[prefix][suffix]['weight'] += weight
        else:
            # Otherwise, create a new edge with the initial weight
            graph.add_edge(prefix, suffix, weight=weight)

    return graph


def remove_paths(
    graph: DiGraph,
    path_list: List[List[str]],
    delete_entry_node: bool,
    delete_sink_node: bool,
) -> DiGraph:
    """Remove a list of path in a graph. A path is set of connected node in
    the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param delete_entry_node: (boolean) True->We remove the first node of a path
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    for path in path_list:
        if delete_entry_node and delete_sink_node:
            # Delete the entire path, including the first and last nodes
            nodes_to_remove = path
        elif delete_entry_node:
            # Delete all nodes except the last one
            nodes_to_remove = path[:-1]
        elif delete_sink_node:
            # Delete all nodes except the first one
            nodes_to_remove = path[1:]
        else:
            # Delete all nodes except the first and last ones
            nodes_to_remove = path[1:-1]

        # Remove the nodes from the graph
        graph.remove_nodes_from(nodes_to_remove)

    return graph


def select_best_path(
    graph: DiGraph,
    path_list: List[List[str]],
    path_length_list: List[int],
    weight_avg_list: List[float],
    delete_entry_node: bool = False,
    delete_sink_node: bool = False,
) -> DiGraph:
    """Select the best path between different paths

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param path_length_list: (list) A list of length of each path
    :param weight_avg_list: (list) A list of average weight of each path
    :param delete_entry_node: (boolean) True->We remove the first node of a path
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    # Compute average weight of each path
    sd_weight = statistics.stdev(weight_avg_list)
    # Select the best path
    if  sd_weight > 0:
        # Select the path with the highest average weight
        best_path_index = weight_avg_list.index(max(weight_avg_list))
    elif sd_weight == 0:
        # Compare the length of each path
        if statistics.stdev(path_length_list) > 0 or len(path_length_list) == 1: # If only one path
            # Select the longest path
            best_path_index = path_length_list.index(max(path_length_list))
        else:
            # Select a random path
            best_path_index = randint(0, len(path_length_list) - 1)
    # Remove all paths except the best one
    paths_to_remove = [path_list[i] for i in range(len(path_list)) if i != best_path_index]
    remove_paths(
        graph,
        paths_to_remove,
        delete_entry_node,
        delete_sink_node,
    )

    return graph


def path_average_weight(graph: DiGraph, path: List[str]) -> float:
    """Compute the weight of a path

    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean(
        [d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)]
    )


def solve_bubble(graph: DiGraph, ancestor_node: str, descendant_node: str) -> DiGraph:
    """Explore and solve bubble issue

    :param graph: (nx.DiGraph) A directed graph object
    :param ancestor_node: (str) An upstream node in the graph
    :param descendant_node: (str) A downstream node in the graph
    :return: (nx.DiGraph) A directed graph object
    """
    # Find all simple paths between the ancestor and descendant nodes
    paths = list(all_simple_paths(graph, ancestor_node, descendant_node))
    # Get the length of each path
    path_length_list = [len(path) for path in paths]
    # Get the average weight of each path
    weight_avg_list = [path_average_weight(graph, path) for path in paths]
    # Select the best path
    select_best_path(
        graph,
        paths,
        path_length_list,
        weight_avg_list,
        delete_entry_node=False,
        delete_sink_node=False,
    )
    return graph


def simplify_bubbles(graph: DiGraph) -> DiGraph:
    """Detect and explode bubbles

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    bubble = False
    for node in graph.nodes():
        # get predecessor list
        predecessors = list(graph.predecessors(node))
        if len(predecessors) > 1:
            # for each combination of predecessors
            for i in range(len(predecessors)):
                for j in range(i + 1, len(predecessors)):
                    ancestor_node = lowest_common_ancestor(graph, predecessors[i], predecessors[j])
                    if ancestor_node is not None:
                        bubble = True
                        break
        if bubble:
            break
    if bubble:
        simplify_bubbles(solve_bubble(graph, ancestor_node, node))
    return graph


def solve_entry_tips(graph: DiGraph, starting_nodes: List[str]) -> DiGraph:
    """Remove entry tips
    :param graph: (nx.DiGraph) A directed graph object
    :param starting_nodes: (list) A list of starting nodes
    :return: (nx.DiGraph) A directed graph object
    """
    while True:
        entry_tips = [node for node in graph if len(list(graph.predecessors(node))) > 1]
        if not entry_tips:
            break
        for node in entry_tips:
            # Find all paths from starting nodes to the entry tip
            paths = [
                list(
                    all_simple_paths(graph, node_start_i, node)
                    ) for node_start_i in starting_nodes
                    ]
            paths = [path[0] for path in paths if len(path) > 0] # Remove empty paths
            lengths = [len(path) - 1 for path in paths]
            # Compute the average weight of the path if the path length is greater than 1
            weights = [path_average_weight(graph, path) if lengths[i] > 1 else graph.get_edge_data(*path)["weight"]
                       for i, path in enumerate(paths)]
            # Remove all paths except the best one
            graph = select_best_path(graph, paths, lengths, weights,
                                     delete_entry_node=True, delete_sink_node=False)
    return graph


def solve_out_tips(graph: DiGraph, ending_nodes: List[str]) -> DiGraph:
    """Remove out tips
    :param graph: (nx.DiGraph) A directed graph object
    :param ending_nodes: (list) A list of ending nodes
    :return: (nx.DiGraph) A directed graph object
    """
    while True:
        found_tip = False
        for node in graph:
            node_success = list(graph.successors(node))
            if len(node_success) > 1:
                paths = [list(all_simple_paths(graph, node, node_end_i))\
                         for node_end_i in ending_nodes]
                paths = [path[0] for path in paths if len(path) > 0]
                lengths = [len(path) - 1 for path in paths]
                weights = [path_average_weight(graph, path) if lengths[i] > 1 else \
                           graph.get_edge_data(*path)["weight"]
                           for i, path in enumerate(paths)]
                graph = select_best_path(graph, paths, lengths, weights,
                                         delete_entry_node=False, delete_sink_node=True)
                found_tip = True
                break
        if not found_tip:
            break
    return graph



def get_starting_nodes(graph: DiGraph) -> List[str]:
    """Get nodes without predecessors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without predecessors
    """
    for node in graph.nodes():
        return [
            node for node in graph.nodes() if len(list(graph.predecessors(node))) == 0
            ]



def get_sink_nodes(graph: DiGraph) -> List[str]:
    """Get nodes without successors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without successors
    """
    for node in graph.nodes():
        return [
            node for node in graph.nodes() if len(list(graph.successors(node))) == 0
            ]


def get_contigs(
    graph: DiGraph, starting_nodes: List[str], ending_nodes: List[str]
) -> List:
    """Extract the contigs from the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param starting_nodes: (list) A list of nodes without predecessors
    :param ending_nodes: (list) A list of nodes without successors
    :return: (list) List of [contiguous sequence and their length]
    """
    list_contigs = []
    for start_node in starting_nodes:
        for end_node in ending_nodes:
            if has_path(graph, start_node, end_node):
                for path in all_simple_paths(graph, start_node, end_node):
                    contig = path[0]
                    for i in range(1, len(path)):
                        contig += path[i][-1]
                    list_contigs.append((contig, len(contig)))
    return list_contigs


def save_contigs(contigs_list: List[str], output_file: Path) -> None:
    """Write all contigs in fasta format

    :param contig_list: (list) List of [contiguous sequence and their length]
    :param output_file: (Path) Path to the output file
    """
    with open(output_file, 'w', encoding='utf-8') as f:
        for i, (contig, length) in enumerate(contigs_list):
            # Write the header
            f.write(f">contig_{i} len={length}\n")
            # Write the contig sequence wrapped to 80 characters per line
            f.write(textwrap.fill(contig, width=80) + "\n")


def draw_graph(graph: DiGraph, graphimg_file: Path) -> None:  # pragma: no cover
    """Draw the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param graphimg_file: (Path) Path to the output file
    """
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d["weight"] > 3]
    # print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d["weight"] <= 3]
    # print(elarge)
    # Draw the graph with networkx
    pos = random_layout(graph)
    draw_networkx_nodes(graph, pos, node_size=6)
    draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    draw_networkx_edges(
        graph, pos, edgelist=esmall, width=6, alpha=0.5, edge_color="b", style="dashed"
    )
    # nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file.resolve())


# ==============================================================
# Main program
# ==============================================================
def main() -> None:  # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Build kmer dictionnary
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)

    # Build the graph
    print("Building graph")
    graph = build_graph(kmer_dict)

    # Remove bubbles
    print("Simplifying bubbles")
    simplify_bubbles(graph)

    # Solve entry tips
    print("Solving entry tips")
    starting_nodes = get_starting_nodes(graph)
    graph = solve_entry_tips(graph, starting_nodes)

    # Solve out tips
    print("Solving out tips")
    ending_nodes = get_sink_nodes(graph)
    graph = solve_out_tips(graph, ending_nodes)

    # Write output contigs
    print("Writing contigs")
    contigs_list = get_contigs(graph, get_starting_nodes(graph), get_sink_nodes(graph))
    save_contigs(contigs_list, args.output_file)
    print("Contigs written in", args.output_file)

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)


if __name__ == "__main__":  # pragma: no cover
    main()
