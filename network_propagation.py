import networkx as nx
import numpy as np
import pandas as pd


def network_propagation(g, lost_genes, t1_exist_genes):
    """
    Prpagate the network- update nodes weight according to the lost genes
    @param g: the graph to propagate
    @param lost_genes: a dictionary of pairs and the genes lost in the progeny only
                       {(t1, t2): [lost genes between t1 to t2],...}
    @param exist_genes: a dictionary of pairs and the genes exist in the progenitor {(t1, t2): [exist genes in t1],...}
    @return: the updated graph
    """

    for strain_by_time in lost_genes:
        for gene in lost_genes[strain_by_time]:  # gene in lost_genes and in t1_exist
            if gene in g.nodes():
                g.node[gene]['lost_genes_per_patient'].append((strain_by_time, gene))  # fot the heterogeneity
                g.node[gene]['calc_weight'][0] = -1
                temp_node = gene
                visited = create_visited(g, lost_genes, strain_by_time, t1_exist_genes)
                queue = [temp_node]
                while queue:
                    temp_node = queue.pop(0)
                    for neighbor in g.neighbors(temp_node):
                        if neighbor in t1_exist_genes[strain_by_time]:
                            path = nx.shortest_path(g, source=gene, target=neighbor)
                            if len(path) == 1:
                                continue
                            reg_type = calc_node_weight(g, neighbor, path)
                            # fot the heterogeneity
                            if (strain_by_time, gene) not in g.node[neighbor]['lost_genes_per_patient']:
                                g.node[neighbor]['lost_genes_effect'].append((gene, reg_type))
                                g.node[neighbor]['lost_genes_per_patient'].append((strain_by_time, gene))
                            update_visited(neighbor, queue, visited)

        clear_nodes_calc_weight(g, t1_exist_genes[strain_by_time], strain_by_time, lost_genes[strain_by_time])
    return g


def clear_nodes_calc_weight(g, time_exist_genes, strain, lost_genes):
    """
    clear the node's calc weight parameter after each run
    @param g: the graph to clean
    @return: the clear graph
    """
    for gene in g.nodes():
        if gene in time_exist_genes:
            if gene not in lost_genes:
                g.node[gene]['weight_vec'].append(g.node[gene]['calc_weight'][0])
                g.node[gene]['weight_dict'][strain] = g.node[gene]['calc_weight'][0]
        g.node[gene]['calc_weight'] = [0]


def update_visited(node, queue, visited):
    """
    update the visited dict that a node has been visited and add it to the queue
    @param node: the node that has been visited
    @param queue: the queue to add to
    @param visited: the visited dict to update
    """
    if visited[node] == 0:
        queue.append(node)
        visited[node] = 1


def create_visited(g, lost_genes, tf, t1_exist_genes):
    """
    create visited nodes dict the mark the lost genes as visited
    @param g: the graph
    @param lost_genes: the lost gene to mark
    @param tf: the current TF
    @:param t1_exist_genes: set of genes that exist in t1
    @return: the visited dic
    """
    visited = {}
    for node in g.nodes():
        if node not in lost_genes[tf] and node in t1_exist_genes[tf]:
            visited[node] = 0
        else:
            visited[node] = 1  # the gene doesn't exist in the network (in lost genes or not in exist)
    return visited


def calc_node_weight(g, node, path):
    """
    calc node weight by its path
    @param g: the graph
    @param node: the node to calc its weight
    @param path: the nodes path from the root
    """
    sign = 1
    for path_index in range(len(path) - 1):
        edge_sign = g[path[path_index]][path[path_index + 1]]['weight']
        sign *= edge_sign
    g.node[node]['calc_weight'].append((len(path) - 1, sign))
    updated_weight = 0
    for tup in g.node[node]['calc_weight'][1:]:
        updated_weight += (1 / float(tup[0])) * (-tup[1])
    g.node[node]['calc_weight'][0] = val_sign(updated_weight)
    return val_sign(updated_weight)


def val_sign(value):
    """
    return -1, 1, or 0 by the value sign
    @param value: the value to asses
    @return: -1 if neg, 1 if pos and 0 otherwise
    """
    if value > 0:
        return 1
    elif value < 0:
        return -1
    else:
        return 0
