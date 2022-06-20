"""
Measure of heterogeneity of the regulatory outcome. Checks if genes that are always up-/down regulated always affected
by the same gene loss events?
"""

import numpy as np
from collections import OrderedDict
import networkx as nx


def update_gene_name(g, cds_dict, node_name):
    """
    This function get's the conversion dictionary and a node name and returns the original gene.
    @param cds_dict: the conversion dictionary
    @param node_name: the node name (such as NP..)
    @return: the original gene name
    """
    for gene, node in cds_dict.items():
        if node == node_name:
            g.node[node]['name'] = gene
            return gene


def net_het(g):
    """
    Calculate the network heterogeneity by finding the lost gene that affected each, how many times each lost gene
    affected in different strain and its distance from the node. Export the output to CSV file named
    "network_heterogeneity_1.csv".
    @param g: the network to estimates its heterogeneity.
    """
    for gene in g.nodes():
        lost_gene_effect = g.node[gene]['lost_genes_effect']
        lost_genes_counter = list(OrderedDict([(el, lost_gene_effect.count(el)) for el in lost_gene_effect]).items())
        down_reg, up_reg = sort_lost_genes_by_reg_type(g, gene, lost_genes_counter)
        g.nodes[gene]['down_reg'] = [item for sublist in down_reg for item in sublist]
        g.nodes[gene]['up_reg'] = [item for sublist in up_reg for item in sublist]
        up_reg_score = calc_het_value_per_gene(up_reg)
        down_reg_score = calc_het_value_per_gene(down_reg)
        main_reg_score, main_reg_type = define_main_reg_type_score(g, gene, up_reg_score, down_reg_score)
        g.node[gene]['het_score'] = main_reg_score
    return g


def define_main_reg_type_score(g, gene, up_reg_score, down_reg_score):
    """
    calculate the main regulation type by each node vector of weights and return the heterogeneity of this reg type.
    @param g: the network
    @param gene: the gene to calc its score
    @param up_reg_score:
    @param down_reg_score:
    @return:
    """
    main_reg_type = ""
    main_reg_score = 0
    if sum(g.node[gene]['weight_vec']) > 0:
        main_reg_type = "up"
        main_reg_score = up_reg_score
    elif sum(g.node[gene]['weight_vec']) < 0:
        main_reg_type = "down"
        main_reg_score = down_reg_score
    return main_reg_score, main_reg_type


def sort_lost_genes_by_reg_type(g, gene, lost_genes_counter):
    """
    Create two lists of up_reg and down_reg from the lost_gene_counter by the regulation type.
    The function checks how each lost gene affected the node by same calculation as the propagation calculation.
    @param g: the network
    @param gene: the gene to sort the lost gene that affected him
    @param lost_genes_counter: the lost gene and number of occurrences for each one of them.
    @return: the divided two list: up_reg and down_reg that contain the lost gene, number of occurrences and its
    distance from the node (distance = number of nodes in the path).
    """
    up_reg = []
    down_reg = []
    for i in range((len(lost_genes_counter))):
        path = nx.shortest_path(g, source=lost_genes_counter[i][0][0], target=gene)
        if lost_genes_counter[i][0][1] == -1:
            down_reg.append([lost_genes_counter[i][0][0], lost_genes_counter[i][1], len(path)-1])
        else:
            up_reg.append([lost_genes_counter[i][0][0], lost_genes_counter[i][1], len(path)-1])
    return down_reg, up_reg


def calc_het_value_per_gene(reg_type_arr):
    """
    score node reg type heterogeneity by the following process:
    1. norm_by_dist : number of occurrences per lost gene * distance from the node
    2. norm_by_occurrence : normalized the score of each node that the the minimal norm_by_dist value will be the
    maximal one (dividing each value by  the minimal value).
    3. Sum these values for each node (per reg type).
    @param reg_type_arr: array of arrays: [lost gene, num of occurrences, distance]
    @return: the calculated score
    """
    if len(reg_type_arr) > 0:
        norm_by_dist_arr = np.empty(len(reg_type_arr))
        for i in range(len(reg_type_arr)):
            val = (float(reg_type_arr[i][1]) / float(reg_type_arr[i][2]))
            norm_by_dist_arr[i] = (float(reg_type_arr[i][1]) / float(reg_type_arr[i][2]))
        norm_by_occurrence = norm_by_dist_arr / np.max(norm_by_dist_arr)
        return sum(norm_by_occurrence)
    else:
        return 0




