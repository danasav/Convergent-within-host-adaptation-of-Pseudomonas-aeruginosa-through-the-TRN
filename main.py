import csv
import matplotlib as mpl
import matplotlib.pyplot as plt
import convert_network_to_refseq
from get_raw_data import get_proper_samples, clean_s1, calc_time_intervals, get_network_df
from network_heterogeneity import net_het
from network_propagation import *
from config import *

NUMBER_OF_STRAINS_BY_TIME = 260

REGULATION_ATTRS = {"-": [-1, DOWN_REG_COLOR],
                    "+": [1, UP_REG_COLOR],
                    "?": [0, 'white'],
                    "d": [0, 'white'],
                    ".": [0, 'white']}


def create_graph(TF_dict, lost_genes=[], add_unknown_edges=False):
    """
    This function creates the network based on the TF dictionary.
    Each TF has an edge towards it's targets.
    Each node has the following attributes:
    - calc_weight - an array used to calc weight in each iteration (for each patient)
    - weight_vec - a vector of the calculated weights of all the iterations.
    - color - should change based of the node weight.
    - weight - the final weight (the mean of the weight_vec)
    @param add_unknown_edges:
    @param lost_genes: list of nodes that we don't want to add to the graph.
    @param TF_dict: a dictionary of {TF: [(target, regulation type), (other target, regulation type)...], .. }
    @return: A directed network with the attributes mentioned above.
    @rtype: networkx directed graph.
    """
    G = nx.DiGraph()
    for TF in TF_dict:
        if TF in lost_genes:
            continue
        for target, reg_type in TF_dict[TF]:
            if target in lost_genes:
                continue
            if REGULATION_ATTRS[reg_type][0] == 0 and not add_unknown_edges:
                continue
            else:
                G.add_node(TF, calc_weight=[0], weight_vec=[], color='blue', lost_genes_effect=[],
                           lost_genes_per_patient=[], weight_dict={}, weight_without_lost_vec=[])
                G.add_node(target, calc_weight=[0], weight_vec=[], color='blue', lost_genes_effect=[],
                           lost_genes_per_patient=[], weight_dict={}, weight_without_lost_vec=[])
                G.add_edge(TF, target, weight=REGULATION_ATTRS[reg_type][0], color=REGULATION_ATTRS[reg_type][1])
    return G


def data_parsing():
    """
    create tf factors dictionary from the raw data
    :return: TF factors dict- each gene has it tf factors and its regulation type (pos or neg)
    """
    data = pd.read_csv(NETWORK_FILE, sep='\t')
    tf_dict = {}
    for i, row in data.iterrows():
        # tf = row['Source']
        # gene = row['Target']
        tf = row['Source_Locus_id']
        gene = row['Target_Locus_id']
        interaction = row['Interaction']
        if interaction == 'Unknown':
            interaction = '?'
        if interaction == '+, +':
            interaction = '+'
        if tf in tf_dict:
            tf_dict[tf].append((gene, interaction))
        else:
            tf_dict[tf] = [(gene, interaction)]
    return tf_dict


def draw_graph(G, max_weight, min_weight, labels={}):
    """
    This function draws the network for each strain, with a color map from blue to red.
    @param G: The strain's network after calculations
    """
    fig = plt.figure(figsize=(19.20, 10.80))
    # fig = plt.figure(figsize=(50, 25))
    G.add_node('toPresent1', excluded_weight=abs(min_weight), weight=abs(min_weight))
    G.add_node('toPresent2', excluded_weight=-max_weight, weight=-max_weight)

    node_list = [node for node in G.nodes() if G.node[node]['weight'] == G.node[node]['excluded_weight']]
    node_colors = [G.node[node]['excluded_weight'] for node in node_list]

    pos = nx.nx_agraph.graphviz_layout(G)
    NODE_SIZE = 20
    nc = nx.draw_networkx_nodes(G, pos=pos, node_size=NODE_SIZE, node_color=node_colors, with_labels=False,
                                cmap=mpl.cm.coolwarm, nodelist=node_list)

    excluded_nodes = [node for node in G.nodes if G.node[node]['weight'] != G.node[node]['excluded_weight']]
    excluded_colors = [G.node[node]['excluded_weight'] for node in excluded_nodes]
    nx.draw_networkx_nodes(G, pos=pos, node_size=NODE_SIZE, with_labels=False,
                           cmap=mpl.cm.coolwarm, nodelist=excluded_nodes, node_color=excluded_colors)

    pos_regulation_edges = [(u, v) for u, v in G.edges() if G[u][v]['color'] == UP_REG_COLOR]
    neg_regulation_edges = [(u, v) for u, v in G.edges() if G[u][v]['color'] == DOWN_REG_COLOR]
    nx.draw_networkx_edges(G, pos=pos, node_size=NODE_SIZE, edgelist=pos_regulation_edges, arrowsize=11,
                           edge_color='darkgray')
    nx.draw_networkx_edges(G, pos=pos, node_size=NODE_SIZE, edgelist=neg_regulation_edges, arrowstyle='-[',
                           arrowsize=7, edge_color='darkgray')

    # handling the labels:
    vmin = min(node_colors)
    vmax = max(node_colors)
    jump = (vmax + abs(vmin)) / 6
    bar_labels = list(np.arange(vmin, vmax + jump, jump))
    cbar = plt.colorbar(nc, format='%.2f', ticks=bar_labels, shrink=0.75)
    plt.axis('off')
    fig.savefig('Network.png', dpi=600)
    plt.show()
    plt.clost()


def calculate_weight_frac(G, node):
    """
    Calculate the weight as a fraction (percentage) instead of the mean calculation 
    @param G: Our regulatory network
    @param node: the current node to calculate it's weight
    @return: the new node's weight
    """
    weight_sum = sum(G.node[node]['weight_vec'])
    num_of_genes = len(G.node[node]['weight_vec'])
    if weight_sum == 0:
        weight = 0
    elif weight_sum > 0:
        sum_up_reg = sum([weight for weight in G.node[node]['weight_vec'] if weight == 1])
        weight = (sum_up_reg / float(num_of_genes)) * 100
    else:
        sum_down_reg = sum([weight for weight in G.node[node]['weight_vec'] if weight == -1])
        weight = (sum_down_reg / float(num_of_genes)) * 100
    G.node[node]['weight_frac'] = weight


def add_plt_bar(cmap, vmax):
    """
    Creates a color map bar to the network plot
    @param cmap: the color map
    @param vmax: the maximum node weight
    """
    vmin = -vmax
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    jump = (vmax + abs(vmin)) / 4
    bar_labels = list(np.arange(vmin, vmax + jump, jump))
    plt.colorbar(sm, pad=0.04, shrink=0.5, format='%.3f', ticks=bar_labels)


def export_to_csv(G, sumed_gene_loss, strain=STRAIN_P):
    """
    Export the network data of each strain into a csv file.
    @param G: the directed graph network with all the weight loss data
    """
    csv_columns = ['Node name', 'Gene name', 'Locus tag', 'Weight', 'Excluded Weight', 'OutDegree', 'NumberOfLosses',
                   'Weight as Fraction', 'InDegree', 'Heterogeneity', 'Weight Dict', 'Weight Vec Length', 'ROS=-1',
                   'ROS=+1', 'ROS=0']
    csv_file = SUPPLEMENTARY_TABLE_1
    try:
        with open(csv_file, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
            writer.writeheader()
            for node in sorted(G.nodes()):
                data = {"Node name": node, "Gene name": G.node[node]['name'], "Locus tag": G.node[node]['locus_tag'],
                        'Weight': G.node[node]['weight'],
                        "Excluded Weight": G.node[node]['excluded_weight'],
                        'OutDegree': G.out_degree(node),
                        'NumberOfLosses': sumed_gene_loss[node], 'Weight as Fraction': G.node[node]['weight_frac'],
                        'InDegree': G.in_degree(node), 'Heterogeneity': G.node[node]['het_score'],
                        'Weight Dict': G.node[node]['weight_dict'],
                        'Weight Vec Length': len(G.node[node]['weight_vec'])}
                writer.writerow(data)
    except IOError:
        print("I/O error")


def sum_lost_genes(g, lost_genes):
    """
    calculates the number of strains in which each node #losses and the out degree in each network
    @param g: our TF network
    @param lost_genes: a dictionary of gene and the number of times it #losses (all "strains" combined) 
    @return: the number of times this gene #losses in all the strains combined in a dictionary
    """
    all_loss = {node: 0 for node in g.nodes()}
    for key in lost_genes:
        for gene in lost_genes[key]:
            if gene not in all_loss:  # happened to YP_008719742.1 for example
                continue
            all_loss[gene] += 1
    return all_loss


def get_gene_by_node(cds_dict, node_name):
    """
    This function get's the conversion dictionary and a node name and returns the original gene.
    @param cds_dict: the conversion dictionary
    @param node_name: the node name (such as NP..)
    @return: the original gene name
    """
    for gene, node in cds_dict.items():
        if node == node_name:
            return gene


def get_lost_and_exist_genes():
    if not os.path.exists(RAW_DATA_P):
        # create progenitor-progeny lost and exist genes ("raw data")
        proper_comparisons, proper_trials = get_proper_samples()
        clean_s1_df = clean_s1()
        comparisons_time_intervals = calc_time_intervals(clean_s1_df, proper_comparisons)
        raw_df = get_network_df(comparisons_time_intervals)
    else:
        raw_df = pd.read_csv(RAW_DATA_P)
    exist_genes = {}
    lost_genes = {}
    for _, row in raw_df.iterrows():
        exist_genes[(row["t1"], row["t2"])] = eval(row["t1_exist_genes"])
        lost_genes[(row["t1"], row["t2"])] = eval(row["lost_genes_t1_to_t2"])
    return lost_genes, exist_genes


def get_basic_network():
    """
    Create basic network from the data, run the propagation on the network
    @return:
    """
    tf_dict = data_parsing()
    converted_tf_dict = convert_network_to_refseq.main(tf_dict, CDS_FILE_P, EXTRA_FILES)  # convert the genes names
    node_to_gene = convert_network_to_refseq.parse_cds(CDS_FILE_P, EXTRA_FILES)
    locus_node_to_gene = convert_network_to_refseq.parse_locus_tag(CDS_FILE_P, EXTRA_FILES)
    # get lost genes {(t1, t2): [genes],...}, and exist gene {(t1, t2): [genes],...}
    lost_genes, exist_genes = get_lost_and_exist_genes()
    g = create_graph(converted_tf_dict)
    g = network_propagation(g, lost_genes, exist_genes)
    g = net_het(g)
    for node in g.nodes():
        calculate_weight_frac(g, node)
        g.node[node]['weight'] = np.mean(g.node[node]['weight_vec'])
        g.node[node]['name'] = get_gene_by_node(node_to_gene, node)
        g.node[node]['locus_tag'] = get_gene_by_node(locus_node_to_gene, node)
    return g, lost_genes, exist_genes


def exclude_gene_by_weight_vec_len(g):
    weights = [g.nodes[node]['weight'] for node in g.nodes() if len(g.nodes[node]['weight_vec']) > 10]
    min_weight = min(weights)
    max_weight = max(weights)
    for node in g.nodes():
        if len(g.node[node]['weight_vec']) < 10:
            if g.node[node]['weight'] < 0:
                g.node[node]['excluded_weight'] = max(min_weight, g.node[node]['weight'])
            else:
                g.node[node]['excluded_weight'] = min(max_weight, g.node[node]['weight'])
        else:
            g.node[node]['excluded_weight'] = g.node[node]['weight']
    return max_weight, min_weight


def main():
    g, lost_genes, exist_genes = get_basic_network()
    max_weight, min_weight = exclude_gene_by_weight_vec_len(g)
    all_loss = sum_lost_genes(g, lost_genes)
    export_to_csv(g, all_loss)
    draw_graph(g, max_weight, min_weight)


if __name__ == '__main__':
    main()
