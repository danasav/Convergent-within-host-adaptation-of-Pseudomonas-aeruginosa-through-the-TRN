"""
Each node in the network is associated with a vector,
indicating whether its function was up- or down-regulated in each strain.
In this file, I'll cluster the nodes,  and try to find if there are some nodes which tend to be lost together.
"""
import csv
import errno
import math
import operator
import scipy.stats
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from sklearn.metrics import silhouette_samples, silhouette_score
import seaborn as sns
import convert_network_to_refseq
from config import *
from main import data_parsing, create_graph
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
import numpy as np
from tqdm import tqdm

CORRELATION_DATA_PATH = "AnalysingClustering/correlated_pairs_by_weight_vector.csv"

NOT_EXIST_VAL = -7


# =========================================================================================================== #
# ========================================= Retrieve Clusters Data  ========================================= #
# =========================================================================================================== #

def get_PAO1_data():
    """
    This function creates a pandas dataframe from a csv file,
    and remove genes that haven't changed at all (weight vector equals to zero)
    """
    relevant_data = pd.read_csv(SUPPLEMENTARY_TABLE_1)[["Node name", "Weight Dict"]]
    strains_db = pd.read_csv(RAW_DATA_P)[["t1", "t2"]]
    all_strains = [(row["t1"], row["t2"]) for _, row in strains_db.iterrows()]

    data = {}
    for idx, row in relevant_data.iterrows():
        d = eval(row["Weight Dict"])
        # remove nodes with weight vector that contains only zeros
        vec = d.values()
        if all(v == 0 for v in vec) or len(vec) < 90:
            continue
        data[row["Node name"]] = {}
        for strain in all_strains:
            data[row["Node name"]][strain] = NOT_EXIST_VAL
            if strain in d:
                data[row["Node name"]][strain] = d[strain]

    return pd.DataFrame(data=data).transpose()


def create_correlation_data():
    """
    This function creates a dictionary that will be used to find all highly
    correlated genes, that are very far from each other.
    @return: a pandas dataframe of correlation, distance in the network, and intersection for each pair of nodes
    """
    relevant_data = pd.read_csv(SUPPLEMENTARY_TABLE_1)[["Node name", "Weight Dict"]]
    clean_df = relevant_data

    for idx, row in clean_df.iterrows():
        d = eval(row["Weight Dict"])
        # remove nodes with weight vector that contains only zeros
        vec = d.values()
        if all(v == 0 for v in vec) or len(vec) < 90:
            clean_df = clean_df.drop(idx, axis=0)

    df_clean = clean_df.set_index('Node name').transpose()
    G = get_undirected_regulatory_network()
    correlation_data = {"NodeA": [], "NodeB": [], "Correlation": [],
                        'Pval': []}  # , 'Distance': [], 'Intersection': [], 'LowestCommonAncestor':[]}
    for gene in tqdm(df_clean.columns):
        for other_gene in df_clean.columns:
            if gene == other_gene or (other_gene in correlation_data["NodeA"] and gene in correlation_data["NodeB"]):
                continue

            # create proper vectors from dictionary:
            gene_dict = eval(df_clean[gene]["Weight Dict"])
            other_gene_dict = eval(df_clean[other_gene]["Weight Dict"])

            strains_intersection = set(gene_dict.keys()).intersection(set(other_gene_dict.keys()))

            gene_vec = []
            other_gene_vec = []

            for comp in strains_intersection:
                gene_vec.append(gene_dict[comp])
                other_gene_vec.append(other_gene_dict[comp])

            rho, pval = scipy.stats.spearmanr(gene_vec, other_gene_vec)
            if math.isnan(rho):
                # i do clean the vectors containing only zeros,
                # but i after we take only the intersected keys, some of them becomes all zeros,
                # which creates a NAN correlation
                continue

            correlation_data["NodeA"].append(gene.strip())
            correlation_data["NodeB"].append(other_gene.strip())
            correlation_data["Correlation"].append(rho)
            correlation_data["Pval"].append(pval)
            # correlation_data["Distance"].append(measure_distance(G, gene, other_gene))
            # correlation_data["Intersection"].append(measure_intersection(G, gene, other_gene))
            # correlation_data["LowestCommonAncestor"].append(measure_lowest_common_ancestor(G, gene, other_gene))

    corr_data_frame = pd.DataFrame.from_dict(correlation_data)
    corr_data_frame.to_csv(CORRELATION_DATA_PATH)
    return corr_data_frame


def measure_distance(g, node_a, node_b):
    """
    @param g - should be undirected graph ( g = orig_g.to_undirected())
    This function measures the distance between node_a and node_b in the network.
    If they are not connected at all (there isn't a path between them), it returns infinity
    """
    try:
        return nx.shortest_path_length(g, source=node_a.strip(), target=node_b.strip())
    except nx.exception.NetworkXNoPath:
        return float('inf')


def measure_lowest_common_ancestor(g, node_a, node_b):
    """
    @param g - should be undirected graph ( g = orig_g.to_undirected())
    This function measures the distance between node_a and node_b in the network.
    If they are not connected at all (there isn't a path between them), it returns infinity
    """
    try:
        LCA_from_A = (None, None)
        LCA_from_B = (None, None)
        node_a_common_ancestors = sorted(
            {an: nx.shortest_path_length(g, source=an, target=node_a) for an in nx.ancestors(g, node_a)}.items(),
            key=operator.itemgetter(1))
        node_b_common_ancestors = sorted(
            {an: nx.shortest_path_length(g, source=an, target=node_b) for an in nx.ancestors(g, node_b)}.items(),
            key=operator.itemgetter(1))
        for ancestor, distance in node_a_common_ancestors:
            node_b_ancestors, node_b_distances = zip(*node_b_common_ancestors)
            if ancestor in node_b_ancestors:
                node_b_dist = node_b_distances[node_b_ancestors.index(ancestor)]
                LCA_from_A = (ancestor, (distance + node_b_dist) / 2.0)
                break
        for ancestor, distance in node_b_common_ancestors:
            node_a_ancestors, node_a_distances = zip(*node_a_common_ancestors)
            if ancestor in node_a_ancestors:
                node_b_dist = node_a_distances[node_a_ancestors.index(ancestor)]
                LCA_from_B = (ancestor, (distance + node_b_dist) / 2.0)
                break
        LCAs = [LCA_from_A, LCA_from_B]
        return LCAs[np.argmin([LCA_from_A[1], LCA_from_B[1]])]
    except ValueError:
        return (None, None)


def measure_intersection(g, node_a, node_b):
    """
    @param g - should be undirected graph ( g = orig_g.to_undirected())
    Node connectivity is equal to the minimum number of nodes that
    must be removed to disconnect G or render it trivial.
    If the nodes are not connected at all (there isn't a path between them), it returns infinity
    """
    try:
        return nx.node_connectivity(g, node_a.strip(), node_b.strip())
    except nx.exception.NetworkXNoPath:
        return float('inf')


def get_undirected_regulatory_network():
    """
    This function assembling the data for the correlated pairs dictionary
    @return: G - undirected regulatory network graph,
    """
    tf_dict = data_parsing()
    converted_tf_dict = convert_network_to_refseq.main(tf_dict, CDS_FILE_P, EXTRA_FILES)  # convert the genes names
    g = create_graph(converted_tf_dict)
    G = g.to_undirected()
    return G


def get_correlation_data():
    """
    This function gets the correlation data  from the correlated pairs csv file
    @return: pandas' dataFrame where each row is : NodeA, NodeB, Correlation, Pval, Distance
    """
    path = CORRELATION_DATA_PATH
    if not os.path.exists(path):
        create_correlation_data()
    return pd.read_csv(path)


# =========================================================================================================== #
# ============================================ Clusters Analysis ============================================ #
# =========================================================================================================== #

def create_distance_correlation_histogram(zoomed=True):
    """
    This function plots for each correlation value, it's distances histogram
    We want to use this function in order to understand better which threshold to choose for the clustering.
    For example, if we will see that most of the pairs with distance > 3 have correlation value of 0.6,
    we would probably like to take that correlation threshold for the clustering
    """
    correlation_data = get_correlation_data()
    if zoomed:
        bins = np.arange(0.39, 0.61, 0.01)
        print "bins", bins
        distance_for_corr = {round(bin, 2): [] for bin in bins}
        print "distance_for_corr", distance_for_corr
    else:
        bins = np.arange(0.1, 1.1, 0.1)
        print "bins", bins
        distance_for_corr = {round(bin, 1): [] for bin in bins}
        print "distance_for_corr", distance_for_corr
    distances_vals = set()
    for _, row in correlation_data.iterrows():
        if zoomed:
            correlation = abs(round(row["Correlation"], 2))
        else:
            correlation = abs(round(row["Correlation"], 1))
        distance = row["Distance"]
        if distance == np.inf:
            continue
        distances_vals.add(distance)
        if correlation in distance_for_corr:
            distance_for_corr[correlation].append(distance)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    cm = plt.get_cmap('Pastel1')
    num_colors = len(distance_for_corr.keys())
    ax.set_color_cycle([cm(1. * i / num_colors) for i in range(num_colors, 0, -1)])
    bins = sorted(distances_vals)
    for corr in sorted(distance_for_corr.keys()):
        x = distance_for_corr[corr]
        ax.hist(x, bins, label="correlation value = {0}".format(corr))
    ax.set_ylabel("Pairs with this correlation value", size=15)
    ax.set_xlabel("Distances in the original network", size=15)
    ax.set_xticks(bins)
    # ax.set_title("Distances histogram for each correlation value")

    lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    if zoomed:
        fig.savefig("AnalysingClustering/Distance histograms for different correlation values zoomed.png",
                    bbox_extra_artists=(lgd,), bbox_inches='tight')
    else:
        fig.savefig("AnalysingClustering/Distance histograms for different correlation values.png",
                    bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.show()
    return distance_for_corr


def draw_cluster_map(correlation_data):
    """
    create seaborn cluster map for visualization
    :param corrlation_data : pandas' dataFrame where each row is : NodeA, NodeB, Correlation, Pval, Distance
    @return: cleaned dataframe (each column represents the gene's weight vector) and the correlation metrix
    """
    rho = create_correlation_matrix(correlation_data).values
    cmap = plt.cm.RdPu
    sns.clustermap(abs(np.nan_to_num(rho)), cmap=cmap)  # clustermap based of the correlation matrix
    plt.savefig('AnalysingClustering/Hierarchical clustering.pdf')
    plt.show()
    return rho


def create_correlation_matrix(correlation_data):
    """
    creates correlation matrix - a dataframe were each cell is the correlation value between the column's node and the row's node
    :param correlation_data: pandas' dataFrame where each row is : NodeA, NodeB, Correlation, Pval, Distance
    :return: correlation matrix, full or diagonal
    """
    pairs_corr = correlation_data[["NodeA", "NodeB", "Correlation"]]
    data = {}  # NodeA : {NodeB: corr,..}
    for _, row in pairs_corr.iterrows():
        if row["NodeA"] not in data:
            data[row["NodeA"]] = {row["NodeA"]: 1}
        if row["NodeB"] not in data:
            data[row["NodeB"]] = {row["NodeB"]: 1}
        data[row["NodeA"]][row["NodeB"]] = row["Correlation"]
        data[row["NodeB"]][row["NodeA"]] = row["Correlation"]
    return pd.DataFrame(data=data)


def show_cluster_in_the_network(clusters, relevant_cluster_id):
    """
    This function is used to visualize only one relevant cluster. It shows in the original network.
    @param clusters: a dictionary of the cluster idx andAnalysingClustering it's nodes
    @param relevant_cluster_id: the relevant cluster id after deletion of irrelevant clusters (contaning only 1 node)
     and after sorting
    """

    tf_dict = data_parsing()
    converted_tf_dict = convert_network_to_refseq.main(tf_dict, CDS_FILE_P, EXTRA_FILES)  # convert the genes names
    g = create_graph(converted_tf_dict)
    node_to_gene = convert_network_to_refseq.parse_cds(CDS_FILE_P, EXTRA_FILES)
    locus_node_to_gene = convert_network_to_refseq.parse_locus_tag(CDS_FILE_P, EXTRA_FILES)
    locus_node_to_gene.update(node_to_gene)
    node_colors = []
    nodes_in_cluster = []
    for node in g.nodes():
        node_cluster = [i for i in clusters if node in clusters[i]]  # get the node's cluster
        if len(node_cluster) == 0:
            node_colors.append('darkgray')
            continue
        node_cluster = node_cluster[0]  # the list contains only one cluster
        if node_cluster == relevant_cluster_id:
            # if node in group_one:
            #     node_colors.append(colors[0])
            # elif node in group_two:
            #     node_colors.append(colors[1])
            # elif node in group_three:
            #     node_colors.append(colors[2])
            # else:
            nodes_in_cluster.append(node)
            node_colors.append('seagreen')
        else:
            node_colors.append('darkgray')

    pos_regulation_edges = []
    neg_regulation_edges = []
    nodes_to_draw = []
    for u, v in g.edges():
        if u in nodes_in_cluster or v in nodes_in_cluster:
            nodes_to_draw.append(u)
            nodes_to_draw.append(v)
            if g[u][v]['color'] == UP_REG_COLOR:
                pos_regulation_edges.append((u, v))
            elif g[u][v]['color'] == DOWN_REG_COLOR:
                pos_regulation_edges.append((u, v))

    f = plt.figure(figsize=(19.20, 10.80))
    ax = f.add_subplot(1, 1, 1)
    pos = nx.nx_agraph.graphviz_layout(g)
    nx.draw_networkx_nodes(g, pos=pos, node_size=19, with_labels=False, node_color=node_colors, nodelist=nodes_to_draw)
    nx.draw_networkx_edges(g, pos=pos, node_size=19, edgelist=pos_regulation_edges, arrowsize=15, edge_color='darkgray')
    nx.draw_networkx_edges(g, pos=pos, node_size=19, edgelist=neg_regulation_edges, arrowstyle='-[',
                           arrowsize=7, edge_color='darkgray')

    plt.axis('off')
    for spine in plt.gca().spines.values():
        spine.set_visible(False)

    # s = lambda m, c: plt.plot([], [], marker=m, color=c, ls="none")[0]
    # handles = [s("s", colors[i]) for i in range(3)]
    # labels = ["Group I", "Group II", "Group III"]
    # plt.legend(handles, labels)

    plt.title("Cluster {0}".format(relevant_cluster_id))
    f.savefig("{0}_cluster.png".format(relevant_cluster_id), dpi=600)
    plt.show()
    plt.close()


def color_clusters_in_network(clusters, corrVal=0.7, type_of_clustering='Hierarchical'):
    """
    This function is used to visualize the clustering in the network.
    Each node will be colored by it's cluster color.
    @param clusters: a dictionary of cluster number and it's nodes
    """
    tf_dict = data_parsing()
    converted_tf_dict = convert_network_to_refseq.main(tf_dict, CDS_FILE_P, EXTRA_FILES)  # convert the genes names
    g = create_graph(converted_tf_dict)

    colors = [
        '#FFFF00', '#1CE6FF', '#FF34FF', '#FF4A46', '#008941', '#006FA6', '#A30059', '#FFDBE5', '#7A4900', '#0000A6',
        '#63FFAC', '#B79762', '#004D43', '#8FB0FF', '#997D87', '#5A0007', '#1B4400', '#4FC601', '#3B5DFF', '#4A3B53',
        '#FF2F80', '#BA0900', '#6B7900', '#00C2A0', '#FFAA92', '#FF90C9', '#B903AA', '#D16100', '#4C325D', '#7B4F4B',
        '#FFB500', '#98D058', '#A1C299', '#300018', '#0AA6D8', '#013349', '#00846F', '#372101']

    colors_labels = {cluster: colors[i] for i, cluster in enumerate(clusters)}
    node_colors = []
    for node in g.nodes():
        node_cluster = [i for i in clusters if node in clusters[i]]  # get the node's cluster
        if len(node_cluster) == 0:
            colors_labels["Not Clustered"] = 'darkgray'
            node_colors.append('darkgray')
            continue
        node_cluster = node_cluster[0]  # the list contains only one cluster
        node_colors.append(colors_labels[node_cluster])

    pos_regulation_edges = [(u, v) for u, v in g.edges() if g[u][v]['color'] == UP_REG_COLOR]
    neg_regulation_edges = [(u, v) for u, v in g.edges() if g[u][v]['color'] == DOWN_REG_COLOR]

    f = plt.figure(figsize=(19.20, 10.80))
    ax = f.add_subplot(1, 1, 1)
    for label in sorted(colors_labels):
        if type(label) != str:
            ax.plot([0], [0], color=colors_labels[label], label="Cluster {0}".format(label))
        else:
            ax.plot([0], [0], color=colors_labels[label], label=label)

    pos = nx.nx_agraph.graphviz_layout(g)
    nx.draw_networkx_nodes(g, pos=pos, node_size=19, with_labels=False, node_color=node_colors)
    nx.draw_networkx_edges(g, pos=pos, node_size=19, edgelist=pos_regulation_edges, arrowsize=15, edge_color='darkgray')
    nx.draw_networkx_edges(g, pos=pos, node_size=19, edgelist=neg_regulation_edges, arrowstyle='-[',
                           arrowsize=7, edge_color='darkgray')

    plt.axis('off')
    for spine in plt.gca().spines.values():
        spine.set_visible(False)

    plt.title("{1} Clustering, using the threshold correlation {0}".format(corrVal, type_of_clustering))
    lgd = ax.legend(bbox_to_anchor=(0.95, 1), loc=2, borderaxespad=0.)
    f.savefig("Visualised_cluster_by_correlation_{0}.png".format(corrVal, type_of_clustering),
              bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=300)
    plt.show()
    plt.close()


def get_gene_name_to_node(node):
    """
    For a given node in the network, returns the gene's name and description.
    """
    node_to_gene_df = pd.read_csv(SUPPLEMENTARY_TABLE_1)
    return node_to_gene_df[node_to_gene_df['Node name'] == node]['Gene name'].values[0]


def create_silhouette(clusters_dict, labels, df_clean, threshold=0.91):
    """
    Create a silhouette for the given clusters.
    @param clusters_dict: a dictionary of the cluster idx and all the nodes it contains
    @param labels: linkage matrix - should be the labels of the clusters, an array in size of the number_of_nodes.
    @param df_clean: the original pandas data frame, where each row is the node's "weight vector"
    @param threshold: current threshold for correlation
    """
    sample_silhouette_values = silhouette_samples(df_clean.values, labels, metric=abs_corr)
    # the silhouette average score.
    # If none of the clusters silhouette score is bigger then the average, we can tell the clustering isn't good
    silhouette_avg = silhouette_score(df_clean.values, labels, metric=abs_corr)
    plot_silhouette(clusters_dict, labels, sample_silhouette_values, silhouette_avg, threshold)


def plot_silhouette(clusters, linkage_mat, sample_silhouette_values, silhouette_avg, threshold):
    """
    Plots the silhouette for the clusters
    @param clusters: a dictionary of the cluster idx and all the nodes it contains
    @param sample_silhouette_values: The Silhouette Coefficient (a measure of how well samples are clustered with
    samples that are similar to themselves) for each sample.
    The best value is 1 and the worst value is -1. Values near 0 indicate overlapping clusters.
    @param silhouette_avg: The mean Silhouette Coefficient of all samples.
    @param threshold: current threshold for correlation
    """
    clusters_far_from_silhouette = 0  # if highest silhouette value in the cluster is smaller then the average - this value should be +1.
    # remove clusters containing only one element:
    clusters = {key: clusters[key] for key in clusters.keys() if len(clusters[key]) > 1}
    y_lower = 10
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cm = plt.get_cmap('Pastel2')
    num_colors = len(clusters.keys())
    ax.set_color_cycle([cm(1.5 * i / num_colors) for i in range(num_colors)])
    for key in clusters.keys():
        # Aggregate silhouette scores for samples belonging to the cluster and sort them
        ith_cluster_silhouette_values = sample_silhouette_values[linkage_mat == key]
        ith_cluster_silhouette_values.sort()
        if all(abs(ith_cluster_silhouette_values) < abs(silhouette_avg)):
            clusters_far_from_silhouette += 1
        y_upper = y_lower + len(clusters[key])
        plt.fill_betweenx(np.arange(y_lower, y_upper),
                          0, ith_cluster_silhouette_values,
                          label=str(key))
        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    # The vertical line for average silhouette score of all the values
    ax.axvline(x=silhouette_avg, color="red", linestyle="--")
    # plt.yticks([])  # Clear the yaxis labels / ticks

    ax.set_ylabel("Cluster labels ({0} clusters)".format(num_colors))
    ax.set_xlabel("The silhouette coefficient values")
    ax.set_yticks([])
    title = "Silhouette plot \n Threshold = {1}, Number of clusters = {2}\n clusters far from silhouette average = {3}"
    ax.set_title(title.format(
        threshold, num_colors, clusters_far_from_silhouette))

    lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    fig.savefig("AnalysingClustering/Silhouette/silhouette_for_clustering_with_threshold_{0}.png".format(
        threshold), bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close()


def export_clustering_to_csv(clusters, th=''):
    """
    This function exports the hierarchical clustering data (number of clusters and information about each cluster)
    into csv file
    @param clusters: the clusters created by the hierarchical clustering
    @param th: the correlation threshold for the clustering (tree cutoff)
    """
    csv_columns = ['Cluster Id', 'Number of Nodes', 'Nodes in Cluster', 'Genes in Cluster',
                   "Cluster Corr Avg", 'Cluster distances', 'Clusters LCAs', 'Number of Positively Correlated',
                   'Number of Negatively Correlated']
    csv_file = \
        "AnalysingClustering/{0}_clustering_over_data_with_corr_metric.csv".format(th)
    try:
        if not os.path.exists(os.path.dirname(csv_file)):
            try:
                os.makedirs(os.path.dirname(csv_file))
            except OSError as exc:  # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise
        with open(csv_file, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
            writer.writeheader()
            for cluster in clusters:
                corr_avg, dist, num_of_pos, num_of_neg, LCAs = get_clustering_identifiers(cluster, clusters)
                data = {"Cluster Id": cluster, "Number of Nodes": len(clusters[cluster]),
                        'Nodes in Cluster': list(clusters[cluster]),
                        'Genes in Cluster': [get_gene_name_to_node(node) for node in clusters[cluster]],
                        "Cluster Corr Avg": corr_avg, 'Cluster distances': dist, 'Clusters LCAs': LCAs,
                        'Number of Positively Correlated': num_of_pos, 'Number of Negatively Correlated': num_of_neg}
                writer.writerow(data)
    except IOError, e:
        print("I/O error")
        print(e)


def get_clustering_identifiers(cluster, clusters):
    """
    This function returns identifiers for the current clustering method
    @return the average correlation, the distances in the original network, number of positively correlated and
    negatively correlated nodes in the cluster
    """
    corr_data = get_correlation_data()
    num_of_pos = 0
    num_of_neg = 0
    pairs_correlation = {}
    pairs_distances = {}
    pairs_LCAs = {}
    for node in clusters[cluster]:
        for neighbor_node in clusters[cluster]:
            if neighbor_node == node or (neighbor_node, node) in pairs_correlation:
                continue
            node_corr = corr_data[corr_data['NodeA'] == node.strip()]
            pair_data = node_corr[node_corr['NodeB'] == neighbor_node.strip()]
            if len(pair_data) != 0:  # there are cases when only the "reversed" pair appear
                corr_val = pair_data.Correlation.iloc[0]
                if corr_val > 0:
                    num_of_pos += 1
                else:
                    num_of_neg += 1
                pairs_correlation[(node, neighbor_node)] = abs(corr_val)
                pairs_distances[(node, neighbor_node)] = pair_data.Distance.iloc[0]
                pairs_LCAs[(node, neighbor_node)] = pair_data.LowestCommonAncestor.iloc[0]
    return np.average(pairs_correlation.values()), pairs_distances, num_of_pos, num_of_neg, pairs_LCAs


def concordance_for_clusters(df_clean, new_clusters):
    node_weights = pd.read_csv(SUPPLEMENTARY_TABLE_1)[["Node name", "Weight"]].set_index('Node name')
    node_weights["Weight"] = np.sign(node_weights["Weight"])

    df = df_clean.replace(-7, np.nan)
    new_df = df.T * df.join(node_weights)["Weight"]
    new_df = new_df.replace(-0.0, 0.0).T

    df = new_df.reset_index()
    all_mean_dfs = []
    all_concordance_dfs = []
    for i in new_clusters.keys():
        c = df[df['index'].isin(new_clusters[i])]
        mean_df = c.mean(skipna=True) / df.mean(skipna=True)
        mean_df = mean_df.to_frame()
        mean_df = mean_df.rename(columns={0: "cluster {0}".format(i)})
        mean_df.to_csv("concordance/clusters_mean_vals_concordance_normalized_cluster_{0}.csv".format(i))  # not abs
        all_mean_dfs.append(mean_df)
        # concordance without normalization!
        concordance_df = c.mean(skipna=True)
        concordance_df = concordance_df.to_frame()
        concordance_df = concordance_df.rename(columns={0: "cluster {0}".format(i)})
        concordance_df.to_csv("concordance/clusters_mean_vals_concordance_cluster_{0}.csv".format(i))
        all_concordance_dfs.append(concordance_df)

        cluster = "cluster {0}".format(i)
        sns.distplot(mean_df[cluster][~mean_df[cluster].isin([np.nan, np.inf, -np.inf])].values,
                     label=cluster, kde=False, bins=200)
        plt.legend()
        plt.savefig("concordance/distplot_concordance_normalized_cluster_{0}.png".format(i))
        plt.close()

        sns.distplot(concordance_df[cluster][~concordance_df[cluster].isin([np.nan, np.inf, -np.inf])].values,
                     label=cluster, kde=False, bins=100)
        plt.legend()
        plt.savefig("concordance/distplot_concordance_cluster_{0}.png".format(i))
        plt.close()

    pd.concat(all_mean_dfs, sort=True).to_csv("clusters_mean_vals_concordance_normalized.csv")  # not abs
    pd.concat(all_concordance_dfs, sort=True).to_csv("clusters_mean_vals_concordance.csv")  # not abs


# =========================================================================================================== #
# ========================================= Hierarchical clustering ========================================= #
# =========================================================================================================== #


def abs_corr(node_a_weights, node_b_weights):
    """
    @return: The absolute spearman correlation between the two nodes by their dictionary
    """
    mask = [(node_a_weights != NOT_EXIST_VAL) & (node_b_weights != NOT_EXIST_VAL)]
    node_a_vec = node_a_weights[mask]
    node_b_vec = node_b_weights[mask]
    rho, _ = scipy.stats.spearmanr(node_a_vec, node_b_vec)
    if math.isnan(rho):
        # print node_a_vec, node_b_vec
        rho = 0
    return 1 - abs(rho)


def hierarchical_clustering(genes, df_clean, threshold=0.91):
    """
    This function preforms hierarchical clustering over the cluster data, with correlation metric and average method
    (UPGMA method - check the nearest neighbor in the tree)
    @:param genes - the genes names for labeling the dendrogram
    @:param df_clean - each row represents the node's "weight" for each "strain". if the node didn't exist before,
     it's value would be -7 for the specific strain
    @:param threshold: the threshold for clustering (Cutting the tree by this threshold)
    @return: created clusters
    """
    linkage_matrix_file = 'linkage_matrix.npy'
    if os.path.isfile(linkage_matrix_file):
        linkage_matrix = np.load(linkage_matrix_file)
    # creating the linkage matrix - a numpy ndarray where each row is the correlation,diagonal matrix
    else:
        linkage_matrix = linkage(df_clean, metric=abs_corr, method='weighted')
        np.save(linkage_matrix_file, linkage_matrix)

    # here the correlation isn't based of absolute values
    plt.figure(figsize=(25, 10))
    plt.title('Hierarchical Clustering Dendrogram - Correlation distance')
    plt.ylabel('Similarity (1 - correlation coefficient)', fontsize=20)
    dendrogram(linkage_matrix, labels=genes)
    clusters = fcluster(linkage_matrix, threshold, criterion='distance')
    plt.savefig("Linkage Over Data, Correlation metric.png")
    plt.tick_params(labelsize=15)
    plt.close()
    return clusters


def determine_tree_cutoff(df_clean):
    """
    Determines the tree best cutoff using silhouette for different cutoff thresholds
    @:param df_clean - each row represents the node's "weight" for each "strain". if the node didn't exist before,
    it's value would be -7 for the specific strain
    """
    corr_vals = np.arange(0.1, 0.98, 0.1)
    # corr_vals = np.arange(0.4, 0.61, 0.01)
    num_of_clusters = []
    # creating the linkage matrix - a numpy ndarray where each row is the correlation,diagonal matrix
    linkage_matrix_file = 'linkage_matrix.npy'
    if os.path.isfile(linkage_matrix_file):
        linkage_matrix = np.load(linkage_matrix_file)
    else:
        linkage_matrix = linkage(df_clean, metric=abs_corr, method='weighted')
        np.save(linkage_matrix_file, linkage_matrix)
    genes = df_clean.index.values.T

    for i in corr_vals:
        labels = fcluster(linkage_matrix, i, criterion='distance')
        clusters_dict = {cluster: genes[labels == cluster] for cluster in set(labels)}
        num_of_clusters.append(len(clusters_dict.keys()))
        print "finished for th: ", i, " with ", len(clusters_dict.keys()), " clusters"
    plt.plot(corr_vals, num_of_clusters)
    plt.xlabel("Different correlation threshold for clustering")
    plt.ylabel("Number of clusters")
    plt.title("Correlation threshold vs. Number of clusters")
    plt.savefig("AnalysingClustering/corr_threshold_vs_numOfClusters.png")
    plt.close()


def determine_silhouette(df_clean):
    """
    Determines the tree best cutoff using silhouette for different cutoff thresholds
    @:param df_clean - each row represents the node's "weight" for each "strain". if the node didn't exist before,
    it's value would be -7 for the specific strain
    """
    corr_vals = np.arange(0.1, 0.98, 0.1)
    # corr_vals = np.arange(0.4, 0.61, 0.01)
    # corr_vals = np.arange(0.21, 0.51, 0.01)
    num_of_clusters = []
    avgs = []

    linkage_matrix_file = 'linkage_matrix.npy'
    if os.path.isfile(linkage_matrix_file):
        linkage_matrix = np.load(linkage_matrix_file)
    # creating the linkage matrix - a numpy ndarray where each row is the correlation,diagonal matrix
    else:
        linkage_matrix = linkage(df_clean, metric=abs_corr, method='weighted')
        np.save(linkage_matrix_file, linkage_matrix)

    f = open("sil_avg.txt", "a+")
    for i in corr_vals:
        labels = fcluster(linkage_matrix, i, criterion='distance')
        genes = df_clean.index.values.T
        clusters_dict = {cluster: genes[labels == cluster] for cluster in set(labels)}
        clusters_len = len([key for key in clusters_dict.keys() if len(clusters_dict[key]) > 1])
        num_of_clusters.append(clusters_len)
        sil_avg_score = silhouette_score(df_clean.values, labels, metric=abs_corr)
        avgs.append(sil_avg_score)
        f.writelines(["{0} : num_of_clusters={1} ,avg={2}\n".format(i, clusters_len, sil_avg_score)])
        print("{0} : num_of_clusters={1} ,avg={2}\n".format(i, clusters_len, sil_avg_score))
    f.close()

    plot_determination()


def plot_determination():
    corr_vals = np.arange(0.1, 0.98, 0.1)
    # corr_vals = np.arange(0.4, 0.61, 0.01)
    num_of_clusters = []
    avgs = []

    f = open("sil_avg.txt", 'r')
    for line in f.readlines():
        line = line.split(":")[1].split(",")
        num_of_clusters.append(int(line[0].split("=")[1].strip()))
        avgs.append(float(line[1].split("=")[1].strip()))
    f.close()

    a = plt.scatter(corr_vals, avgs, c=num_of_clusters, label=num_of_clusters, cmap="winter")
    plt.colorbar(a, label='number of classes')
    plt.xlabel("Threshold values", size="20")
    plt.ylabel("Silhouette averages", size="20")
    plt.savefig("sil_avgs.png")
    plt.show()

    a = plt.scatter(corr_vals, num_of_clusters, c=avgs, cmap="cool")
    plt.colorbar(a, label='Silhouette averages')
    plt.xlabel("Threshold values", size="20")
    plt.ylabel("Number of clusters", size="20")
    plt.savefig("elbow_plot.png")
    plt.show()


def cluster_by_comparisons(df_clean):
    f = plt.figure(figsize=(10, 12))
    comparisons = list(df_clean.columns)
    linkage_matrix = linkage(df_clean.T, metric=abs_corr, method='weighted')
    dendrogram(linkage_matrix, labels=comparisons, orientation='left', leaf_font_size=8, truncate_mode='level')
    high_concordance_score_14 = [
        ("ERR1189823", "ERR1189824"), ("ERR1538256", "ERR1538257"),  ("ERR2286126", "ERR2286128"),
        ("ERR2322162", "ERR2322163"), ("ERR2322170", "ERR2322171"),  ("ERR2322211", "ERR2322202"),
        ("ERR2322211", "ERR2322203"), ("ERR2322228", "ERR2322233"),  ("ERR2322235", "ERR2322236"),
        ("ERR2322240", "ERR2322241"), ("ERR2322243", "ERR2322245"),  ("ERR2322245", "ERR2322246"),
        ("ERR431085", "ERR431096"), ("ERR431094", "ERR431092"), ("ERR431094", "ERR431093"), ("ERR431106", "ERR431104"),
        ("ERR431106", "ERR431105"), ("ERR431113", "ERR431110"), ("ERR431113", "ERR431112"), ("ERR431186", "ERR431183"),
        ("ERR431188", "ERR431189"), ("ERR431350", "ERR431364"), ("ERR431371", "ERR431350"), ("ERR705703", "ERR705706"),
        ("SRR3121771", "SRR3121774"), ("SRR3121805", "SRR3121809"), ("SRR4097555", "SRR4097742"),
        ("SRR4097555", "SRR4097763"), ("SRR3121707", "SRR3121710"), ("ERR2322228", "ERR2322229"),
        ("ERR2322214", "ERR2322216"), ("ERR2322235", "ERR2322237"), ("r6p18t7",	"Sommer_setup_r6p18t6"),
        ("ERR2286126", "ERR2286129"), ("ERR431373", "ERR431351"), ("ERR2322166", "ERR2322167"),
        ("ERR2322218", "ERR2322219"), ("SRR3121707", "SRR3121761"), ("ERR1959757", "ERR1959753"),
        ("ERR2322207", "ERR2322208"), ("SRR3121707", "SRR3121772"), ("ERR953519", "ERR1959748"),
        ("ERR953525", "ERR1959756"), ("ERR953525", "ERR1959757"), ("SRR3121731", "SRR3121732"),
        ("ERR1189818", "ERR1189831")
    ]
    colors = ['g']
    ax = plt.gca()
    xlbls = ax.get_ymajorticklabels()
    for lbl in xlbls:
        if eval(lbl.get_text()) in high_concordance_score_14:
            lbl.set_color(colors[0])
        else:
            lbl.set_color('black')
    s = lambda m, c: plt.plot([], [], marker=m, color=c, ls="none")[0]
    handles = [s("s", colors[i]) for i in range(len(colors))]
    labels = ["High Concordance Score for cluster 14"]
    plt.legend(handles, labels)
    plt.tight_layout()
    f.savefig("Clustering comparisons (concordance 0.75)-truncated.png", dpi=500, bbox_inches='tight')
    plt.show()


# =========================================================================================================== #
# ============================================== Main functions ============================================= #
# =========================================================================================================== #


def choose_clustering_threshold(correlation_data, df_clean):
    draw_cluster_map(correlation_data)
    determine_tree_cutoff(df_clean)
    determine_silhouette(df_clean)
    create_distance_correlation_histogram(zoomed=False)
    create_distance_correlation_histogram(zoomed=True)


def cluster_specific_threshold(df_clean, threshold):
    genes = list(df_clean.index)
    labels = hierarchical_clustering(genes, df_clean, threshold=threshold)
    genes = df_clean.index.values.T
    clusters_dict = {cluster: genes[labels == cluster] for cluster in set(labels)}
    sorted_clusters = sorted(clusters_dict, key=lambda key: len(clusters_dict[key]))
    i = 1
    clusters = {}
    for cluster in sorted_clusters:
        if len(clusters_dict[cluster]) > 1:
            clusters[i] = clusters_dict[cluster]
            i += 1
    export_clustering_to_csv(clusters, th=threshold)
    create_silhouette(clusters_dict, labels, df_clean, threshold=threshold)
    color_clusters_in_network(clusters, threshold)
    show_cluster_in_the_network(clusters, threshold)
    return clusters


def main():
    df_clean = get_PAO1_data()

    # cluster by the comparisons
    cluster_by_comparisons(df_clean)

    # choose threshold:
    correlation_data = get_correlation_data()
    choose_clustering_threshold(correlation_data, df_clean)

    # after choosing the threshold, analyse the results
    new_clusters = cluster_specific_threshold(df_clean, 0.6)
    concordance_for_clusters(df_clean, new_clusters)


if __name__ == '__main__':
    main()
