import convert_network_to_refseq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection
from config import *


def go_converter():
    """
    Convert the go numbers to go descriptions
    @return: the converter dictionary
    """
    convert_data = open(CONVERT_GO_PATH, "r")
    go_convert_dict = {}
    data = convert_data.readlines()
    i = 0
    while i < len(data):
        if len(data[i].split()) > 0:
            if data[i].split()[0] == "id:":
                go = data[i].split()[1]
                i += 1
                go_convert_dict[go] = data[i].split(": ")[1][:-1]
        i += 1
    return go_convert_dict


def go_data_parsing():
    """
    Parse the GO data file, the  and the gene convert names file
    @return: go_dict: a dictionary of pathways and their gene numbers, go_dict_with_names: a dictionary of annotations
    and their gene names, genes_in_pathways: a set of all the genes that take part in these go annotation.
    """
    gene_convert_dict = convert_network_to_refseq.parse_cds(CDS_FILE_P, EXTRA_FILES)
    go_convert_dict = go_converter()
    go_dict = {}
    data = open(GO_PATH, "r")
    for row in data.readlines():
        row = row.split()
        for gene_name in gene_convert_dict:
            if gene_convert_dict[gene_name] == row[0]:
                for go in row[1].split(","):
                    if go in go_convert_dict:
                        go = go_convert_dict[go]
                    if go in go_dict:
                        if row[0] not in go_dict[go]:
                            go_dict[go].append(row[0])
                    else:
                        go_dict[go] = [row[0]]
                break
    return go_dict


def calc_scores(go_dict):
    """
    Calc the scores (mean heterogeneity, mean mROS) of each GO annotation (and the genes that are not in the pathway).
    """
    df = pd.DataFrame(columns=go_dict.keys()).T
    pao1_df = pd.read_csv(SUPPLEMENTARY_TABLE_1)
    for go in go_dict:
        go_het_scores = []
        non_go_het_scores = []
        go_reg_scores = []
        non_go_reg_scores = []
        for i, gene in enumerate(pao1_df['Node name']):
            if gene in go_dict[go]:
                if pao1_df.loc[i, 'Weight'] != 0:
                    go_het_scores.append(pao1_df.loc[i, 'Heterogeneity'])
                if str(pao1_df.loc[i, 'Weight']) != 'nan':
                    go_reg_scores.append(pao1_df.loc[i, 'Weight'])
            else:
                if pao1_df.loc[i, 'Weight'] != 0:
                    non_go_het_scores.append(pao1_df.loc[i, 'Heterogeneity'])
                if str(pao1_df.loc[i, 'Weight']) != 'nan':
                    non_go_reg_scores.append(pao1_df.loc[i, 'Weight'])
        df.loc[go, 'go_mean_het'] = np.mean(go_het_scores)
        df.loc[go, 'non_go_mean_het'] = np.mean(non_go_het_scores)
        df.loc[go, 'num_of_genes'] = len(go_het_scores)
        df.loc[go, 'go_mROS'] = np.mean(go_reg_scores)
        df.loc[go, 'non_go_mROS'] = np.mean(non_go_reg_scores)
        _, het_p_value = mannwhitneyu(go_het_scores, non_go_het_scores, alternative='greater')
        _, reg_p_value = mannwhitneyu(go_reg_scores, non_go_reg_scores, alternative='two-sided')
        df.loc[go, 'het_p_value'] = het_p_value
        df.loc[go, 'reg_p_value'] = reg_p_value
    df = df[df["num_of_genes"] > 1]
    df['het_rejected'], df['het_pvalue_corrected'] = fdrcorrection(df['het_p_value'].tolist(), alpha=0.1)
    df['reg_rejected'], df['reg_pvalue_corrected'] = fdrcorrection(df['reg_p_value'].tolist(), alpha=0.1)
    df = df[df['go_mean_het'].notna()]
    return df


def plot_go_het(df):
    """
    Plot the mean heterogeneity score of each GO annotation vs. the genes that are not in the pathway.
    Add star and red color to the significant pathways.
    :param df: the df to plot
    """
    df = df.sort_values(by='go_mean_het')
    idx = df.index.values.tolist()
    fig, ax = plt.subplots()
    fig.set_size_inches(19, 11)
    fig.subplots_adjust(bottom=0.55)
    ax.scatter(np.arange(len(idx)), df["go_mean_het"], s=10, c='darkgreen', vmin=3)
    ax.scatter(np.arange(len(idx)), df["non_go_mean_het"], s=10, vmin=3, c='lightgray', zorder=-1)
    plt.xticks([])
    plt.xlabel("GO annotation", fontsize=13)
    plt.ylabel("Mean heterogeneity score per GO annotation", fontsize=13)
    plt.title("GO annotations heterogeneity scores")
    plt.legend(["Other genes mean heterogeneity score", "Pathway mean heterogeneity score"])
    plt.legend(["GO annotation mean heterogeneity score", "Other genes mean heterogeneity score"])
    plt.savefig("go_het.png")
    plt.show()
    plt.close()


def plot_go_reg(df):
    """
    Plot the mean mROS of each GO annotation vs. the genes that are not in the pathway.
    Add star and red color to the significant pathways.
    :param df: the df to plot
    """
    df = df.sort_values(by='go_mROS')
    idx = df.index.values.tolist()
    fig, ax = plt.subplots()
    fig.set_size_inches(19, 11)
    fig.subplots_adjust(bottom=0.55)
    ax.scatter(np.arange(len(idx)), df["go_mROS"], s=10, c='darkgreen', vmin=3)
    ax.scatter(np.arange(len(idx)), df["non_go_mROS"], s=10, vmin=3, c='lightgray', zorder=-1)
    plt.xticks([])
    plt.xlabel("GO annotation", fontsize=13)
    plt.ylabel("mROS score per GO annotation", fontsize=13)
    plt.title("GO annotations mROS")
    plt.legend(["GO annotation mROS", "Other genes mROS"])
    plt.savefig("go_reg.png")
    plt.show()
    plt.close()


def plot_go_het_partial(df, n):
    """
    Plot the mean heterogeneity score of the highest n GO annotations (by the mean score) vs. the genes
    that are not in the pathway.
    Add star and red color to the significant pathways.
    :param n: number of pathways to plot.
    :param df: the df to plot
    """
    df = df.sort_values(by='go_mean_het')
    df = df.tail(n)
    idx = df.index.values.tolist()
    fig, ax = plt.subplots()
    fig.set_size_inches(19, 11)
    fig.subplots_adjust(bottom=0.55)
    ax.scatter(np.arange(len(idx)), df["go_mean_het"], s=10, c='darkgreen', vmin=3)
    ax.scatter(np.arange(len(idx)), df["non_go_mean_het"], s=10, vmin=3, c='lightgray', zorder=-1)
    df = df.reset_index()
    x_tick_colors = []
    for i, row in df.iterrows():
        corrected_p_value = row['het_pvalue_corrected']
        if float(corrected_p_value) <= 0.1:
            plt.text(i-0.16, int(row['go_mean_het'])+0.5, '*')
            x_tick_colors.append('r')
        else:
            x_tick_colors.append('black')
    ax.xaxis.set_ticks(np.arange(len(idx)))
    ax.set_xticklabels(idx, rotation=45, size=10, ha='right')
    for ticklabel, tickcolor in zip(plt.gca().get_xticklabels(), x_tick_colors):
        ticklabel.set_color(tickcolor)
    plt.xlabel("GO annotation", fontsize=13)
    plt.ylabel("Mean heterogeneity score per GO annotation", fontsize=13)
    plt.title("GO annotations heterogeneity scores")
    plt.legend(["Other genes mean heterogeneity score", "Pathway mean heterogeneity score"])
    plt.legend(["GO annotation mean heterogeneity score", "Other genes mean heterogeneity score"])
    plt.savefig("go_partial_n{0}_het.png".format(n))
    plt.show()
    plt.close()


def plot_go_reg_partial(df, n):
    """
    Plot the mean mROS of the highest n GO annotations (by the mean score) vs. the genes
    that are not in the pathway.
    Add star and red color to the significant pathways.
    :param n: number of pathways to plot.
    :param df: the df to plot
    """
    df = df.sort_values(by='go_mROS')
    head_df = df.head(n)
    idx = head_df.index.values.tolist()
    fig, ax = plt.subplots()
    fig.set_size_inches(19, 11)
    fig.subplots_adjust(bottom=0.55)
    ax.scatter(np.arange(len(idx)), head_df["go_mROS"], s=10, c='darkgreen', vmin=3)
    ax.scatter(np.arange(len(idx)), head_df["non_go_mROS"], s=10, vmin=3, c='lightgray', zorder=-1)
    head_df = head_df.reset_index()
    x_tick_colors = []
    for i, row in head_df.iterrows():
        corrected_p_value = row['reg_pvalue_corrected']
        if float(corrected_p_value) <= 0.1:
            plt.text(i-0.16, row['go_mROS']+0.005, '*')
            x_tick_colors.append('r')
        else:
            x_tick_colors.append('black')
    ax.xaxis.set_ticks(np.arange(len(idx)))
    ax.set_xticklabels(idx, rotation=45, size=10, ha='right')
    for ticklabel, tickcolor in zip(plt.gca().get_xticklabels(), x_tick_colors):
        ticklabel.set_color(tickcolor)
    plt.xlabel("GO annotation", fontsize=13)
    plt.ylabel("mROS score per GO annotation", fontsize=13)
    plt.title("GO annotations mROS - down regulation")
    plt.legend(["Other genes mROS", "Pathway mROS"])
    plt.legend(["GO annotation mROS", "Other genes mROS"])
    plt.savefig("go_partial_n{0}_reg_down.png".format(n))
    plt.show()
    plt.close()

    tail_df = df.tail(n)
    idx = tail_df.index.values.tolist()
    fig, ax = plt.subplots()
    fig.set_size_inches(19, 11)
    fig.subplots_adjust(bottom=0.55)
    ax.scatter(np.arange(len(idx)), tail_df["go_mROS"], s=10, c='darkgreen', vmin=3)
    ax.scatter(np.arange(len(idx)), tail_df["non_go_mROS"], s=10, vmin=3, c='lightgray', zorder=-1)
    tail_df = tail_df.reset_index()
    x_tick_colors = []
    for i, row in tail_df.iterrows():
        corrected_p_value = row['reg_pvalue_corrected']
        if float(corrected_p_value) <= 0.1:
            plt.text(i-0.16, row['go_mROS']+0.005, '*')
            x_tick_colors.append('r')
        else:
            x_tick_colors.append('black')
    ax.xaxis.set_ticks(np.arange(len(idx)))
    ax.set_xticklabels(idx, rotation=45, size=10, ha='right')
    for ticklabel, tickcolor in zip(plt.gca().get_xticklabels(), x_tick_colors):
        ticklabel.set_color(tickcolor)
    plt.xlabel("GO annotation", fontsize=13)
    plt.ylabel("mROS score per GO annotation", fontsize=13)
    plt.title("GO annotations mROS - up regulation")
    plt.legend(["Other genes mROS", "Pathway mROS"])
    plt.legend(["GO annotation mROS", "Other genes mROS"])
    plt.savefig("go_partial_n{0}_reg_up.png".format(n))
    plt.show()
    plt.close()


def main():
    go_dict = go_data_parsing() # go_dict = {k: go_dict[k] for k in go_dict.keys()[:50]}
    df = calc_scores(go_dict)
    df.to_csv("go.csv")
    plot_go_het(df)
    plot_go_reg(df)
    plot_go_het_partial(df, 50)
    plot_go_reg_partial(df, 50)


if __name__ == '__main__':
    main()
