import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import convert_network_to_refseq
from config import *
from statsmodels.stats.multitest import fdrcorrection


def kegg_data_parsing():
    """
    Parse the KEGG data file, the  and the gene convert names file
    @return: kegg_dict: a dictionary of pathways and their gene numbers, kegg_dict_with_names: a dictionary of pathways
    and their gene names, genes_in_pathways: a set of all the genes that take part in these pathways.
    """
    gene_convert_dict = convert_network_to_refseq.parse_cds(CDS_FILE_P, EXTRA_FILES)
    kegg_convert_dict = kegg_converter()
    kegg_dict = {}
    data = open(KEGG_PATH, "r")
    for row in data.readlines():
        row = row.split()
        pathway = (row[0]).split(":")[1][3:]
        if pathway in kegg_convert_dict:
            pathway = kegg_convert_dict[pathway]# + " (" + pathway + ")"
        gene = (row[1]).split(":")[1].lower()
        if gene in gene_convert_dict:
            if pathway in kegg_dict:
                kegg_dict[pathway].append(gene_convert_dict[gene])
            else:
                kegg_dict[pathway] = [gene_convert_dict[gene]]
    return kegg_dict


def kegg_converter():
    """
    Convert the kegg numbers to kegg descriptions
    @return:the converter dictionary
    """
    convert_data = open(CONVERT_KEGG_PATH, "r")
    kegg_convert_dict = {}
    for row in convert_data.readlines():
        kegg_convert_dict[row.split()[0]] = " ".join(row.split()[1:])
    return kegg_convert_dict


def calc_scores(kegg_dict):
    """
    Calc the scores (mean heterogeneity, mean mROS) of each KEGG pathway (and the genes that are not in the pathway).
    @return: the updated data frame
    """
    df = pd.DataFrame(columns=kegg_dict.keys()).T
    pao1_df = pd.read_csv(SUPPLEMENTARY_TABLE_1)
    for kegg in kegg_dict:
        kegg_het_scores = []
        non_kegg_het_scores = []
        kegg_reg_scores = []
        non_kegg_reg_scores = []
        for i, gene in enumerate(pao1_df['Node name']):
            if gene in kegg_dict[kegg]:
                if pao1_df.loc[i, 'Weight'] != 0:
                    kegg_het_scores.append(pao1_df.loc[i, 'Heterogeneity'])
                if str(pao1_df.loc[i, 'Weight']) != 'nan':
                    kegg_reg_scores.append(pao1_df.loc[i, 'Weight'])
            else:
                if pao1_df.loc[i, 'Weight'] != 0:
                    non_kegg_het_scores.append(pao1_df.loc[i, 'Heterogeneity'])
                if str(pao1_df.loc[i, 'Weight']) != 'nan':
                    non_kegg_reg_scores.append(pao1_df.loc[i, 'Weight'])
        df.loc[kegg, 'kegg_mean_het'] = np.mean(kegg_het_scores)
        df.loc[kegg, 'non_kegg_mean_het'] = np.mean(non_kegg_het_scores)
        df.loc[kegg, 'num_of_genes'] = len(kegg_het_scores)
        df.loc[kegg, 'kegg_mROS'] = np.mean(kegg_reg_scores)
        df.loc[kegg, 'non_kegg_mROS'] = np.mean(non_kegg_reg_scores)
        _, het_p_value = mannwhitneyu(kegg_het_scores, non_kegg_het_scores, alternative='greater')
        _, reg_p_value = mannwhitneyu(kegg_reg_scores, non_kegg_reg_scores, alternative='two-sided')
        df.loc[kegg, 'het_p_value'] = het_p_value
        df.loc[kegg, 'reg_p_value'] = reg_p_value
    df = df[df["num_of_genes"] > 1]
    df['het_rejected'], df['het_pvalue_corrected'] = fdrcorrection(df['het_p_value'].tolist(), alpha=0.1)
    df['reg_rejected'], df['reg_pvalue_corrected'] = fdrcorrection(df['reg_p_value'].tolist(), alpha=0.1)
    df = df[df['kegg_mean_het'].notna()]
    return df

def plot_kegg_het(df):
    """
    Plot the mean heterogeneity score of each KEGG pathway vs. the genes that are not in the pathway.
    Add star and red color to the significant pathways.
    :param df: the df to plot
    """
    df = df.sort_values(by='kegg_mean_het')
    idx = df.index.values.tolist()
    fig, ax = plt.subplots()
    fig.set_size_inches(19, 11)
    fig.subplots_adjust(bottom=0.55)
    ax.scatter(np.arange(len(idx)), df["kegg_mean_het"], s=10, c='darkgreen', vmin=3)
    ax.scatter(np.arange(len(idx)), df["non_kegg_mean_het"], s=10, vmin=3, c='lightgray', zorder=-1)
    df = df.reset_index()
    for i, row in df.iterrows():
        corrected_p_value = row['het_pvalue_corrected']
        if float(corrected_p_value) <= 0.1:
            plt.text(i-0.28, row['kegg_mean_het']+0.005, '*')
    plt.xticks([])
    plt.xlabel("Pathway", fontsize=13)
    plt.ylabel("Mean heterogeneity score per KEGG pathway", fontsize=13)
    plt.title("KEGG pathways heterogeneity scores")
    plt.legend(["Other genes mean heterogeneity score", "Pathway mean heterogeneity score"])
    plt.legend(["Pathway mean heterogeneity score", "Other genes mean heterogeneity score"])
    plt.savefig("kegg_het.png")
    plt.show()
    plt.close()


def plot_kegg_reg(df):
    """
    Plot the mean mROS of each KEGG pathway vs. the genes that are not in the pathway.
    Add star and red color to the significant pathways.
    :param df: the df to plot
    """
    df = df.sort_values(by='kegg_mROS')
    idx = df.index.values.tolist()
    fig, ax = plt.subplots()
    fig.set_size_inches(19, 11)
    fig.subplots_adjust(bottom=0.55)
    ax.scatter(np.arange(len(idx)), df["kegg_mROS"], s=10, c='darkgreen', vmin=3)
    ax.scatter(np.arange(len(idx)), df["non_kegg_mROS"], s=10, vmin=3, c='lightgray', zorder=-1)
    df = df.reset_index()
    for i, row in df.iterrows():
        corrected_p_value = row['reg_pvalue_corrected']
        if float(corrected_p_value) <= 0.1:
            plt.text(i-0.28, row['kegg_mROS']+0.001, '*')
    plt.xticks([])
    plt.xlabel("Pathway", fontsize=13)
    plt.ylabel("mROS score per KEGG pathway", fontsize=13)
    plt.title("KEGG pathways mROS")
    plt.legend(["Other genes mROS", "Pathway mROS"])
    plt.legend(["Pathway mROS", "Other genes mROS"])
    plt.savefig("kegg_reg.png")
    plt.show()
    plt.close()


def plot_kegg_het_partial(df, n=50):
    """
    Plot the mean heterogeneity score of the highest n KEGG pathways (by the mean score) vs. the genes
    that are not in the pathway.
    Add star and red color to the significant pathways.
    :param n: number of pathways to plot.
    :param df: the df to plot
    """
    df = df.sort_values(by='kegg_mean_het')
    df = df.tail(n)
    idx = df.index.values.tolist()
    fig, ax = plt.subplots()
    fig.set_size_inches(19, 6.5)
    fig.subplots_adjust(bottom=0.55)
    ax.scatter(np.arange(len(idx)), df["kegg_mean_het"], s=10, c='darkgreen', vmin=3)
    ax.scatter(np.arange(len(idx)), df["non_kegg_mean_het"], s=10, vmin=3, c='lightgray', zorder=-1)
    df = df.reset_index()
    x_tick_colors = []
    for i, row in df.iterrows():
        corrected_p_value = row['het_pvalue_corrected']
        if float(corrected_p_value) <= 0.1:
            plt.text(i-0.14, row['kegg_mean_het']+0.005, '*')
            x_tick_colors.append('r')
        else:
            x_tick_colors.append('black')
    ax.xaxis.set_ticks(np.arange(len(idx)))
    ax.xaxis.set_ticklabels(idx, rotation=45, size=10, ha='right')
    for ticklabel, tickcolor in zip(plt.gca().get_xticklabels(), x_tick_colors):
        ticklabel.set_color(tickcolor)
    plt.xlabel("Pathway", fontsize=13)
    plt.ylabel("Mean heterogeneity score", fontsize=13)
    plt.title("KEGG pathways heterogeneity scores")
    plt.legend(["Other genes mean heterogeneity score", "Pathway mean heterogeneity score"])
    plt.legend(["Pathway mean heterogeneity score", "Other genes mean heterogeneity score"])
    plt.savefig("kegg_partial_n{0}_het.png".format(n))
    plt.show()
    plt.close()


def plot_kegg_reg_partial(df, n):
    """
    Plot the mean heterogeneity score of the highest n KEGG pathways (by the mean score) vs. the genes
    that are not in the pathway.
    Add star and red color to the significant pathways.
    :param n: number of pathways to plot.
    :param df: the df to plot
    """
    df = df.sort_values(by='kegg_mROS')
    head_df = df.head(n)
    idx = head_df.index.values.tolist()
    fig, ax = plt.subplots()
    fig.set_size_inches(10.2, 7)
    fig.subplots_adjust(bottom=0.55)
    ax.scatter(np.arange(len(idx)), head_df["kegg_mROS"], s=10, c='darkgreen', vmin=3)
    ax.scatter(np.arange(len(idx)), head_df["non_kegg_mROS"], s=10, vmin=3, c='lightgray', zorder=-1)
    head_df = head_df.reset_index()
    x_tick_colors = []
    for i, row in head_df.iterrows():
        corrected_p_value = row['reg_pvalue_corrected']
        if float(corrected_p_value) <= 0.1:
            plt.text(i-0.14, row['kegg_mROS']+0.0005, '*')
            x_tick_colors.append('r')
        else:
            x_tick_colors.append('black')
    ax.xaxis.set_ticks(np.arange(len(idx)))
    ax.xaxis.set_ticklabels(idx, rotation=45, size=10, ha='right')
    for ticklabel, tickcolor in zip(plt.gca().get_xticklabels(), x_tick_colors):
        ticklabel.set_color(tickcolor)
    plt.xlabel("Pathway", fontsize=6.5)
    plt.ylabel("Mean mROS", fontsize=13)
    plt.title("KEGG pathways mROS - down regulation")
    plt.legend(["Other genes mROS", "Pathway mROS"])
    plt.legend(["Pathway mROS", "Other genes mROS"])
    fig.tight_layout(w_pad=2)
    plt.savefig("kegg_partial_n{0}_reg_down.png".format(n))
    plt.show()
    plt.close()

    tail_df = df.tail(n)
    idx = tail_df.index.values.tolist()
    fig, ax = plt.subplots()
    fig.set_size_inches(9.5, 6.5)
    fig.subplots_adjust(bottom=0.55)
    ax.scatter(np.arange(len(idx)), tail_df["kegg_mROS"], s=10, c='darkgreen', vmin=3)
    ax.scatter(np.arange(len(idx)), tail_df["non_kegg_mROS"], s=10, vmin=3, c='lightgray', zorder=-1)
    tail_df = tail_df.reset_index()
    x_tick_colors = []
    for i, row in tail_df.iterrows():
        corrected_p_value = row['reg_pvalue_corrected']
        if float(corrected_p_value) <= 0.1:
            plt.text(i-0.14, row['kegg_mROS']+0.0005, '*')
            x_tick_colors.append('r')
        else:
            x_tick_colors.append('black')
    ax.xaxis.set_ticks(np.arange(len(idx)))
    ax.xaxis.set_ticklabels(idx, rotation=45, size=10, ha='right')
    for ticklabel, tickcolor in zip(plt.gca().get_xticklabels(), x_tick_colors):
        ticklabel.set_color(tickcolor)
    plt.xlabel("Pathway", fontsize=13)
    plt.ylabel("Mean mROS", fontsize=13)
    plt.title("KEGG pathways mROS - up regulation")
    plt.legend(["Other genes mROS", "Pathway mROS"])
    plt.legend(["Pathway mROS", "Other genes mROS"])
    fig.tight_layout()
    plt.savefig("kegg_partial_n{0}_reg_up.png".format(n))
    plt.show()
    plt.close()


def main():
    kegg_dict = kegg_data_parsing()
    df = calc_scores(kegg_dict)
    df.to_csv("kegg.csv")
    plot_kegg_het(df)
    plot_kegg_reg(df)
    plot_kegg_het_partial(df, 40)
    plot_kegg_reg_partial(df, 20)


if __name__ == '__main__':
    main()
