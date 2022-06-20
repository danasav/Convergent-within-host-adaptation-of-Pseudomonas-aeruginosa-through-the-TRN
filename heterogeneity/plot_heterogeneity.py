"""
Create the graph of the heterogeneity as function of the weight per each gene
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from config import *


def plot_heterogeneity_vs_mROS():
    """
    Create the graph of the heterogeneity as function of the weight per each gene
    """
    df = pd.read_csv(SUPPLEMENTARY_TABLE_1)
    df["color"] = "gray"
    df.loc[(df["Heterogeneity"] > 2.5) & (df["Excluded Weight"] > 0.15), "color"] = 'lightcoral'
    df.loc[(df["Heterogeneity"] > 2.5) & (df["Excluded Weight"] < -0.15),"color"] = 'royalblue'
    weight_jitter = rand_jitter(df["Excluded Weight"].values.tolist())
    heterogeneity_jitter = rand_jitter(df["Heterogeneity"].values.tolist())
    # plt.scatter(df["Excluded Weight"].values.tolist(), df["Heterogeneity"].values.tolist(), color=df["color"], s=6)
    plt.scatter(weight_jitter, heterogeneity_jitter, color=df["color"], s=6)
    plt.title("Heterogeneiy of the regulatory network")
    plt.xlabel("Regulatory outcome score", fontsize=13)
    plt.ylabel("Heterogeneity score", fontsize=13)
    plt.savefig("Heterogeneity_of_the_regulatory_network.png")
    plt.show()
    plt.close()


def jitter(x, y, s=20, c='b', marker='o', cmap=None, norm=None, vmin=None, vmax=None, alpha=None, linewidths=None,
           verts=None, hold=None, **kwargs):
    """
    Add random jitter in order to avoid overlapping data points.
    """
    return plt.scatter(rand_jitter(x), rand_jitter(y), s=s, c=c, marker=marker, cmap=cmap, norm=norm, vmin=vmin,
                       vmax=vmax, alpha=alpha, linewidths=linewidths, verts=verts, hold=hold, **kwargs)


def rand_jitter(arr):
    """
    Add the random jitter to the data array.
    """
    stdev = .003 * (max(arr) - min(arr))
    return arr + np.random.randn(len(arr)) * stdev


def main():
    plot_heterogeneity_vs_mROS()


if __name__ == '__main__':
    main()
