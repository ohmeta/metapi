#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')

import scipy.cluster.hierarchy

import pickle



import drep
import drep.d_cluster
import drep.d_filter
from drep.d_analyze import plot_MASH_dendrogram


import matplotlib
matplotlib.use('Agg')

import os

import pandas as pd
import seaborn as sns
import scipy.cluster.hierarchy

import numpy as np

from matplotlib import pyplot as plt

def _x_fig_size(points, factor= .07, min= 8):
    '''
    Calculate how big the x of the figure should be
    '''
    size = points * factor
    return max([size,min])

def gen_color_dict(names, name2cluster):
    '''
    Make the dictionary name2color

    Args:
        names: key in the returned dictionary
        name2cluster: a dictionary of name to it's cluster

    Returns:
        dict: name -> color
    '''
    #cm = _rand_cmap(len(set(name2cluster.values()))+1,type='bright')
    vals = np.linspace(0,1,len(set(name2cluster.values()))+1)
    np.random.shuffle(vals)
    cm = plt.cm.colors.ListedColormap(plt.cm.jet(vals))

    # 1. generate cluster to color
    cluster2color = {}
    clusters = set(name2cluster.values())
    NUM_COLORS = len(clusters)
    for cluster in clusters:
        try:
            cluster2color[cluster] = cm(1.*int(cluster)/NUM_COLORS)
        except:
            cluster2color[cluster] = cm(1.*int(str(cluster).split('_')[1])/NUM_COLORS)

    #2. name to color
    name2color = {}
    for name in names:
        name2color[name] = cluster2color[name2cluster[name]]

    return name2color

def plot_mash_dendrogram(Mdb, Cdb, linkage, threshold=False, plot_dir=False):
    '''
    Make a dendrogram of the primary clustering

    Args:
        Mdb: DataFrame of Mash comparison results
        Cdb: DataFrame of Clustering results
        linkage: Result of scipy.cluster.hierarchy.linkage
        threshold (optional): Line to plot on x-axis
        plot_dir (optional): Location to store plot

    Returns:
        Makes and shows plot
    '''
    sns.set_style('white',{'axes.grid': False})

    db = Mdb.pivot("genome1","genome2","similarity")
    names = list(db.columns)
    name2cluster = Cdb.set_index('genome')['primary_cluster'].to_dict()
    name2color = gen_color_dict(names, name2cluster)

    # Make the dendrogram
    g = fancy_dendrogram_(linkage,names,name2color,threshold=threshold)
    plt.title('MASH clustering')
    plt.xlabel('MASH Average Nucleotide Identity (ANI)')
    #plt.xlim([0,.4])

    sns.despine(left=True,top=True,right=True,bottom=False)

    # Adjust the figure size
    fig = plt.gcf()
    fig.set_size_inches(10,_x_fig_size(len(names),factor=.2))
    plt.subplots_adjust(left=0.3)

    # Adjust the x labels
    plt.tick_params(axis='both', which='major', labelsize=8)
    axes = plt.gca()
    labels = axes.xaxis.get_majorticklocs()
    for i, label in enumerate(labels):
        labels[i] = (1 - float(label)) * 100
    axes.set_xticklabels(labels)

    # Add cluster to the y axis
    g2c = Cdb.set_index('genome')['secondary_cluster'].to_dict()
    axes = plt.gca()
    labels = [item.get_text() for item in axes.get_yticklabels()]
    for i, label in enumerate(labels):
        labels[i] = "{0} ({1})".format(label, g2c[label])
    axes.set_yticklabels(labels)

    # Save the figure
    if plot_dir != None:
        plt.savefig(os.path.join(plot_dir, 'Primary_clustering_dendrogram_truncate.pdf'),\
            format="pdf", transparent=True, bbox_inches='tight')
    plt.show()
    plt.close('all')

def fancy_dendrogram_(linkage,names,name2color=False,threshold=False,self_thresh=False):
    '''
    Make a fancy dendrogram
    '''
    # Make the dendrogram
    if threshold == False:
        scipy.cluster.hierarchy.dendrogram(linkage,truncate_mode="lastp", labels=names,orientation='right')
    else:
        scipy.cluster.hierarchy.dendrogram(linkage,truncate_mode="lastp", labels=names, color_threshold=threshold,\
                                            orientation='right')

    # Color the names
    if name2color != False:
        ax = plt.gca()
        xlbls = ax.get_ymajorticklabels()
        for lbl in xlbls:
            lbl.set_color(name2color[lbl.get_text()])

    # Add the threshold
    if threshold:
        plt.axvline(x=threshold, c='k', linestyle='dotted')
    if self_thresh != False:
        plt.axvline(x=self_thresh, c='red', linestyle='dotted', lw=1)

    g = plt.gcf()
    return g

def main():
    data_tables = "/home/alienzj/projects/bgimeta_ro/assay/07.drep/output_292S/data_tables"
    clustering_files = "/home/alienzj/projects/bgimeta_ro/assay/07.drep/output_292S/data/Clustering_files"
    bdb_csv = os.path.join(data_tables, "Bdb.csv")
    cdb_csv = os.path.join(data_tables, "Cdb.csv")
    mdb_csv = os.path.join(data_tables, "Mdb.csv")
    ndb_csv = os.path.join(data_tables, "ndb.csv")
    primary_linkage = os.path.join(clustering_files, "primary_linkage.pickle")

    mdb = pd.read_csv(mdb_csv)
    cdb = pd.read_csv(cdb_csv)

    cluster_handle = open(primary_linkage, "rb")
    linkage = pickle.load(cluster_handle)
    linkage_db = pickle.load(cluster_handle)
    arguments = pickle.load(cluster_handle)

    plot_dir = "/home/alienzj/workspace/pycharm/bgimeta_ro/mash_clustering"
    plot_mash_dendrogram(mdb, cdb, linkage, threshold=0.1, plot_dir=plot_dir)

if __name__ == '__main__':
    main()