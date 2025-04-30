# Author: Xiang Zhang
# Date created: 06/09/2021
# Modified: 12/04/2024 for outputing node-specific stats after merging Golgi, endosome and Lysosome 

import sys
import os
from os.path import expanduser
from pathlib import Path

# Add path to folder containing safepy
sys.path.append(expanduser('~') + '/Documents/Research_CB/GIN/safepy_network_annotation/safepy/')

import safe

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

#%matplotlib inline


def get_node_coordinates(graph):

    x = dict(graph.nodes.data('x'))
    y = dict(graph.nodes.data('y'))

    ds = [x, y]
    pos = {}
    for k in x:
        pos[k] = np.array([d[k] for d in ds])

    node_xy = np.vstack(list(pos.values())).astype(float)

    return node_xy

def plot_network_revise1(sf, background_color='#000000', node_list=None, save_fig_path=None):
    ax = plot_network_revise2(sf.graph, background_color=background_color, node_list=node_list, save_fig_path=save_fig_path)
    
def plot_network_revise2(G, ax=None, background_color='#000000', node_list=None, save_fig_path=None):

    foreground_color = '#ffffff'
    if background_color == '#ffffff':
        foreground_color = '#000000'

    node_xy = get_node_coordinates(G)

    if ax is None:
        fig, ax = plt.subplots(figsize=(20, 10), facecolor=background_color, edgecolor=foreground_color)
        fig.set_facecolor(background_color)

    # Randomly sample a fraction of the edges (when network is too big)
    edges = tuple(G.edges())
    if len(edges) > 30000:
        edges = random.sample(edges, int(len(edges)*0.1))

    nx.draw(G, ax=ax, pos=node_xy, edgelist=edges,
            node_color=foreground_color, edge_color=foreground_color, node_size=10, width=1, alpha=0.2)
        
    # Add code for highlighting the genes in the node list
    for i in range(len(G.nodes)):
        if G.nodes[i]['shared name'] in node_list:
            ax.plot(float(G.nodes[i]['x']), float(G.nodes[i]['y']), 'ro', markersize=4)
        
    ax.set_aspect('equal')
    ax.set_facecolor(background_color)

    ax.grid(False)
    ax.invert_yaxis()
    ax.margins(0.1, 0.1)

    ax.set_title('Network', color=foreground_color)
      
    plt.axis('off')

    try:
        fig.set_facecolor(background_color)
    except NameError:
        pass

    plt.savefig(save_fig_path, bbox_inches='tight', dpi=500, facecolor=fig.get_facecolor())
    plt.close(fig)
    
    return ax

def safe_test_plot(input_network_file, input_attribute_file, output_dir, metric='euclidean', nerbor_radius=0.15, dist_thre=0.65, min_neighbor=10, enrichment_threshold=0.05, multiple_testing_bool=True):
    sf = safe.SAFE()
    
    sf.load_network(network_file=input_network_file, node_key_attribute="shared name")

    #sf.plot_network()

    sf.load_attributes(attribute_file=input_attribute_file)
    
    sf.define_neighborhoods(node_distance_metric=metric, neighborhood_radius=nerbor_radius)

    sf.compute_pvalues(multiple_testing=multiple_testing_bool) # May visualize right after this point

    print(sf.pvalues_pos.shape)
    #print(sf.pvalues_pos)
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Assume the order follows the nodes
    for i in range(sf.attributes.shape[0]):
        node_list = []
        sub_node_list = []
        pval_list = []
        nes_list = []
        node_opacity_list = []
        for j in range(len(sf.graph.nodes)):
            pval_ = sf.pvalues_pos[j][i]
            nes_ = sf.nes[j][i]
            alpha_ = -np.log10(sf.pvalues_pos[j][i])/10
            node_list.append(sf.graph.nodes[j]['shared name'])
            pval_list.append(pval_)
            nes_list.append(nes_)
            node_opacity_list.append(alpha_)

            if sf.pvalues_pos[j][i] < enrichment_threshold:
                sub_node_list.append(sf.graph.nodes[j]['shared name'])
        df_i = pd.DataFrame(data={"Gene": node_list, "FDR": pval_list,
                "NES(-log10 FDR)": nes_list, "Opacity(alpha)": node_opacity_list})
        output_dir_pval = output_dir + '/' + sf.attributes.iloc[i]['name'] + "_enrich_stats.csv"
        df_i.to_csv(output_dir_pval, index=None)
        plot_network_revise1(sf, node_list = sub_node_list, save_fig_path=(output_dir + '/' + sf.attributes.iloc[i]['name'] + '.pdf'))
    '''
        for j in range(len(sf.graph.nodes)):
            if sf.pvalues_pos[j][i] < enrichment_threshold:
                sub_node_list.append(sf.graph.nodes[j]['shared name'])
        plot_network_revise1(sf, node_list = sub_node_list, save_fig_path=(output_dir + '/' + sf.attributes.iloc[i]['name'] + '.pdf'))
    '''

    return sf

if __name__ == '__main__':
    #input_file = expanduser('~') + '/Documents/Research_CB/GIN/safepy_network_annotation/network_generation/max_centroid_rm/mc_core_12122022/expanded_c12rm_no507Var2_041_30_min010.cys'
    # 09082023: switch to the updated coordinates (Perox area rotated)
    input_file = "/Users/zhangxiang/Documents/Research_CB/GIN/safepy_network_annotation/network_generation/new_expanded_Dec2022/rotated perox network_new_coordinates.cys"

    #input_attribute_file_compartment_4_5 = '/Users/zhangxiang/Documents/Research_CB/GIN/safepy_network_annotation/human_compartment/annotation_human_compartment_7terms_thre4_5_12042024.txt'
    input_attribute_file_compartment_4_5 = '/Users/zhangxiang/Documents/Research_CB/GIN/safepy_network_annotation/human_compartment/annotation_human_compartment_7terms_thre4_5_EXPANDED_12102024.txt'

    #output_dir_compartment_015_001 = expanduser('~') + '/Documents/Research_CB/GIN/safepy_network_annotation/network_generation/new_expanded_Dec2022/localization/expanded_30_041_c12rm_no507Var2_human_compartment_7terms_thre4_5_radius0_15_enrich0_001_12042024'
    output_dir_compartment_015_001 = expanduser('~') + '/Documents/Research_CB/GIN/safepy_network_annotation/network_generation/new_expanded_Dec2022/localization/expanded_30_041_c12rm_no507Var2_human_compartment_7terms_thre4_5_EXPANDED_radius0_15_enrich0_001_12102024'
    
    sf_compartment_015_005 = safe_test_plot(input_file, input_attribute_file_compartment_4_5, output_dir_compartment_015_001, metric='shortpath_weighted_layout', nerbor_radius=0.15, dist_thre=0.65, min_neighbor=10, enrichment_threshold=0.001)

