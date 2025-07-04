import re
import os

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import scipy.io as spio
import pandas as pd
import zipfile
import random
import shutil
import ast

from os.path import expanduser
from scipy.spatial import ConvexHull
from scipy.spatial.distance import pdist, squareform
from scipy.stats import gaussian_kde
from scipy.optimize import fmin
from collections import Counter
from xml.dom import minidom


def load_network_from_txt(filename, layout='spring_embedded', node_key_attribute='key', verbose=True):
    """
    Loads network from tab-delimited text file and applies a network layout.

    :param str filename: Path to the network file
    :param str layout: Name of the network layout to be applied
    :param str node_key_attribute: Name of the node attribute to be used as key for matching with additional attributes
    :param bool verbose:
    :return:

    """

    filename = re.sub('~', expanduser('~'), filename)

    # Get the number of columns first
    with open(filename, 'r') as f:
        first_line = f.readline()
        num_cols = len(first_line.split('\t'))

    if num_cols == 3:

        data = pd.read_table(filename, sep='\t', header=None, dtype={0: str, 1: str, 2: float})
        data = data.rename(columns={0: 'node_key1', 1: 'node_key2', 2: 'edge_weight'})
        data['node_label1'] = data['node_key1']
        data['node_label2'] = data['node_key2']

    elif num_cols == 5:

        data = pd.read_table(filename, sep='\t', header=None, dtype={0: str, 1: str, 2: str, 3: str, 4: float})
        data = data.rename(
            columns={0: 'node_label1', 1: 'node_key1', 2: 'node_label2', 3: 'node_key2', 4: 'edge_weight'})

    else:

        raise ValueError('Unknown network file format. 3 or 5 columns are expected.')

    # Merge nodes1 and nodes2 and drop duplicates
    nodes = data[['node_label1', 'node_key1']] \
        .append(data[['node_label2', 'node_key2']].rename(columns={'node_label2': 'node_label1', 'node_key2': 'node_key1'})) \
        .drop_duplicates()

    # Re-number the node index
    nodes = nodes.reset_index(drop=True)

    # Add the node index to network data
    nodes = nodes.reset_index().set_index('node_label1')
    data['node_index1'] = nodes.loc[data['node_label1'], 'index'].values
    data['node_index2'] = nodes.loc[data['node_label2'], 'index'].values

    # Create the graph
    G = nx.Graph()

    # Add the nodes & their attributes
    nodes = nodes.reset_index().set_index('index')
    G.add_nodes_from(nodes.index.values)

    for n in G:
        G.nodes[n]['label'] = nodes.loc[n, 'node_label1']
        G.nodes[n][node_key_attribute] = nodes.loc[n, 'node_key1']

    # Add the edges between the nodes
    edges = [tuple(x) for x in data[['node_index1', 'node_index2']].values]
    G.add_edges_from(edges)

    G = apply_network_layout(G, layout=layout)
    G = calculate_edge_lengths(G, verbose=verbose)

    return G


def load_network_from_gpickle(filename, verbose=True):

    filename = re.sub('~', expanduser('~'), filename)
    G = nx.read_gpickle(filename)

    return G


def load_network_from_mat(filename, verbose=True):

    filename = re.sub('~', expanduser('~'), filename)

    if verbose:
        print('Loading the mat file...')

    mat = load_mat(filename)
    G = nx.Graph(mat['layout']['edges'])

    for n in G:
        G.nodes[n]['label'] = mat['layout']['label'][n]
        G.nodes[n]['key'] = mat['layout']['label_orf'][n]
        G.nodes[n]['x'] = mat['layout']['x'][n]
        G.nodes[n]['y'] = mat['layout']['y'][n]

    G = calculate_edge_lengths(G, verbose=verbose)

    return G


def load_network_from_cys(filename, verbose=True):

    filename = re.sub('~', expanduser('~'), filename)

    if verbose:
        print('Loading the cys file %s...' % filename)

    # Unzip CYS file
    zip_ref = zipfile.ZipFile(filename, 'r')
    files = zip_ref.namelist()
    zip_ref.extractall('./')
    zip_ref.close()

    # Keep the name of the top directory (to remove the unzipped files after we're done)
    top_dirs = list(set([f.split('/')[0] for f in files]))

    # Get node positions (from the view)
    viewfile = [f for f in files if '/views/' in f][0]

    if verbose:
        print('Loading the first view: %s' % viewfile)

    mydoc = minidom.parse(viewfile)
    nodes = mydoc.getElementsByTagName('node')

    node_labels = dict()
    node_xs = dict()
    node_ys = dict()

    for node in nodes:
        node_id = int(node.attributes['cy:nodeId'].value)
        node_labels[node_id] = node.attributes['label'].value
        for child in node.childNodes:
            if (child.nodeType == 1) and (child.tagName == 'graphics'):
                node_xs[node_id] = float(child.attributes['x'].value)
                node_ys[node_id] = float(child.attributes['y'].value)

    # Get edges (from the network)
    networkfile = [f for f in files if '/networks/' in f][0]

    if verbose:
        print('Loading the first network: %s' % networkfile)

    mydoc = minidom.parse(networkfile)
    edges = mydoc.getElementsByTagName('edge') # Xiang: this tag does not have any info about pcc values. Need to switch to the "table" in CytoscapeSession
    
    # Added by Xiang on 12/16/2021 - Read the edge attributes (from /tables/)
    [file_name, file_extension] = os.path.splitext(os.path.basename(networkfile))
    contains = ['/tables/', file_name, 'SHARED_ATTRS', 'edge.cytable']
    attributefile = [f for f in files if all(c in f for c in contains)]

    attributes = pd.read_csv(attributefile[0], sep=',', header=None, skiprows=1)
    
    col_headers = []
    row_start = 0
    for ix_row in np.arange(7):
        val = attributes.iloc[ix_row, 0]
        if val == 'SUID':
            col_headers = list(attributes.iloc[ix_row, :])
        elif str(val).isnumeric():
            row_start = ix_row
            break

    attributes.columns = col_headers
    attributes = attributes.iloc[row_start:, :]

    attributes['SUID'] = attributes['SUID'].astype(int)
    
#print(attributes[attributes['shared name']==edges.item(2).getAttribute('label')]['value'].values[0])

    # Integrate the weight into the edges
    edge_list = []
    for edge in edges:
        e1 = int(edge.attributes['source'].value)
        e2 = int(edge.attributes['target'].value)
        e3 = float(attributes[attributes['shared name']==edge.getAttribute('label')]['weight'].values[0])
        edge_list.append((e1, e2, e3))
    
    # Build the graph
    G = nx.Graph()
    G.add_weighted_edges_from(edge_list) # changed from add_edges_from by Xiang on 12/16/2021
	
    nodes_to_remove = []
    for node in G.nodes:
        if node in node_labels.keys():
            G.nodes[node]['label'] = node_labels[node]
            G.nodes[node]['x'] = float(node_xs[node])
            G.nodes[node]['y'] = float(node_ys[node])
        else:
            nodes_to_remove.append(node)

    for node in nodes_to_remove:
        G.remove_node(node)

    # Read the node attributes (from /tables/)
    [file_name, file_extension] = os.path.splitext(os.path.basename(networkfile))
    contains = ['/tables/', file_name, 'SHARED_ATTRS', 'node.cytable']
    attributefile = [f for f in files if all(c in f for c in contains)]

    attributes = pd.read_csv(attributefile[0], sep=',', header=None, skiprows=1)

    col_headers = []
    row_start = 0
    for ix_row in np.arange(7):
        val = attributes.iloc[ix_row, 0]
        if val == 'SUID':
            col_headers = list(attributes.iloc[ix_row, :])
        elif str(val).isnumeric():
            row_start = ix_row
            break

    attributes.columns = col_headers
    attributes = attributes.iloc[row_start:, :]

    attributes['SUID'] = attributes['SUID'].astype(int)

    for ix_row, row in attributes.iterrows():
        if row['SUID'] in G.nodes:
            for c in col_headers[1:]:
                G.nodes[row['SUID']][c] = row[c]

    # Relabel the node ids to sequential numbers to make calculations faster
    mapping = dict()
    for ix_node, node in enumerate(G.nodes):
        mapping[node] = ix_node

    G = nx.relabel_nodes(G, mapping)

#G = apply_network_layout(G, layout='spring_embedded', verbose=True) # Added by Xiang on 12/13/2021 for testing the spring embedded layout. Uncomment this if you are going to apply the spring embedded layout

    G = calculate_edge_lengths(G, verbose=verbose)

    # Remove unzipped files/directories
    for top_dir in top_dirs:
        shutil.rmtree(top_dir)

    return G


def apply_network_layout(G, layout='spring_embedded', verbose=True):

    if layout == 'kamada_kawai':

        if verbose:
            print('Applying the Kamada-Kawai network layout... (may take several minutes)')

        pos = nx.kamada_kawai_layout(G)

    elif layout == 'spring_embedded':

        if verbose:
            print('Applying the spring-embedded network layout... (may take several minutes)')

        #pos = nx.spring_layout(G, k=0.2, iterations=100) # modified by Xiang on 12/13/2021
        pos = nx.spring_layout(G, k=0.2, iterations=100, seed=42)

    for n in G:
        G.nodes[n]['x'] = pos[n][0]
        G.nodes[n]['y'] = pos[n][1]

    return G


def calculate_edge_lengths(G, verbose=True):

    # Calculate the lengths of the edges

    if verbose:
        print('Calculating edge lengths...')

    x = np.matrix(G.nodes.data('x'))[:, 1]
    y = np.matrix(G.nodes.data('y'))[:, 1]

    node_coordinates = np.concatenate([x, y], axis=1)
    node_distances = squareform(pdist(node_coordinates, 'euclidean'))
    #print(node_distances) # added by Xiang

    adjacency_matrix = np.array(nx.adjacency_matrix(G).todense())
    adjacency_matrix = adjacency_matrix.astype(float)
    adjacency_matrix[adjacency_matrix == 0] = np.nan

    edge_lengths = np.multiply(node_distances, adjacency_matrix)
    
    edge_attr_dict = {index: v for index, v in np.ndenumerate(edge_lengths) if ~np.isnan(v)}
    nx.set_edge_attributes(G, edge_attr_dict, 'length')

    return G


def load_attributes(attribute_file='', node_label_order=None, mask_duplicates=False, fill_value=np.nan, verbose=True):

    node2attribute = pd.DataFrame()
    attributes = pd.DataFrame()

    if isinstance(attribute_file, str):

        file_name = re.sub('~', expanduser('~'), attribute_file)
        [_, file_extension] = os.path.splitext(file_name)

        if file_extension == '.mat':
            mat = load_mat(file_name)

            node2attribute = pd.DataFrame(data=np.transpose(mat['go']['term2orf']),
                                          index=mat['go']['orfs'],
                                          columns=mat['go']['term_ids'])
            node2attribute = node2attribute.apply(pd.to_numeric, downcast='unsigned')

            data = {'id': mat['go']['term_ids'], 'name': mat['go']['term_names']}
            attributes = pd.DataFrame(data=data)

        elif (file_extension == '.txt') or (file_extension == '.gz'):

            node2attribute = pd.read_csv(file_name, sep='\t', dtype={0: str})
            node2attribute.set_index(node2attribute.columns[0], drop=True, inplace=True)
            node2attribute = node2attribute.apply(pd.to_numeric, downcast='float', errors='coerce')

            data = {'id': np.arange(len(node2attribute.columns)), 'name': node2attribute.columns}
            attributes = pd.DataFrame(data=data)

            node2attribute.columns = np.arange(len(node2attribute.columns))

    elif isinstance(attribute_file, pd.DataFrame):

        node2attribute = attribute_file
        data = {'id': np.arange(len(node2attribute.columns)), 'name': node2attribute.columns}
        attributes = pd.DataFrame(data=data)

    # Force all values to numeric
    node2attribute = node2attribute.apply(pd.to_numeric, errors='coerce')

    # Force attribute names to be strings
    attributes['name'] = attributes['name'].astype(str)

    # Averaging duplicate rows (with notification)
    if not node2attribute.index.is_unique:
        print('\nThe attribute file contains multiple values for the same labels. Their values will be averaged.')
        node2attribute = node2attribute.groupby(node2attribute.index, axis=0).mean()
        
    if not node_label_order:
        node_label_order = node2attribute.index.values

    node_label_in_file = node2attribute.index.values
    node_label_not_mapped = [x for x in node_label_in_file if x not in node_label_order]

    node2attribute = node2attribute.reindex(index=node_label_order, fill_value=fill_value)

    # If necessary, de-duplicate the network nodes (leave one node per gene)
    if mask_duplicates:

        # Keep a random node every time
        idx = np.random.permutation(np.arange(len(node2attribute)))
        mask_dups = node2attribute.iloc[idx].index.duplicated(keep='first')

        num_dups = mask_dups.sum()
        print('\nThe network contains %d nodes with duplicate labels. '
              'Only one random node per label will be considered. '
              'The attribute values of all other nodes will be set to NaN.' % num_dups)
        node2attribute.iloc[idx[mask_dups], :] = np.nan

    node2attribute = node2attribute.values

    if verbose:
        print('\nAttribute data provided: %d labels x %d attributes' % (len(node_label_in_file), attributes.shape[0]))

        # Notification about labels **not** mapped onto the network
        n = np.min([len(node_label_not_mapped), 3])
        m = len(node_label_not_mapped) - n
        if n > 0:
            msg1 = ', '.join(node_label_not_mapped[:n])
            msg2 = format(' and %d other labels in the attribute file were not found in the network.' % m)
            print(msg1 + msg2)

        n_nlm = len(node_label_in_file) - len(node_label_not_mapped)
        print('\nAttribute data mapped onto the network: %d labels x %d attributes' % (n_nlm, attributes.shape[0]))
        print('Values: %d NaNs' % np.sum(np.isnan(node2attribute)))
        print('Values: %d zeros' % np.sum(node2attribute[~np.isnan(node2attribute)] == 0))
        print('Values: %d positives' % np.sum(node2attribute[~np.isnan(node2attribute)] > 0))
        print('Values: %d negatives' % np.sum(node2attribute[~np.isnan(node2attribute)] < 0))

    return attributes, node_label_order, node2attribute


def plot_network(G, ax=None, background_color='#000000'):

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

    return ax


def plot_network_contour(graph, ax, background_color='#000000'):

    foreground_color = '#ffffff'
    if background_color == '#ffffff':
        foreground_color = '#000000'

    x = dict(graph.nodes.data('x'))
    y = dict(graph.nodes.data('y'))

    ds = [x, y]
    pos = {}
    for k in x:
        pos[k] = np.array([d[k] for d in ds]).astype('float64')

    # Compute the convex hull to delineate the network
    hull = ConvexHull(np.vstack(list(pos.values())))

    vertices_x = [pos.get(v)[0] for v in hull.vertices]
    vertices_y = [pos.get(v)[1] for v in hull.vertices]

    vertices_x = np.array(vertices_x).astype('float64')
    vertices_y = np.array(vertices_y).astype('float64')

    # Find center of mass and radius to approximate the hull with a circle
    xm = np.nanmean(vertices_x)
    ym = np.nanmean(vertices_y)

    rm = np.nanmean(np.sqrt((vertices_x - xm) ** 2 + (vertices_y - ym) ** 2))

    # Best fit a circle to these points
    def err(x0):
        [w, v, r] = x0
        pts = [np.linalg.norm([x - w, y - v]) - r for x, y in zip(vertices_x, vertices_y)]
        return (np.array(pts) ** 2).sum()

    [xf, yf, rf] = fmin(err, [xm, ym, rm], disp=False)

    circ = plt.Circle((xf, yf), radius=rf*1.01, color=foreground_color, linewidth=1, fill=False)
    ax.add_patch(circ)

    return xf, yf, rf


def plot_costanzo2016_network_annotations(graph, ax, path_to_data, colors=True, clabels=False,
                                          background_color='#000000'):

    foreground_color = '#ffffff'
    if background_color == '#ffffff':
        foreground_color = '#000000'

    path_to_network_annotations = 'other/Data File S5_SAFE analysis_Gene cluster identity and functional enrichments.xlsx'
    filename = os.path.join(path_to_data, path_to_network_annotations)

    costanzo2016 = pd.read_excel(filename, sheet_name='Global net. cluster gene list')
    processes = costanzo2016['Global Similarity Network Region name'].unique()
    processes = processes[pd.notnull(processes)]

    process_colors = pd.read_csv(os.path.join(path_to_data, 'other/costanzo_2016_colors.txt'), sep='\t')
    if colors:
        process_colors = process_colors[['R', 'G', 'B']].values / 256
    else:
        if foreground_color == '#ffffff':
            process_colors = np.ones((process_colors.shape[0], 3))
        else:
            process_colors = np.zeros((process_colors.shape[0], 3))

    labels = nx.get_node_attributes(graph, 'label')
    labels_dict = {k: v for v, k in labels.items()}

    x = list(dict(graph.nodes.data('x')).values())
    y = list(dict(graph.nodes.data('y')).values())

    pos = {}
    for idx, k in enumerate(x):
        pos[idx] = np.array([x[idx], y[idx]])

    for n_process, process in enumerate(processes):
        nodes = costanzo2016.loc[costanzo2016['Global Similarity Network Region name'] == process, 'Gene/Allele']
        nodes_indices = [labels_dict[node] for node in nodes if node in labels_dict.keys()]

        pos3 = {idx: pos[node_index] for idx, node_index in enumerate(nodes_indices)}
        pos3 = np.vstack(list(pos3.values()))

        kernel = gaussian_kde(pos3.T)
        [X, Y] = np.mgrid[np.min(x):np.max(x):100j, np.min(y):np.max(y):100j]
        positions = np.vstack([X.ravel(), Y.ravel()])
        Z = np.reshape(kernel(positions).T, X.shape)

        C = ax.contour(X, Y, Z, [1e-6], colors=[tuple(process_colors[n_process, :])], alpha=1)

        if clabels:
            C.levels = [n_process+1]
            plt.clabel(C, C.levels, inline=True, fmt='%d', fontsize=16)
            print('%d -- %s' % (n_process+1, process))


def plot_labels(labels, graph, ax):

    node_labels = nx.get_node_attributes(graph, 'label')
    node_labels_dict = {k: v for v, k in node_labels.items()}

    x = list(dict(graph.nodes.data('x')).values())
    y = list(dict(graph.nodes.data('y')).values())

    # x_offset = (np.nanmax(x) - np.nanmin(x))*0.01

    idx = [node_labels_dict[x] for x in labels if x in node_labels_dict.keys()]
    labels_idx = [x for x in labels if x in node_labels_dict.keys()]
    x_idx = [x[i] for i in idx]
    y_idx = [y[i] for i in idx]

    # ax.plot(x_idx, y_idx, 'r*')
    for i in np.arange(len(idx)):
        ax.text(x_idx[i], y_idx[i], labels_idx[i], fontdict={'color': 'white', 'size': 14, 'weight': 'bold'},
                bbox={'facecolor': 'black', 'alpha': 0.5, 'pad': 3},
                horizontalalignment='center', verticalalignment='center')

    # Print out labels not found
    labels_missing = [x for x in labels if x not in node_labels_dict.keys()]
    if labels_missing:
        labels_missing_str = ', '.join(labels_missing)
        print('These labels are missing from the network (case sensitive): %s' % labels_missing_str)


def get_node_coordinates(graph):

    x = dict(graph.nodes.data('x'))
    y = dict(graph.nodes.data('y'))

    ds = [x, y]
    pos = {}
    for k in x:
        pos[k] = np.array([d[k] for d in ds])

    node_xy = np.vstack(list(pos.values())).astype(float)

    return node_xy


def load_mat(filename):
    """
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    """

    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)


def _check_keys(dict):

    """
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    """

    for key in dict:
        if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict


def _todict(matobj):

    """
    A recursive function which constructs from matobjects nested dictionaries
    """

    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict


def chop_and_filter(s):
    single_str = s.str.cat(sep=' ')
    single_list = re.findall(r"[\w']+", single_str)

    single_list_count = dict(Counter(single_list))
    single_list_count = [k for k in sorted(single_list_count, key=single_list_count.get, reverse=True)]

    to_exclude = ['of', 'a', 'the', 'an', ',', 'via', 'to', 'into', 'from']
    single_list_words = [w for w in single_list_count if w not in to_exclude]

    return ', '.join(single_list_words[:5])




