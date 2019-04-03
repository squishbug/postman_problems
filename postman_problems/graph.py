import warnings
import networkx as nx
import pandas as pd
import numpy as np
import itertools
from math import inf

def read_edgelist(edgelist_filename, keep_optional=False):
    """
    Read an edgelist table into a pandas dataframe
    Args:
        edgelist_filename (str): filename of edgelist.  See cpp.py for more details.
        keep_optional (Boolean): keep or discard optional edges (used for RPP)

    Returns:
        pandas dataframe of edgelist
    """
    el = pd.read_csv(edgelist_filename, dtype={0: str, 1: str})  # node_ids as strings makes life easier
    el = el.dropna(how='all')  # drop rows with all NAs... as I find CSVs created w Numbers annoyingly do.

    if (not keep_optional) & ('required' in el.columns):
        el = el[el['required'] == 1]

    assert 'augmented' not in el.columns, \
        'Edgelist cannot contain a column named "augmented", sorry. This will cause computation problems'

    if 'id' in el.columns:
        warnings.warn("Edgelist contains field named 'id'.  This is a field that will be assigned to edge attributes "
                      "with the `create_networkx_graph_from_edgelist function.  That is OK though.  We'll use your 'id'"
                      "field if it is unique.")
        assert el['id'].nunique() == len(el), 'Provided edge "id" field is not unique.  Please drop "id" or try again.'
    return el


def create_networkx_graph_from_edgelist(edgelist, edge_id='id'):
    """
    Create a networkx MultiGraph object from an edgelist (pandas dataframe).
    Used to create the user's starting graph for which a CPP solution is desired.

    Args:
        edgelist (pandas dataframe): output of `read_edgelist` function.
            The first two columns are treated as source and target node names.
            The following columns are treated as edge attributes.
        edge_id (str): name of edge attribute which will be used in `create_eulerian_circuit`.

    Returns:
        networkx.MultiGraph:
            Returning a MultiGraph rather than Graph to support parallel edges
    """
    g = nx.MultiGraph()
    if edge_id in edgelist.columns:
        warnings.warn('{} is already an edge attribute in `edgelist`.  We will try to use it, but recommend '
                      'renaming this column in your edgelist to allow this function to create it in a standardized way'
                      'where it is guaranteed to be unique'.format(edge_id))

    for i, row in enumerate(edgelist.iterrows()):
        edge_attr_dict = row[1][2:].to_dict()
        if edge_id not in edge_attr_dict:
            edge_attr_dict[edge_id] = i
        g.add_edge(row[1][0], row[1][1], **edge_attr_dict)
    return g


def read_graphml(filename, edge_weight="distance", max_degree_connect=None):
    g_full = nx.read_graphml(filename)
    g_full = nx.MultiGraph(g_full) # convert Graph to MultiGraph (adds "keys")

    # every edge is "required"
    for n1,n2,k in g_full.edges(keys=True):
        g_full.edges[(n1,n2,k)]["required"] = 1
        g_full.edges[(n1,n2,k)][edge_weight] = float(g_full.edges[(n1,n2,k)][edge_weight])

    # add optional edges
    if max_degree_connect is not None:
        # get nodes with degree < max_degree_connect
        if max_degree_connect==-1:
            nodes = list(g_full.nodes())
        else:
            nodes = [node for (node, degree) in g_full.degree() if degree <= max_degree_connect]

        # connect all indentified nodes to each other, all-to-all-style
        pairs = itertools.combinations(nodes, 2)
        for pair in pairs:
            dist = great_circle_vec(
                float(g_full.nodes[pair[0]]['x']),
                float(g_full.nodes[pair[0]]['y']),
                float(g_full.nodes[pair[1]]['x']),
                float(g_full.nodes[pair[1]]['y']))
            g_full.add_edge(pair[0], pair[1], required=0, distance=dist)

    return g_full;


def _get_even_or_odd_nodes(graph, mod):
    """
    Helper function for get_even_nodes.  Given a networkx object, return names of the odd or even nodes
    Args:
        graph (networkx graph): determine the degree of nodes in this graph
        mod (int): 0 for even, 1 for odd

    Returns:
        list[str]: list of node names of odd or even degree
    """
    degree_nodes = []
    for v, d in graph.degree():
        if d % 2 == mod:
            degree_nodes.append(v)
    return degree_nodes


def get_odd_nodes(graph):
    """
    Given a networkx object, return names of the odd degree nodes

    Args:
        graph (networkx graph): graph used to list odd degree nodes for

    Returns:
        list[str]: names of nodes with odd degree
    """
    return _get_even_or_odd_nodes(graph, 1)


def get_even_nodes(graph):
    """
    Given a networkx object, return names of the even degree nodes

    Args:
        graph (networkx graph): graph used to list even degree nodes for

    Returns:
        list[str]: names of nodes with even degree

    """
    return _get_even_or_odd_nodes(graph, 0)

def great_circle_vec(lat1, lng1, lat2, lng2, earth_radius=6371009): # meters
    phi1 = np.deg2rad(90 - lat1)
    phi2 = np.deg2rad(90 - lat2)

    theta1 = np.deg2rad(lng1)
    theta2 = np.deg2rad(lng2)

    cos = (np.sin(phi1) * np.sin(phi2) * np.cos(theta1 - theta2) + np.cos(phi1) * np.cos(phi2))
    arc = np.arccos(cos)

    distance = arc * earth_radius
    return distance

def filter_by_haversine_distance(graph, pairs, max_distance=100):
    """
    Filter node pairs by a max distance, discard any pairs too far apart

    Args:
        graph (networkx graph)
        pairs (list[2tuple]): List of length 2 tuples containing node pairs
        max_distance (str): Defines the max distance node pairs are allowed to be apart

    Returns:
        list[2tuple]: The input pairs filtered by a max distance constraint.
    """
    if max_distance is None:
        return pairs
    else:
        filtered_list = []
        for pair in pairs:
            if great_circle_vec(graph.nodes[pair[0]]['x'], graph.nodes[pair[0]]['y'], graph.nodes[pair[1]]['x'], graph.nodes[pair[1]]['y']) < max_distance:
                filtered_list.append(pair)

        print(len(filtered_list), 'down from', len(pairs))
        return filtered_list


def get_shortest_paths_distances(graph, pairs, edge_weight_name='distance'):
    """
    Calculate shortest distance between each pair of nodes in a graph

    Args:
        graph (networkx graph)
        pairs (list[2tuple]): List of length 2 tuples containing node pairs to calculate shortest path between
        edge_weight_name (str): edge attribute used for distance calculation

    Returns:
        dict: mapping each pair in `pairs` to the shortest path using `edge_weight_name` between them.
    """
    distances = {}
    for pair in pairs:
        distances[pair] = nx.dijkstra_path_length(graph, pair[0], pair[1], weight=edge_weight_name)
    return distances


def create_complete_graph(pair_weights, flip_weights=True):
    """
    Create a perfectly connected graph from a list of node pairs and the distances between them.

    Args:
        pair_weights (dict): mapping between node pairs and distance calculated in `get_shortest_paths_distances`.
        flip_weights (Boolean): True negates the distance in `pair_weights`.  We negate whenever we want to find the
         minimum weight matching on a graph because networkx has only `max_weight_matching`, no `min_weight_matching`.

    Returns:
        complete (fully connected graph) networkx graph using the node pairs and distances provided in `pair_weights`
    """
    g = nx.Graph()
    for k, v in pair_weights.items():
        wt_i = -v if flip_weights else v
        g.add_edge(k[0], k[1], **{'distance': v, 'weight': wt_i})
    return g


def dedupe_matching(matching):
    """
    Remove duplicates node pairs from the output of networkx.algorithms.max_weight_matching since we don't care about order.

    Args:
        matching (dict): output from networkx.algorithms.max_weight_matching.  key is "from" node, value is "to" node.

    Returns:
        list[2tuples]: list of node pairs from `matching` deduped (ignoring order).
    """
    matched_pairs_w_dupes = [tuple(sorted([k, v])) for k, v in matching.items()]
    return list(set(matched_pairs_w_dupes))


def add_augmenting_path_to_graph(graph, min_weight_pairs, edge_weight_name='weight'):
    """
    Add the min weight matching edges to the original graph
    Note the resulting graph could (and likely will) have edges that didn't exist on the original graph.  To get the
    true circuit, we must breakdown these augmented edges into the shortest path through the edges that do exist.  This
    is done with `create_eulerian_circuit`.

    Args:
        graph (networkx graph):
        min_weight_pairs (list[2tuples): output of `dedupe_matching` specifying the odd degree nodes to link together
        edge_weight_name (str): edge attribute used for distance calculation

    Returns:
        networkx graph: `graph` augmented with edges between the odd nodes specified in `min_weight_pairs`
    """
    graph_aug = graph.copy()  # so we don't mess with the original graph
    for pair in min_weight_pairs:
        graph_aug.add_edge(pair[0],
                           pair[1],
                           **{'distance': nx.dijkstra_path_length(graph, pair[0], pair[1], weight=edge_weight_name),
                              'augmented': True}
                           )
    return graph_aug


def create_eulerian_circuit(graph_augmented, graph_original, start_node=None, edge_weight='distance'):
    """
    networkx.eulerian_circuit only returns the order in which we hit each node.  It does not return the attributes of the
    edges needed to complete the circuit.  This is necessary for the postman problem where we need to keep track of which
    edges have been covered already when multiple edges exist between two nodes.
    We also need to annotate the edges added to make the eulerian to follow the actual shortest path trails (not
    the direct shortest path pairings between the odd nodes for which there might not be a direct trail)

    Args:
        graph_augmented (networkx graph): graph w links between odd degree nodes created from `add_augmenting_path_to_graph`.
        graph_original (networkx graph): orginal graph created from `create_networkx_graph_from_edgelist`
        start_node (str): name of starting (and ending) node for CPP solution.

    Returns:
        networkx graph (`graph_original`) augmented with edges directly between the odd nodes
    """

    euler_circuit = list(nx.eulerian_circuit(graph_augmented, source=start_node, keys=True))
    assert len(graph_augmented.edges()) == len(euler_circuit), 'graph and euler_circuit do not have equal number of edges.'

    for edge in euler_circuit:
        aug_path = nx.shortest_path(graph_original, edge[0], edge[1], weight='distance')
        edge_attr = graph_augmented[edge[0]][edge[1]][edge[2]]
        if not edge_attr.get('augmented'):
            yield edge + (edge_attr,)
        else:
            for edge_aug in list(zip(aug_path[:-1], aug_path[1:])):
                # find edge with shortest distance (if there are two parallel edges between the same nodes)
                edge_aug_dict = graph_original[edge_aug[0]][edge_aug[1]]
                edge_key = min(edge_aug_dict.keys(), key=(lambda k: edge_aug_dict[k][edge_weight]))  # index with min distance
                edge_aug_shortest = edge_aug_dict[edge_key]
                edge_aug_shortest['augmented'] = True
                edge_aug_shortest['id'] = edge_aug_dict[edge_key]['id']
                yield edge_aug + (edge_key, edge_aug_shortest, )


def create_required_graph(graph):
    """
    Strip a graph down to just the required nodes and edges.  Used for RPP.  Expected edge attribute "required" with
     True/False or 0/1 values.

    Args:
        graph (networkx MultiGraph):

    Returns:
        networkx MultiGraph with optional nodes and edges deleted
    """

    graph_req = nx.MultiGraph(graph.copy())  # preserve original structure

    # remove optional edges
    for e in list(graph_req.edges(data=True, keys=True)):
        if not e[3]['required']:
            graph_req.remove_edge(e[0], e[1], key=e[2])

    # remove any nodes left isolated after optional edges are removed (no required incident edges)
    for n in list(nx.isolates(graph_req)):
        graph_req.remove_node(n)

    return graph_req


def assert_graph_is_connected(graph):
    """
    Ensure that the graph is still a connected graph after the optional edges are removed.

    Args:
        graph (networkx MultiGraph):

    Returns:
        True if graph is connected
    """

    assert nx.algorithms.connected.is_connected(graph), "Sorry, the required graph is not a connected graph after " \
                                                        "the optional edges are removed.  This is a requirement for " \
                                                        "this implementation of the RPP here which generalizes to the " \
                                                        "CPP."
    return True


def is_connected(graph):
    """
    Ensure that the graph is still a connected graph after the optional edges are removed.

    Args:
        graph (networkx MultiGraph):

    Returns:
        True if graph is connected, False o.w.
    """

    return nx.algorithms.connected.is_connected(graph)




def make_connected(graph, graph_full, edge_weight=None):
    """
    Graph must be a subgraph of graph_full. Add edges to graph to create a single connected component.
    Please note: graph_full keys must all be 0. This is bad, I will try to fix it later...
    Args:
        graph (networkx MultiGraph)
        graph_full (networkx MultiGraph)
        edge_weight (str): name of the graph edge attribute to use as length measurement

    Returns:
        graph_conn (networkx MultiGraph)
    """

    # Test assumptions
    assert nx.algorithms.connected.is_connected(graph_full), "Sorry, graph_full is not a connected -- please make sure the full graph is connected." # full graph must be connected, or this won't work
    assert len([x for x in graph.edges() if x not in graph_full.edges()])==0, "Sorry, the graph is not a subgraph of graph_full!" # graph must be a subgraph of graph_full

    # Connect all subgraph components
    subgraphs =  list(enumerate(graph.subgraph(c) for c in nx.algorithms.connected.connected_components(graph)))
    subgraph_pairs = itertools.combinations(subgraphs, 2)

    subgraph_graph = nx.Graph() # each node is a subgraph of graph. This is to track the shortest path that connects all the subgraphs...

    for G1, G2 in subgraph_pairs:
        # prefer to connect odd nodes...
        nodes1 = get_odd_nodes(G1[1])
        if len(nodes1)==0:
            nodes1 = G1[1].nodes()
        nodes2 = get_odd_nodes(G2[1])
        if len(nodes2)==0:
            nodes2 = G2[1].nodes()

        # track shortest path between nodes
        path_segment = {}
        for n1, n2 in itertools.product(nodes1,nodes2):
            shortest_path_length = nx.algorithms.shortest_paths.generic.shortest_path_length(graph_full, n1, n2, edge_weight)
            if shortest_path_length < path_segment.get('length', inf):
                path_segment['length'] = shortest_path_length
                path_segment['path'] = nx.algorithms.shortest_paths.generic.shortest_path(graph_full, n1, n2, edge_weight)

        subgraph_graph.add_edge(G1[0], G2[0], **path_segment)

    shortest_spanning_tree = nx.minimum_spanning_tree(subgraph_graph, 'length')
    for p in shortest_spanning_tree.edges(data=True):
        # each edge corresponds to a path in g_full that connects two disjoint subgraphs of graph
        for n1, n2  in zip(p[2]["path"], p[2]["path"][1:]):
            # iterate over the edges in the path and add them to graph
            graph.add_edge(n1, n2, **graph_full[n1][n2][0]) # this is where graph_full keys HAVE to all be 0



"""


import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from graph import read_edgelist
from graph import create_networkx_graph_from_edgelist
from graph import get_odd_nodes, get_shortest_paths_distances
from graph import  create_required_graph
from graph import make_connected
df = read_edgelist('edgelist_B.csv', keep_optional=True)
A = create_networkx_graph_from_edgelist(df)
B = create_required_graph(A)
#nx.draw(B); plt.show()
#nx.draw(A); plt.show()
make_connected(B,A,'distance')
arc_wt = nx.get_edge_attributes(B,'distance')
arc_wt2 = dict([(key[:2],val) for key,val in arc_wt.items()])
npos = nx.spring_layout(B)
nx.draw_networkx(B,pos=npos, node_size=450)
nx.draw_networkx_edges(B, pos=npos)
nx.draw_networkx_edge_labels(B, pos = npos, edge_labels=arc_wt2)
plt.show()


"""
