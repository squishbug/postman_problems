import itertools
import logging
import networkx as nx
import time
import warnings
import pdb

from postman_problems.graph import read_edgelist, create_networkx_graph_from_edgelist, read_graphml, create_required_graph, \
    assert_graph_is_connected, get_odd_nodes, get_shortest_paths_distances, create_complete_graph, dedupe_matching, \
    add_augmenting_path_to_graph, create_eulerian_circuit, filter_by_haversine_distance, great_circle_vec, is_connected, make_connected

logger_rpp = logging.getLogger('{0}.{1}'.format(__name__, 'rpp'))
logger_cpp = logging.getLogger('{0}.{1}'.format(__name__, 'cpp'))


def rpp(edgelist_filename=None, start_node=None, edge_weight='distance', verbose=False, graphml=False, max_distance=None, max_degree_connect=None, g_full=None):
    """
    Solving the RPP from beginning (load network data) to end (finding optimal route).  This optimization makes a
     relatively strong assumption: the starting graph must stay a connected graph when optional edges are removed.
    If this is not so, an assertion is raised.  This class of RPP generalizes to the CPP strategy.

    Args:
        edgelist_filename (str): filename of edgelist.  See cpp.py for more details
        start_node (str or can be cast to str): name of starting node.  See cpp.py for more details
        edge_weight (str): name edge attribute that indicates distance to minimize in CPP
        verbose (boolean): log info messages?
        graphml (boolean): is edgelist filename a in graphml format?
        max_distance (double): NOT IMPLEMENTED
        max_degree_connect (int): min degree of a node in the full graph -- nodes with smaller degree are connected with all-to-all optional edges. Use -1 for all-to-all graph.
        g_full (networkx multigraph): pre-loaded networkx MultiGraph. Either g_full or edgelist_filename must be specified. If both are given, filename will be used.

    Returns:
        tuple(list[tuple(str, str, dict)], networkx.MultiGraph]:
        Each tuple is a direction (from one node to another) from the CPP solution route.
          The first element is the starting ("from") node.
          The second element is the end ("to") node.
          The third element is the dict of edge attributes for that edge.
        The original graph is returned as well.  This is needed for visualization
    """

    print("Running RPP solver!")

    #pdb.set_trace()
    
    logger_rpp.disabled = not verbose
    logger_rpp.info('initialize full graph')

    reset_ids = False

    if edgelist_filename is not None:
        # edgelist filename is given - load graph from file

        if graphml:
            # read in the graph
            g_full = read_graphml(edgelist_filename, edge_weight, max_degree_connect)

            # make sure edge id exists and is unique
            shared_keys = set.intersection(*[set(z.keys()) for x,y,z in list(g_full.edges(data=True))])
            if 'id' not in shared_keys:
                reset_ids = True
            else:
                # id is already specified - ensure that it is unique
                if len({edg[3]['id'] for edg in g_full.edges(keys=True, data=True)}) != g_full.number_of_edges():
                    warnings.warn("Edgelist contains field named 'id' but the values provided are not unique."
                                  "Replacing id field with uniquely defined values.")
                    #raise ValueError("If id is specified on edges of g_full it must be unique!")
                    reset_ids = True

        else:
            # regular csv file format...
            el = read_edgelist(edgelist_filename, keep_optional=True)
            g_full = create_networkx_graph_from_edgelist(el)
    elif g_full is None:
        # none of edgelist filename or g_full is given - no graph specified
        raise TypeError("One of edgelist_filename or g_full must be given!")
    else:
        # use g_full - must ensure that format matches the expected format
        g_full = nx.MultiGraph(g_full)
        # check for all needed fields - if id is not set it will be set manually
        shared_keys = set.intersection(*[set(z.keys()) for x,y,z in list(g_full.edges(data=True))])
        if not all([x in shared_keys for x in {'required',edge_weight}]):
            raise ValueError("g_full must include values for 'required' and '{}' for every edge".format(edge_weight))
        if 'id' not in shared_keys:
            # not every edge has a defined edge id - create a new one.
            reset_ids = True
        else:
            # id is already specified - ensure that it is unique
            if len({edg[3]['id'] for edg in g_full.edges(keys=True, data=True)}) != g_full.number_of_edges():
                warnings.warn("Edgelist contains field named 'id' but the values provided are not unique."
                              "Replacing id field with uniquely defined values.")
                reset_ids = True

    # if needed, create new id
    if reset_ids:
        for ii, edg in enumerate(g_full.edges(keys=True)):
            g_full.edges[edg]['id'] = str(ii)

    # if start node is given, make sure it's a string!
    if start_node is not None:
        start_node = str(start_node)

    # if required graph is not connected, use additional edges from g_full to make it connected
    logger_rpp.info('create required graph')
    g_req = create_required_graph(g_full)
    if not is_connected(g_req):
        make_connected(g_req, g_full, edge_weight) # THIS STEP COULD BE SLOW

    logger_rpp.info('getting odd node pairs')
    odd_nodes = get_odd_nodes(g_req)
    odd_node_pairs = list(itertools.combinations(odd_nodes, 2))

    start = time.time()
    logger_rpp.info('get shortest paths between odd nodes')
    odd_node_pairs_shortest_paths = get_shortest_paths_distances(g_full, odd_node_pairs, edge_weight)

    logger_rpp.info('Find min weight matching using blossom algorithm')
    g_odd_complete = create_complete_graph(odd_node_pairs_shortest_paths, flip_weights=True)
    odd_matching = dedupe_matching(nx.algorithms.max_weight_matching(g_odd_complete, True))

    logger_rpp.info('add the min weight matching edges to g')
    g_aug = add_augmenting_path_to_graph(g_req, odd_matching)

    logger_rpp.info('get eulerian circuit route')

    #pdb.set_trace();

    circuit = list(create_eulerian_circuit(g_aug, g_full, start_node, edge_weight=edge_weight))
    end = time.time()
    print('matching and augment time:', end - start)

    # Remove already visited nodes starting from the back (since we dont care about the "full circuit")
    new_ending_idx = len(circuit) - 1;
    for idx in range(0, len(circuit), 1):
        end_offset_idx = len(circuit) - 1 - idx
        if circuit[idx][0] == circuit[end_offset_idx][0] or circuit[idx][0] == circuit[end_offset_idx][1] or circuit[idx][1] == circuit[end_offset_idx][0] or circuit[idx][1] == circuit[end_offset_idx][1]:
            new_ending_idx = end_offset_idx
        else:
            break

    circuit = circuit[idx+1:]
    print('Removed', idx, 'edges from the circuit start')

    return circuit, g_full


def cpp(edgelist_filename, start_node=None, edge_weight='distance', verbose=False, graphml=False, max_distance=None, max_degree_connect=0, g=None):
    """
    Solving the CPP from beginning (load network data) to end (finding optimal route).
    Can be run from command line with arguments from cpp.py, or from an interactive Python session (ex jupyter notebook)

    Args:
        edgelist_filename (str): filename of edgelist.  See cpp.py for more details
        start_node (str): name of starting node.  See cpp.py for more details
        edge_weight (str): name edge attribute that indicates distance to minimize in CPP
        verbose (boolean): log info messages?
        graphml (boolean): is edgelist filename a in graphml format?
        max_distance (double): NOT IMPLEMENTED
        max_degree_connect (int): NOT IMPLEMENTED
        g (networkx multigraph): pre-loaded networkx MultiGraph. Either g or edgelist_filename must be specified. If both are given, filename will be used.

    Returns:
        tuple(list[tuple(str, str, dict)], networkx.MultiGraph]:
        Each tuple is a direction (from one node to another) from the CPP solution route.
          The first element is the starting ("from") node.
          The second element is the end ("to") node.
          The third element is the dict of edge attributes for that edge.
        The original graph is returned as well.  This is needed for visualization
    """
    logger_cpp.disabled = not verbose

    reset_ids = False;
        
    logger_cpp.info('initialize graph')
    if edgelist_filename is not None:
        # edgelist filename is given - load graph from file
        if graphml:
            g = read_graphml(edgelist_filename, edge_weight=edge_weight, max_degree_connect=max_degree_connect);

            # make sure edge id exists and is unique
            shared_keys = set.intersection(*[set(z.keys()) for x,y,z in list(g.edges(data=True))])
            if 'id' not in shared_keys:
                reset_ids = True
            else:
                # id is already specified - ensure that it is unique
                if len({edg[3]['id'] for edg in g.edges(keys=True, data=True)}) != g.number_of_edges():
                    warnings.warn("Edgelist contains field named 'id' but the values provided are not unique."
                                  "Replacing id field with uniquely defined values.")
                    #raise ValueError("If id is specified on edges of g_full it must be unique!")
                    reset_ids = True

        else:
            el = read_edgelist(edgelist_filename, keep_optional=False)
            g = create_networkx_graph_from_edgelist(el)
    elif g is None:
        # none of edgelist filename or g is given - no graph specified
        raise TypeError("One of edgelist_filename or g must be given!")
    else:
        # use g - must ensure that format matches the expected format
        g = nx.MultiGraph(g)
        # check for all needed fields - if id is not set it will be set manually
        shared_keys = set.intersection(*[set(z.keys()) for x,y,z in list(g.edges(data=True))])
        if edge_weight not in shared_keys:
            raise ValueError("g must include value for '{}' for every edge".format(edge_weight))
        if 'id' not in shared_keys:
            # create new id
            reset_ids = True
        else:
            # id is already specified - ensure that it is unique
            if len({edg[3]['id'] for edg in g.edges(keys=True, data=True)}) != g.number_of_edges():
                warnings.warn("Edgelist contains field named 'id' but the values provided are not unique."
                              "Replacing id field with uniquely defined values.")
                reset_ids = True

    # if needed, create new id
    if reset_ids:
        for ii, edg in enumerate(g.edges(keys=True)):
            g.edges[edg]['id'] = str(ii)

    # if start node is given, make sure it's a string!
    if start_node is not None:
        start_node = str(start_node)

    logger_cpp.info('get augmenting path for odd nodes')
    odd_nodes = get_odd_nodes(g)
    odd_node_pairs = list(itertools.combinations(odd_nodes, 2))

    # 'x' and 'y' is not in the generated graphml file, so this filtering is not supported until x and y is added
    # odd_node_pairs = filter_by_haversine_distance(g, odd_node_pairs, max_distance=max_distance)

    start = time.time()
    odd_node_pairs_shortest_paths = get_shortest_paths_distances(g, odd_node_pairs, edge_weight)
    g_odd_complete = create_complete_graph(odd_node_pairs_shortest_paths, flip_weights=True)

    logger_cpp.info('Find min weight matching using blossom algorithm')
    odd_matching = dedupe_matching(nx.algorithms.max_weight_matching(g_odd_complete, True))

    logger_cpp.info('add the min weight matching edges to g')
    g_aug = add_augmenting_path_to_graph(g, odd_matching)

    print(len(get_odd_nodes(g)), ' odd nodes, now', len(get_odd_nodes(g_aug)), nx.is_connected(g_aug))
    logger_cpp.info('get eulerian circuit route')

    #pdb.set_trace();
    
    circuit = list(create_eulerian_circuit(g_aug, g, start_node))
    end = time.time()
    print('matching and augment time:', end - start)

    # Remove already visited nodes starting from the back (since we dont care about the "full circuit")
    new_ending_idx = len(circuit) - 1
    for idx in range(0, len(circuit), 1):
        end_offset_idx = len(circuit) - 1 - idx
        if circuit[idx][0] == circuit[end_offset_idx][0] or circuit[idx][0] == circuit[end_offset_idx][1] or circuit[idx][1] == circuit[end_offset_idx][0] or circuit[idx][1] == circuit[end_offset_idx][1]:
            new_ending_idx = end_offset_idx
        else:
            break

    circuit = circuit[idx+1:]
    print('Removed', idx, 'edges from the circuit start')

    return circuit, g
