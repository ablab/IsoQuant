import flowpaths as fp
import networkx as nx
from collections import defaultdict
from src.intron_graph import VERTEX_polya, VERTEX_polyt, VERTEX_read_end, VERTEX_read_start
import logging
from itertools import chain

logger = logging.getLogger('IsoQuant')

'''
# Create a simple graph
graph = nx.DiGraph()
graph.add_edge("s", "a", flow=2)
graph.add_edge("a", "t", flow=2)
graph.add_edge("s", "b", flow=5)
graph.add_edge("b", "t", flow=5)
# ...

# Create a Minimum Flow Decomposition solver
mfd_solver = fp.MinFlowDecomp(graph, flow_attr="flow")

mfd_solver.solve() # We solve it

if mfd_solver.is_solved(): # We get the solution
    print(mfd_solver.get_solution())
    # {'paths': [['s', 'b', 't'], ['s', 'a', 't']], 'weights': [5, 2]}
'''
def node2string(nodetup):
    (x,y)=nodetup
    return ""+x+","+y


def Intron2Nx(
        intron_graph, 
        skip_isolated_nodes=True,
        skip_terminal_nodes=True,
    ):

    G = nx.DiGraph()

    # We create the full networkx graph

    if not skip_isolated_nodes:
        # We add all the nodes
        for intron in intron_graph.intron_collector.clustered_introns.keys(): #this only adds the internal vertices, we still need to add the start vertices and end vertices
            G.add_node(str(intron))
        for intron in chain.from_iterable(chain(intron_graph.outgoing_edges.values(), intron_graph.incoming_edges.values())):
            if intron not in G:
                G.add_node(str(intron))

    # We add all the edges
    additional_starts = []
    additional_ends = []
    edges_to_ignore = []

    for intron in intron_graph.incoming_edges.keys():
        for preceding_intron in intron_graph.incoming_edges[intron]:
            edge_weight = intron_graph.edge_weights[(preceding_intron, intron)]
            if skip_terminal_nodes:
                if preceding_intron[0] in [VERTEX_polyt, VERTEX_read_start]:
                    additional_starts.append(str(intron))
                else:
                    G.add_edge(str(preceding_intron), str(intron), flow=edge_weight)
            else:
                G.add_edge(str(preceding_intron), str(intron), flow=edge_weight)
                if preceding_intron[0] in [VERTEX_polyt, VERTEX_read_start]:
                    additional_starts.append(str(preceding_intron))
                    edges_to_ignore.append((str(preceding_intron), str(intron)))

    for intron in intron_graph.outgoing_edges.keys():
        for subsequent_intron in intron_graph.outgoing_edges[intron]:
            edge_weight = intron_graph.edge_weights[(intron, subsequent_intron)]
            if skip_terminal_nodes:
                if subsequent_intron[0] in [VERTEX_polya, VERTEX_read_end]:
                    additional_ends.append(str(intron))
                else:
                    G.add_edge(str(intron), str(subsequent_intron), flow=edge_weight)
            else:
                G.add_edge(str(intron), str(subsequent_intron), flow=edge_weight)
                if subsequent_intron[0] in [VERTEX_polya, VERTEX_read_end]:
                    additional_ends.append(str(subsequent_intron))
                    edges_to_ignore.append((str(intron), str(subsequent_intron)))
    
    logger.debug(G.edges(data=True))

    return G, additional_starts, additional_ends, edges_to_ignore

def Intron2Nx_old(intron_graph):
    intron2vertex = dict()
    vertex2intron = dict()
    vertex_id = 1
    G = nx.DiGraph()

    #TODO add id to nxGraph
    #graph.graph["id"] = "simple_graph"
    # add all intron vertices
    for intron in intron_graph.intron_collector.clustered_introns.keys(): #this only adds the internal vertices, we still need to add the start vertices and end vertices
        intron2vertex[intron] = vertex_id
        vertex2intron[vertex_id] = intron
        vertex_id += 1
    #we add the introns that edges lead to
    for intron_set in intron_graph.outgoing_edges.values():
        for intron in intron_set:
            if intron in intron2vertex: continue
            intron2vertex[intron] = vertex_id
            vertex2intron[vertex_id] = intron
            vertex_id += 1
    #we add the edges that edges come from
    for intron_set in intron_graph.incoming_edges.values():
        for intron in intron_set:
            if intron in intron2vertex: continue
            intron2vertex[intron] = vertex_id
            vertex2intron[vertex_id] = intron
            vertex_id += 1

    source = "s"
    target = "t"

    nodes = set()
    #print("incoming",intron_graph.incoming_edges)
    # create edges
    edge_list = []
    starting_introns = defaultdict(int)
    for intron in intron_graph.incoming_edges.keys():
        for preceding_intron in intron_graph.incoming_edges[intron]:
            if intron_graph.edge_weights[(preceding_intron, intron)] > 0:
                edge_list.append((str(intron2vertex[preceding_intron]), str(intron2vertex[intron])))
                if preceding_intron[0] in [VERTEX_polyt, VERTEX_read_start]:
                    starting_introns[preceding_intron] += intron_graph.edge_weights[(preceding_intron, intron)]

    terminal_introns = defaultdict(int)
    for intron in intron_graph.outgoing_edges.keys():
        for subsequent_intron in intron_graph.outgoing_edges[intron]:
            if subsequent_intron[0] in [VERTEX_polya, VERTEX_read_end]:
                terminal_introns[subsequent_intron] += intron_graph.edge_weights[(intron, subsequent_intron)]

    flow_dict = dict()
    for intron in intron_graph.incoming_edges.keys():
        for preceding_intron in intron_graph.incoming_edges[intron]:
            if intron_graph.edge_weights[(preceding_intron, intron)]>0:
                u = str(intron2vertex[preceding_intron])
                v = str(intron2vertex[intron])
                nodes.add(u)
                nodes.add(v)
                flow_dict[(u, v)] = intron_graph.edge_weights[(preceding_intron, intron)]

    # add connection to super source and total weight
    for starting_intron in starting_introns.keys():
        starting_vertex = intron2vertex[starting_intron]
        edge_list.append((str(source), str(starting_vertex)))
        flow_dict[(str(source), str(starting_vertex))] = starting_introns[starting_intron]

    for terminal_intron in terminal_introns.keys():
        terminal_vertex = intron2vertex[terminal_intron]
        edge_list.append((str(terminal_vertex), str(target)))
        flow_dict[(str(terminal_vertex), str(target))] = terminal_introns[terminal_intron]
    #TODO: there seems to be a bug in how the weights are calculated. We should get weights>0
    assert len(edge_list) == len(flow_dict)

    # add nodes and edges to networkx graph
    for node in nodes:
        #print("Adding ", node, ", ",type(node))
        G.add_node(node)
    for (u,v) in edge_list:
        G.add_edge(u,v,flow=flow_dict[(u,v)])
    #print(G.nodes)
    vertex_id += 1 #vertex_id equals the number of nods in G.nx
    logger.debug(vertex2intron)
    logger.debug(flow_dict)
    return G


def ILP_Solver(intron_graph, transcripts_constraints=[], epsilon=0.25, timeout=300, threads=5):
    #for key in intron_graph.edge_weights.keys():
        #print(key,", ",intron_graph.edge_weights[key])
    # graph = Intron2Nx_old(intron_graph)
    graph, additional_starts, additional_ends, edges_to_ignore = Intron2Nx(intron_graph)
    #print(graph.edges(data=True))
    fp.utils.draw_solution_basic(
        graph=graph,
        flow_attr="flow",
        id="this_graph",  # this will be used as filename
        draw_options={
            "show_graph_edges": True,
            "show_edge_weights": True,
            "show_path_weights": False,
            "show_path_weight_on_first_edge": True,
            "pathwidth": 2,
        })
    #graph.graph["id"] = "graph" + str(id(graph))
    min_path_error_model = fp.NumPathsOptimization(
        model_type = fp.kMinPathError,
        stop_on_first_feasible=True,
        G=graph, 
        flow_attr="flow",
        additional_starts=additional_starts,
        additional_ends=additional_ends,
        edges_to_ignore=edges_to_ignore,
        )
    #print("Attempting to solve mpe_model with k=",str(this_k))
    min_path_error_model.solve()
    
    # this_k = 5
    # mpe_model = fp.kMinPathError(graph, flow_attr="flow", k=this_k, weight_type=float)
    # mpe_model.solve()
    process_solution(graph, min_path_error_model)


def process_solution(graph: nx.DiGraph, model: fp.kMinPathError):
    if model.is_solved():
        solution = model.get_solution()
        print(solution)
        print(model.solve_statistics)
        print("model.is_valid_solution()", model.is_valid_solution())
        fp.utils.draw_solution_basic(
            graph=graph,
            flow_attr="flow",
            paths=solution["paths"],
            weights=solution["weights"],
            id=graph.graph["id"], # this will be used as filename
            draw_options={
            "show_graph_edges": True,
            "show_edge_weights": True,
            "show_path_weights": False,
            "show_path_weight_on_first_edge": True,
            "pathwidth": 2,
        })
    else:
        print("Model could not be solved.")