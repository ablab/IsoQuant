import flowpaths as fp
import networkx as nx
from collections import defaultdict
from src.intron_graph import VERTEX_polya, VERTEX_polyt, VERTEX_read_end, VERTEX_read_start
import logging

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


def Intron2Nx(intron_graph):
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
    graph = Intron2Nx(intron_graph)
    #print(graph.edges(data=True))
    this_k = 5
    mpe_model = fp.kMinPathError(graph, flow_attr="flow", k=this_k, weight_type=float)
    print("Attempting to solve mpe_model with k=",str(this_k))
    mpe_model.solve()
    process_solution(mpe_model)


def process_solution(model: fp.kMinPathError):
    if model.is_solved():
        print(model.get_solution())
        print(model.solve_statistics)
        print("model.is_valid_solution()", model.is_valid_solution())
    else:
        print("Model could not be solved.")