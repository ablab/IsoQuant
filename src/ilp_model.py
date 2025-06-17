import flowpaths as fp
import networkx as nx
from collections import defaultdict
from src.intron_graph import VERTEX_polya, VERTEX_polyt, VERTEX_read_end, VERTEX_read_start
import logging
from itertools import chain
import pickle
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

def clean_graph(G,skip_isolated_nodes):
    if skip_isolated_nodes:
        nodes_to_remove = []
        for node in G.nodes():
            if G.in_degree(node) == 0 and G.out_degree(node) == 0:
                nodes_to_remove.append(node)
        G.remove_nodes_from(nodes_to_remove)

def Intron2Nx_Node(
        intron_graph,
        skip_isolated_nodes=True,
        skip_terminal_nodes=True,
):
    G = nx.DiGraph()
    # We create the full networkx graph
        # We add all the nodes
    for intron, this_flow in intron_graph.intron_collector.clustered_introns.items(): #this only adds the internal vertices, we still need to add the start vertices and end vertices
         G.add_node(str(intron),flow=this_flow)

    for intron in chain.from_iterable(chain(intron_graph.outgoing_edges.values(), intron_graph.incoming_edges.values())):
        if intron not in G:
            G.add_node(str(intron))


        # We add all the edges
    additional_starts = []
    additional_ends = []
    edges_to_ignore = []

    for intron in intron_graph.incoming_edges.keys():
        for preceding_intron in intron_graph.incoming_edges[intron]:
            print("Prec_introns",preceding_intron)
            if skip_terminal_nodes:
                if preceding_intron[0] in [VERTEX_polyt, VERTEX_read_start]:
                    additional_starts.append(str(intron))
                else:
                    G.add_edge(str(preceding_intron), str(intron))
            else:
                G.add_edge(str(preceding_intron), str(intron))
                if preceding_intron[0] in [VERTEX_polyt, VERTEX_read_start]:
                    additional_starts.append(str(preceding_intron))
                    edges_to_ignore.append((str(preceding_intron), str(intron)))
    for intron in intron_graph.outgoing_edges.keys():
        for subsequent_intron in intron_graph.outgoing_edges[intron]:
            print("Subs_introns", subsequent_intron)
            if skip_terminal_nodes:
                if subsequent_intron[0] in [VERTEX_polya, VERTEX_read_end]:
                    additional_ends.append(str(intron))
                else:
                    G.add_edge(str(intron), str(subsequent_intron))
            else:
                G.add_edge(str(intron), str(subsequent_intron))
                if subsequent_intron[0] in [VERTEX_polya, VERTEX_read_end]:
                    additional_ends.append(str(subsequent_intron))
                    edges_to_ignore.append((str(intron), str(subsequent_intron)))
    clean_graph(G, skip_isolated_nodes)


    logger.debug(G.edges(data=True))

    return G, additional_starts, additional_ends, edges_to_ignore


def transfer_constraints(transcripts_constrints,graph):
    pc_transformed = []
    for constraint in transcripts_constrints:
        this_constraint=[]
        #all_edges_in_graph = True
        for first, second in zip(constraint, constraint[1:]):
            this_constraint.append((str(first),str(second)))
            #for edge in this_constraint:
            #    if not (edge in graph.edges()):
            #        all_edges_in_graph = False
        if len(this_constraint) > 0:  # and all_edges_in_graph:
            pc_transformed.append(this_constraint)
    print(pc_transformed)
    return pc_transformed

import time

def export_data(graph, additional_starts,additional_ends,edges_to_ignore,constraints):
    nx.write_gml(graph, "graph.gml")
    add_startsfile = open('add_starts', 'wb')
    pickle.dump(additional_starts, add_startsfile)
    add_startsfile.close()
    add_endsfile = open('add_ends', 'wb')
    pickle.dump(additional_ends, add_endsfile)
    add_endsfile.close()
    edges_ignore = open('edges_ignored', 'wb')
    pickle.dump(edges_to_ignore, edges_ignore)
    edges_ignore.close()
    add_constraintsfile = open('constraints','wb')
    pickle.dump(constraints,add_constraintsfile)
    add_constraintsfile.close()

def remove_nodes(nodes_to_remove, datastructure):
    for node in nodes_to_remove:
        datastructure.remove(node)
        
def filter_constraints(graph, additional_starts,additional_ends):
    startnodes_missing=[]
    endnodes_missing=[]
    for startnode in additional_starts:
        if not( startnode in graph.nodes()):
            startnodes_missing.append(startnode)
    for endnode in additional_ends:
        if not(endnode in graph.nodes()):
            endnodes_missing.append(endnode)
    remove_nodes(startnodes_missing,additional_starts)
    remove_nodes(endnodes_missing,additional_ends)



def ILP_Solver_Nodes(intron_graph,chr_id, gene_id,transcripts_constraints=[],ground_truth_isoforms=[], epsilon=0.25, timeout=300, threads=5):
    print("constraints", transcripts_constraints)
    print("Running ILP part")
    export = False

    graph, additional_starts, additional_ends, edges_to_ignore = Intron2Nx_Node(intron_graph)
    constraints = transfer_constraints(transcripts_constraints, graph)
    if len(constraints) == 0:
        constraints = None
    #filter_constraints(graph, additional_starts,additional_ends)
    for key in intron_graph.edge_weights.keys():
        print(key,", ",intron_graph.edge_weights[key])
    # graph = Intron2Nx_old(intron_graph)
    print(chr_id, " ", gene_id)
    print("constraints", constraints)
    print("Nodes",graph.nodes())
    if not(len(graph.nodes()) == 0 or len(graph.edges())== 0):
        #print("Edges with data")
        #print(graph.edges(data=True))
        if export:
            export_data(graph, additional_starts,additional_ends,edges_to_ignore,constraints)
        fp.utils.draw(
            G=graph,
            flow_attr="flow",
            filename=chr_id+"_"+gene_id+"_"+str(id(graph)) + "graph.png",  # this will be used as filename
            draw_options={
                "show_graph_edges": True,
                "show_edge_weights": True,
                "show_path_weights": False,
                "show_node_weights": True,
                "show_path_weight_on_first_edge": True,
                "pathwidth": 2,

            },
            additional_starts = additional_starts,
            additional_ends = additional_ends,
            subpath_constraints = constraints,)
        print(graph.nodes())
        print("Running MinErrorFlow")
        correction_model = fp.MinErrorFlow(
            G=graph,
            flow_attr="flow",
            flow_attr_origin="node",
            weight_type=int,
            additional_starts=additional_starts,
            additional_ends=additional_ends,
        )
        correction_model.solve()
        corrected_graph = correction_model.get_corrected_graph()
        optimization_options = {
            "optimize_with_safe_paths": False,
            "optimize_with_safe_sequences": True,
            "optimize_with_safety_from_largest_antichain": True,
        }
        #solver_options = {
        #    "threads": 1,
        #    "time_limit": 600,
        #}




        print("Running MinFlowDecomp")
        mfd_model = fp.MinFlowDecomp(
            G=corrected_graph,
            flow_attr="flow",
            flow_attr_origin="node",
            additional_starts=additional_starts,
            additional_ends=additional_ends,
            subpath_constraints=constraints,
            optimization_options=optimization_options,
            #solver_options=solver_options,
        )

        # draw the ground truth isoforms , might yield bugs (if partaking nodes are not part of the current graph)!!
        gtweights = [1] * len(ground_truth_isoforms)
        fp.utils.draw(
            G=graph,
            flow_attr="flow",
            paths=ground_truth_isoforms,
            weights=gtweights,
            filename=chr_id+"_"+gene_id+"_"+str(id(graph)) + "groundtruth.png",  # this will be used as filename
            draw_options={
                "show_graph_edges": True,
                "show_edge_weights": False,
                "show_path_weights": False,
                "show_path_weight_on_first_edge": True,
                "pathwidth": 2,
            },
            additional_starts=additional_starts,
            additional_ends=additional_ends, )
        start = time.time()
        mfd_model.solve()
        end = time.time()
        print("Solution found! in ", end - start, " seconds")
        solution = process_solution(graph, mfd_model,additional_starts,additional_ends)
        print("solution",solution)
        # Condensing the paths in the expanded graph to paths in the the original graph
        original_paths = solution["paths"]
        weights = solution["weights"]
        print("original paths",original_paths)
        if len(original_paths) != len(weights):
            raise ValueError("The number of paths and weights must be the same.")

        res = list(zip(original_paths, weights))
    else:
        res=[]
    # this_k = 5
    # mpe_model = fp.kMinPathError(graph, flow_attr="flow", k=this_k, weight_type=float)
    # mpe_model.solve()
    # process_solution(
    #    graph=graph,
    #    model=min_path_error_model,
    #    additional_starts=additional_starts,
    #    additional_ends=additional_ends)
    return res


def process_solution(
        graph: nx.DiGraph, 
        model: fp.kMinPathError,

        additional_starts: list = [],
        additional_ends: list = []):
    if model.is_solved():
        solution = model.get_solution()
        print(solution)
        print(model.solve_statistics)
        print("model.is_valid_solution()", model.is_valid_solution())
        fp.utils.draw(
            G=graph,
            flow_attr="flow",
            paths = solution["paths"],
            weights = solution["weights"],
            filename = str(id(graph))+"solution.png", # this will be used as filename
            draw_options = {
            "show_graph_edges": True,
            "show_edge_weights": True,
            "show_path_weights": False,
            "show_path_weight_on_first_edge": True,
            "pathwidth": 2,
        },
        additional_starts = additional_starts,
        additional_ends = additional_ends,)

        return solution
    else:
        print("Model could not be solved.")