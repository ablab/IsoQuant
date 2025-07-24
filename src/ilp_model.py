import time

import flowpaths as fp
import networkx as nx
from collections import defaultdict
from src.intron_graph import VERTEX_polya, VERTEX_polyt, VERTEX_read_end, VERTEX_read_start
import logging
from itertools import chain
import pickle
logger = logging.getLogger('IsoQuant')

from .cell_type_tree import CellTypeTree
from .models.minflowcelldcomp import MinFlowCellDecomp
from .models.minerror_celltypeflow import MinErrorCellTypeFlow
from .models.cellflowcorrection import MinErrorCellTypeFlowCorrection

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

if mfd_solver.is_solved(): # We get the solutionFalse
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


def Constraints_Transfer_Format(input_constraints, skip_isolated_nodes = True, skip_terminal_nodes = True,):
    transferred_constraints = []
    for pc_entry in input_constraints:
        #print("pcent",pc_entry)
        constraint_list = []
        for node in pc_entry:
            if skip_terminal_nodes:
                if not node[0] in [VERTEX_polyt, VERTEX_read_start] and not node[0] in [VERTEX_polya, VERTEX_read_end]:
                    constraint_list.append(node)
            else:
                constraint_list.append(node)
        #print("c_list",constraint_list)
        if len(constraint_list) > 1:
            this_constraint = []
            for first, second in zip(constraint_list, constraint_list[1:]):
                this_constraint.append((str(first), str(second)))
            transferred_constraints.append(this_constraint)
    #for constraint in transferred_constraints:
        #if len(constraint)==0:

    return transferred_constraints


def Intron2Nx_Node(intron_graph, skip_isolated_nodes = True, skip_terminal_nodes = True,):

    G = nx.DiGraph()
    # We create the full networkx graph
        # We add all the nodes
    cell_types = set()
    # First construct the tree to allow right cell type array assignement
    for intron, this_flow in intron_graph.intron_collector.clustered_introns.items(): #this only adds the internal vertices, we still need to add the start vertices and end vertices
        introns_cell_types = [cell_type for cell_type in intron_graph.intron_collector.clustered_introns_by_cell_type[intron]]
        cell_types.update(introns_cell_types)

    #print(cell_types)
    cell_type_tree = CellTypeTree(cell_types)

    for intron, this_flow in intron_graph.intron_collector.clustered_introns.items(): #this only adds the internal vertices, we still need to add the start vertices and end vertices
        #print(intron_graph.intron_collector.clustered_introns_by_cell_type[intron])
        G.add_node(str(intron),
                    flow = this_flow,
                    cell_flow = cell_type_tree.transform_counts(intron_graph.intron_collector.clustered_introns_by_cell_type[intron])
                    )
        

    for intron in chain.from_iterable(chain(intron_graph.outgoing_edges.values(), intron_graph.incoming_edges.values())):
        if intron not in G:
            G.add_node(str(intron))

    # We add all the edges
    additional_starts = []
    additional_ends = []
    edges_to_ignore = []

    for intron in intron_graph.incoming_edges.keys():
        for preceding_intron in intron_graph.incoming_edges[intron]:
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

    return G, cell_type_tree, additional_starts, additional_ends, edges_to_ignore



def transfer_constraints(transcripts_constrints, graph):
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
    #print(pc_transformed)
    return pc_transformed


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
    remove_nodes(startnodes_missing, additional_starts)
    remove_nodes(endnodes_missing, additional_ends)


def ILP_Solver_Nodes(intron_graph, chr_id, gene_id, index, transcripts_constraints: list = [], ground_truth_isoforms: list = [], epsilon: float = 0.25, timeout: float = 300, threads: int = 5, draw_graphs: bool = True):
    #print("constraints", transcripts_constraints)
    #print("Running ILP part")
    export = False

    graph, cell_type_tree, additional_starts, additional_ends, edges_to_ignore = Intron2Nx_Node(intron_graph)
    constraints = Constraints_Transfer_Format(transcripts_constraints)
    #print(constraints)
    #constraints = transfer_constraints(transcripts_constraints, graph)
    
    if not(len(graph.nodes()) == 0 or len(graph.edges())== 0):

        if export:
            export_data(graph, additional_starts,additional_ends,edges_to_ignore,constraints)

        if draw_graphs:
            # Draw the graph with normal flows    
            fp.utils.draw(
                G = graph,
                flow_attr = "flow",
                filename = "graphs/" + chr_id + "_" + gene_id + "_" + str(index) + "_" + str(id(graph))+ "graph.FL.png",
                draw_options = {
                    "show_graph_edges": True,
                    "show_edge_weights": True,
                    "show_path_weights": False,
                    "show_node_weights": True,
                    "show_path_weight_on_first_edge": True,
                    "pathwidth": 2,
                },
                additional_starts = additional_starts,
                additional_ends = additional_ends,
                subpath_constraints = [], #constraints,
                )

            # Draw the graph with cell type flows
            fp.utils.draw(
                G = graph,
                flow_attr = "cell_flow",
                filename = "graphs/" + chr_id + "_" + gene_id + "_" + str(index) + "_" + str(id(graph)) + "graph.CT.png",  # this will be used as filename
                draw_options = {
                    "show_graph_edges": True,
                    "show_edge_weights": True,
                    "show_path_weights": False,
                    "show_node_weights": True,
                    "show_path_weight_on_first_edge": True,
                    "pathwidth": 2,
                },
                additional_starts = additional_starts,
                additional_ends = additional_ends,
                subpath_constraints = [], #constraints,
                )
            
        #print(graph.nodes())
        logger.info("Running MinErrorCellTypeFlow")
        
        additional_starts_pruned = []
        additional_ends_pruned = []
        subpath_constaints_pruned = []
        for node in additional_starts:
            if node in graph.nodes:
                additional_starts_pruned.append(node)
        for node in additional_ends:
            if node in graph.nodes:
                additional_ends_pruned.append(node)
        for path in constraints:
            include = True
            for edge in path:
                if edge not in graph.edges:
                    include = False; break
            if include: subpath_constaints_pruned.append(path)

        correction_model = MinErrorCellTypeFlowCorrection(
            G = graph,
            cell_tree = cell_type_tree,
            flow_attr = "flow",
            cell_flow_attr = "cell_flow",
            flow_attr_origin = "node",
            weight_type = int,
            additional_starts = additional_starts,
            additional_ends = additional_ends,
        )

        correction_model.solve()
        corrected_graph = correction_model.get_corrected_graph()
        
         # Draw the graph with cell type flows
        if draw_graphs:
            fp.utils.draw(
                G =corrected_graph,
                flow_attr = "cell_flow",
                filename = "graphs/" + chr_id + "_" + gene_id + "_" +  str(index) + "_" + str(id(graph))+ "graph.corrected.CT.png",  # this will be used as filename
                draw_options = {
                    "show_graph_edges": True,
                    "show_edge_weights": True,
                    "show_path_weights": False,
                    "show_node_weights": True,
                    "show_path_weight_on_first_edge": True,
                    "pathwidth": 2,
                },
                additional_starts = additional_starts,
                additional_ends = additional_ends,
                #subpath_constraints = constraints,
                )
            
            fp.utils.draw(
                G =corrected_graph,
                flow_attr = "flow",
                filename = "graphs/" + chr_id + "_" + gene_id + "_" +  str(index) + "_" + str(id(graph)) + "graph.corrected.FL.png",  # this will be used as filename
                draw_options = {
                    "show_graph_edges": True,
                    "show_edge_weights": True,
                    "show_path_weights": False,
                    "show_node_weights": True,
                    "show_path_weight_on_first_edge": True,
                    "pathwidth": 2,
                },
                additional_starts = additional_starts,
                additional_ends = additional_ends,
                #subpath_constraints = constraints,
                )
            

        optimization_options = {
            "optimize_with_safe_paths": True,
            "optimize_with_safe_sequences": False,
            "optimize_with_safety_from_largest_antichain": True,
        }
    
        logger.info("Running MinFlowDecomp")
        
    
        mcd_model = MinFlowCellDecomp(
            G = corrected_graph,
            flow_attr = "flow",
            cell_flow_attr = "cell_flow",
            cell_tree = cell_type_tree,
            flow_attr_origin = "node",
            additional_starts = additional_starts,
            additional_ends = additional_ends,
            subpath_constraints = subpath_constaints_pruned,
            optimization_options = optimization_options,
        )

        start = time.time()
        mcd_model.solve()
        end = time.time()
        logger.info("Solution found! in %d seconds", end - start)
        solution = process_solution(graph, mcd_model, additional_starts,additional_ends)
        #print("solution",solution)
        # Condensing the paths in the expanded graph to paths in the the original graph
        
        if solution is None: return []
        
        original_paths = solution["paths"]
        weights = solution["weights"]
        #print("original paths",original_paths)

        solution = mcd_model.get_solution()

        if solution is None: return []

        paths = solution["paths"]
        weights = solution["weights"]
        ct_weights = [cell_type_tree.transfrom_count_array_to_dict(sol) for sol in solution["ct_weights"]]

        if draw_graphs:
            fp.utils.draw(
                G = graph,
                flow_attr = "flow",
                paths = paths,
                weights = weights,
                filename = "graphs/" + chr_id + "_" + gene_id + "_" +  str(index) + "_" + str(id(graph))  + "graph.solution.CT.png",
                draw_options= {
                        "show_graph_edges": True,
                        "show_edge_weights": False,
                        "show_path_weights": False,
                        "show_path_weight_on_first_edge": True,
                        "pathwidth": 2,
                    },
            )


        if len(original_paths) != len(weights):
            raise ValueError("The number of paths and weights must be the same.")

        return list(zip(original_paths, weights, ct_weights))
    
    return []


def process_solution(graph: nx.DiGraph, model: MinFlowCellDecomp, additional_starts: list = [], additional_ends: list = []):
    
    if model.is_solved():
    
        solution = model.get_solution()
    
        #fp.utils.draw(
        #    G=graph,
        #    flow_attr="cell_type",
        #    paths = solution["paths"],
        #    weights = solution["ct_weights"],
        #    filename = str(id(graph))+"solution.png", # this will be used as filename
        #   draw_options = {
        #   "show_graph_edges": True,
        #    "show_edge_weights": True,
        #    "show_path_weights": False,
        #    "show_path_weight_on_first_edge": True,
        #    "pathwidth": 2,
        #},
        #additional_starts = additional_starts,
        #additional_ends = additional_ends,)

        return solution
    else:
        print("Model could not be solved.")