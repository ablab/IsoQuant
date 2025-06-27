import numpy as np
import flowpaths as fp
import networkx as nx

import flowpaths.utils.solverwrapper as sw
import flowpaths.stdigraph as stdigraph
import flowpaths.utils as utils
import flowpaths.nodeexpandeddigraph as nedg
import networkx as nx
from copy import deepcopy

import flowpaths.utils as utils

from copy import deepcopy
from src.cell_type_tree import CellTypeTree

class MinErrorCellTypeFlowCorrection(fp.MinErrorFlow):

    def __init__(
            self, 
            G: nx.DiGraph,
            cell_tree: CellTypeTree,
            flow_attr: str,
            cell_flow_attr: str,
            flow_attr_origin: str = "edge",
            weight_type: type = float,
            sparsity_lambda: float = 0,
            few_flow_values_epsilon: float = None,
            elements_to_ignore: list = [],
            error_scaling: dict = {},
            additional_starts: list = [],
            additional_ends: list = [],
            solver_options: dict = {},   
            ):
        
        self.cell_flow_attr = cell_flow_attr
        self.cell_tree = cell_tree
        self.cell_types = cell_tree.get_leaf_types()
        self.cell_groups = cell_tree.get_cell_types()

       # Handling node-weighted graphs
        self.flow_attr_origin = flow_attr_origin
        if self.flow_attr_origin == "node":
            if G.number_of_nodes() == 0:
                utils.logger.error(f"{__name__}: The input graph G has no nodes. Please provide a graph with at least one node.")
                raise ValueError(f"The input graph G has no nodes. Please provide a graph with at least one node.")
            
            self.G_internal = nedg.NodeExpandedDiGraph(G, node_flow_attr=flow_attr, node_cell_flow_attr=cell_flow_attr)
            additional_starts_internal = self.G_internal.get_expanded_additional_starts(additional_starts)
            additional_ends_internal = self.G_internal.get_expanded_additional_ends(additional_ends)

            edges_to_ignore_internal = self.G_internal.edges_to_ignore
            if not all(isinstance(node, str) for node in elements_to_ignore):
                utils.logger.error(f"elements_to_ignore must be a list of nodes, i.e. strings, not {elements_to_ignore}")
                raise ValueError(f"elements_to_ignore must be a list of nodes, i.e. strings, not {elements_to_ignore}")
            edges_to_ignore_internal += [self.G_internal.get_expanded_edge(node) for node in elements_to_ignore]
            edges_to_ignore_internal = list(set(edges_to_ignore_internal))

            error_scaling_internal = {self.G_internal.get_expanded_edge(node): error_scaling[node] for node in error_scaling}

        elif self.flow_attr_origin == "edge":
            if G.number_of_edges() == 0:
                utils.logger.error(f"{__name__}: The input graph G has no edges. Please provide a graph with at least one edge.")
                raise ValueError(f"The input graph G has no edges. Please provide a graph with at least one edge.")
            
            self.G_internal = G
            additional_starts_internal = additional_starts
            additional_ends_internal = additional_ends

            if not all(isinstance(edge, tuple) and len(edge) == 2 for edge in elements_to_ignore):
                utils.logger.error(f"elements_to_ignore must be a list of edges (i.e. tuples of nodes), not {elements_to_ignore}")
                raise ValueError(f"elements_to_ignore must be a list of edges (i.e. tuples of nodes), not {elements_to_ignore}")
            edges_to_ignore_internal = elements_to_ignore

            error_scaling_internal = error_scaling
        else:
            utils.logger.error(f"flow_attr_origin must be either 'node' or 'edge', not {self.flow_attr_origin}")
            raise ValueError(f"flow_attr_origin must be either 'node' or 'edge', not {self.flow_attr_origin}")

        self.original_graph_copy = deepcopy(self.G_internal)
        self.sparsity_lambda = sparsity_lambda
        
        if nx.is_directed_acyclic_graph(self.G_internal):
            self.is_acyclic = True
            self.G = stdigraph.stDiGraph(self.G_internal, additional_starts=additional_starts_internal, additional_ends=additional_ends_internal)
            self.edges_to_ignore = set(edges_to_ignore_internal).union(self.G.source_sink_edges)
        else:
            self.G = self.G_internal
            self.is_acyclic = False
            self.edges_to_ignore = set(edges_to_ignore_internal)
            if self.sparsity_lambda != 0:
                utils.logger.error(f"{__name__}: You cannot set sparsity_lambda != 0 for a graph with cycles.")
                raise ValueError(f"You cannot set sparsity_lambda != 0 for a graph with cycles.")
        self.edge_error_scaling = error_scaling_internal
        # If the error scaling factor is 0, we ignore the edge
        self.edges_to_ignore |= {edge for edge, factor in self.edge_error_scaling.items() if factor == 0}
        
        # Checking that every entry in self.error_scaling is between 0 and 1
        for key, value in error_scaling.items():
            if value < 0 or value > 1:
                utils.logger.error(f"{__name__}: Error scaling factor for {key} must be between 0 and 1.")
                raise ValueError(f"Error scaling factor for {key} must be between 0 and 1.")

        self.flow_attr = flow_attr
        if weight_type not in [int, float]:
            utils.logger.error(f"{__name__}: weight_type must be either int or float, not {weight_type}")
            raise ValueError(f"weight_type must be either int or float, not {weight_type}")
        self.weight_type = weight_type
        self.solver_options = solver_options

        # Checking that every entry in self.edge_error_scaling is between 0 and 1
        for key, value in self.edge_error_scaling.items():
            if value < 0 or value > 1:
                utils.logger.error(f"{__name__}: Error scaling factor for {key} must be between 0 and 1.")
                raise ValueError(f"Error scaling factor for {key} must be between 0 and 1.")


        self.different_flow_values_epsilon = few_flow_values_epsilon
        if few_flow_values_epsilon is not None:
            if few_flow_values_epsilon < 0:
                utils.logger.error(f"{__name__}: different_flow_values_epsilon must be greater than or equal to 0, not {few_flow_values_epsilon}")
                raise ValueError(f"different_flow_values_epsilon must be greater than or equal to 0, not {few_flow_values_epsilon}")
            if few_flow_values_epsilon == 0:
                self.different_flow_values_epsilon = None        

        self._solution = None
        self._is_solved = None
        self.solve_statistics = dict()

        self.edge_vars = {}
        self.edge_error_vars = {}
        self.edge_sol = {}

        self.w_max = max(
            [
                self.G[u][v].get(self.flow_attr, 0)
                for (u, v) in self.G.edges() 
            ]
        )
        self.ub = self.w_max * self.G.number_of_edges()

        self._create_solver()

        self._encode_flow()

        self._encode_min_sum_errors_objective()  


    def _encode_flow(self):

        # Creating the edge variables
        self.edge_indexes = [(u, v, ct) for (u, v) in self.G.edges() for ct in self.cell_types]
        
        self.edge_vars = self.solver.add_variables(
            self.edge_indexes, 
            name_prefix="edge_vars", 
            lb=0, 
            ub=self.ub, 
            var_type="integer" if self.weight_type == int else "continuous",
        )
        
        self.edge_error_vars = self.solver.add_variables(
            self.edge_indexes, 
            name_prefix="edge_error_vars", 
            lb = -self.ub, 
            ub = self.ub, 
            var_type="integer" if self.weight_type == int else "continuous",
        )

        self.edge_absolute_error_vars = self.solver.add_variables(
            self.edge_indexes, 
            name_prefix="absolute_error_vars", 
            lb = 0, 
            ub = self.ub, 
            var_type="integer" if self.weight_type == int else "continuous",
        )

        # Adding flow conservation constraints
        for node in self.G.nodes():

            if self.G.in_degree(node) == 0 or self.G.out_degree(node) == 0:
                continue

            # Flow conservation constraint
            for ct in self.cell_types:
                self.solver.add_constraint(
                    self.solver.quicksum(
                        self.edge_vars[(u, v, ct)]
                        for (u, v) in self.G.in_edges(node)
                    )
                    - self.solver.quicksum(
                        self.edge_vars[(u, v, ct)]
                        for (u, v) in self.G.out_edges(node)
                    )
                    == 0,
                    name=f"flow_conservation_{node}_{ct}",
                )
        
        for u, v, data in self.G.edges(data=True):
            
            if (u, v) in self.edges_to_ignore:
                for ct in self.cell_types:
                    self.solver.add_constraint(
                        self.edge_error_vars[(u, v, ct)] == 0,
                        name=f"edge_error_u={u}_v={v}",
                    )
                continue
            
            if self.flow_attr not in data:
                utils.logger.error(f"{__name__}: Flow attribute '{self.flow_attr}' not found in edge data for edge {str((u, v))}, and this edge is not in the edges_to_ignore list.")
                raise ValueError(f"Flow attribute '{self.flow_attr}' not found in edge data for edge {str((u, v))}, and this edge is not in the edges_to_ignore list.")

            # Add absolute error calculation
            for ct in self.cell_types:

                self.solver.add_constraint(
                    self.edge_error_vars[(u, v, ct)] <= self.edge_absolute_error_vars[(u, v, ct)]
                )

                self.solver.add_constraint(
                    -self.edge_error_vars[(u, v, ct)] <= self.edge_absolute_error_vars[(u, v, ct)]
                )

            f_u_v = data[self.cell_flow_attr]

            for cg in self.cell_groups:

                if cg == "ROOT":
                    continue
                
                f_u_v_ct = int(f_u_v[self.cell_tree.get_cell_type_index(cg)])
                #print(cg, self.cell_tree.get_child_leafs(cg), self.cell_tree.children[cg])

                self.solver.add_constraint(
                    self.solver.quicksum(
                        self.edge_error_vars[(u, v, ct)] + self.edge_vars[(u, v, ct)]
                        for ct in self.cell_tree.get_child_leafs(cg)
                    ) >= f_u_v_ct,
                    name=f"edge_error_u={u}_v={v}_cg={cg}",
                )

            f_u_v_ct = int(f_u_v[self.cell_tree.get_cell_type_index("ROOT")])

            # Root equality constraint
            self.solver.add_constraint(
                self.solver.quicksum(
                    self.edge_error_vars[(u, v, ct)] + self.edge_vars[(u, v, ct)]
                    for ct in self.cell_types
                ) == f_u_v_ct, #f_u_v_ct,
                name=f"edge_error_u={u}_v={v}",
            )



    def _encode_min_sum_errors_objective(self):
        
        # Objective function: minimize the sum of the edge error variables
        # plus the sparsity of the solution (i.e. sparsity_lambda * sum of the corrected flow going out of the source)
        self.solver.set_objective(
            self.solver.quicksum(
                self.edge_absolute_error_vars[(u, v, ct)] * self.edge_error_scaling.get((u, v), 1)
                for (u, v) in self.G.edges() if (u, v) not in self.edges_to_ignore
                for ct in self.cell_types
            ) + (self.sparsity_lambda * self.solver.quicksum(
                self.edge_vars[(u, v, ct)]
                for ct in self.cell_types
                for (u, v) in self.G.out_edges(self.G.source)) if self.sparsity_lambda > 0 else 0
            ),
            sense="minimize",
        )


    def get_solution(self):
        """
        Returns the solution to the problem, if the model was solved, as a dictionary containing the following keys:

        - `graph`: the corrected graph, as a networkx DiGraph.
        - `error`: the error of the solution, i.e. the sum of the absolute differences between the original weights and the corrected weights.
        - `objective_value`: the value of the objective function.
        
        !!! warning "Warning"
            Call the `solve` method first.
        """
        if self._solution is not None:
            return self._solution
        
        self._check_is_solved()

        edge_sol_dict = self.solver.get_variable_values("edge_vars", [str, str, str])

        for edge in edge_sol_dict.keys():
            self.edge_sol[edge] = (
                round(edge_sol_dict[edge])
                if self.weight_type == int
                else float(edge_sol_dict[edge])
            )
        
        edge_error_sol_dict = self.solver.get_variable_values("edge_error_vars", [str, str, str])
        error = sum(edge_error_sol_dict.values())


        # The vectorization of results
        edge_sol_dict = self.cell_tree.transform_counts_to_arrays(edge_sol_dict)
        corrected_graph = deepcopy(self.original_graph_copy)
        
        for u, v in corrected_graph.edges():
            #print(u, v)
            if self.flow_attr in corrected_graph[u][v]:
                corrected_graph[u][v][self.cell_flow_attr] = edge_sol_dict[(u, v)]
                corrected_graph[u][v][self.flow_attr] = corrected_graph[u][v][self.cell_flow_attr][0]

        if self.flow_attr_origin == "edge":
            self._solution = {
                "graph": corrected_graph,
                "error": error,
                "objective_value": self.solver.get_objective_value(),
            }
        elif self.flow_attr_origin == "node":
            self._solution = {
                "graph": corrected_graph.get_condensed_graph(),
                "error": error,
                "objective_value": self.solver.get_objective_value(),
            }
        
        return self._solution  