import numpy as np
import flowpaths as fp
import networkx as nx

import flowpaths.utils as utils

from copy import deepcopy
from src.cell_type_tree import CellTypeTree

class ExtendedMinErrorCellTypeFlow(fp.MinErrorFlow):

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

        print(self.cell_types)
        print(self.cell_groups)

        super().__init__(
            G = G,
            flow_attr = flow_attr,
            flow_attr_origin = flow_attr_origin,
            weight_type = weight_type,
            sparsity_lambda = sparsity_lambda,
            few_flow_values_epsilon = few_flow_values_epsilon,
            elements_to_ignore = elements_to_ignore,
            error_scaling = error_scaling,
            additional_starts = additional_starts,
            additional_ends = additional_ends,
            solver_options = solver_options,  
        )


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
            name_prefix = "error_vars", 
            lb = 0, 
            ub = self.ub, 
            var_type = "integer" if self.weight_type == int else "continuous",
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
                    name=f"flow_conservation_{node}",
                )
        
        # Encoding the edge error variables
        for u, v, data in self.G.edges(data=True):
            
            if (u, v) in self.edges_to_ignore:
                # Making sure the error of the edges to ignore gets set to 0
                for ct in self.cell_types:
                    self.solver.add_constraint(
                        self.edge_error_vars[(u, v, ct)] == 0,
                        name=f"edge_error_u={u}_v={v}",
                    )
                continue
            
            # If the edge is not in the edges_to_ignore list, we need to check if it has a flow attribute
            if self.flow_attr not in data:
                utils.logger.error(f"{__name__}: Flow attribute '{self.flow_attr}' not found in edge data for edge {str((u, v))}, and this edge is not in the edges_to_ignore list.")
                raise ValueError(f"Flow attribute '{self.flow_attr}' not found in edge data for edge {str((u, v))}, and this edge is not in the edges_to_ignore list.")
            
            # Getting the flow value of the edge            
            f_u_v = data[self.cell_flow_attr]
            
            for cg in self.cell_groups:

                f_u_v_ct = f_u_v[self.cell_tree.get_cell_type_index(cg)]
                #print(cg, f_u_v_ct)

                # Top level constraint
                if cg == "ROOT":

                    print("here", f_u_v_ct, u, v)

                    self.solver.add_constraint(
                        self.solver.quicksum(self.edge_error_vars[(u, v, ct)] for ct in self.cell_types) >= 
                        self.solver.quicksum(self.edge_vars[(u, v, ct)] for ct in self.cell_types) - f_u_v_ct
                    )

                    self.solver.add_constraint(
                        self.solver.quicksum(self.edge_error_vars[(u, v, ct)] for ct in self.cell_types) >= 
                        f_u_v_ct - self.solver.quicksum(self.edge_vars[(u, v, ct)] for ct in self.cell_types)
                    )

                #self.solver.add_constraint(
                #    self.solver.quicksum(self.edge_vars[(u, v, ct)] for ct in self.cell_tree.get_child_leafs(cg)) >=
                #    f_u_v_ct - self.solver.quicksum(self.edge_error_vars[(u, v, ct)] for ct in self.cell_tree.get_child_leafs(cg))
                #)

    def _encode_min_sum_errors_objective(self):
        
        # Objective function: minimize the sum of the edge error variables
        # plus the sparsity of the solution (i.e. sparsity_lambda * sum of the corrected flow going out of the source)
        self.solver.set_objective(
            self.solver.quicksum(
                self.edge_error_vars[(u, v, ct)] * self.edge_error_scaling.get((u, v), 1)
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

        corrected_graph = deepcopy(self.original_graph_copy)
        
        for u, v in corrected_graph.edges():
            if self.flow_attr in corrected_graph[u][v]:
                corrected_graph[u][v][self.cell_flow_attr] = np.array([self.edge_sol[(u, v, ct)] for ct in self.cell_types])
                corrected_graph[u][v][self.flow_attr] = corrected_graph[u][v][self.cell_flow_attr].sum()

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