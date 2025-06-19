import numpy as np
import flowpaths as fp
import networkx as nx

import flowpaths.utils as utils

from cell_type_tree import CellTypeTree

class MinErrorCellTypeFlow(fp.MinErrorFlow):

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
        self.edge_indexes = [(u, v) for (u, v) in self.G.edges()]
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
            lb=0, 
            ub=self.ub, 
            var_type="integer" if self.weight_type == int else "continuous",
        )

        # Adding flow conservation constraints
        for node in self.G.nodes():
            if self.G.in_degree(node) == 0 or self.G.out_degree(node) == 0:
                continue
            # Flow conservation constraint
            self.solver.add_constraint(
                self.solver.quicksum(
                    self.edge_vars[(u, v)]
                    for (u, v) in self.G.in_edges(node)
                )
                - self.solver.quicksum(
                    self.edge_vars[(u, v)]
                    for (u, v) in self.G.out_edges(node)
                )
                == 0,
                name=f"flow_conservation_{node}",
            )
        
        # Encoding the edge error variables
        for u, v, data in self.G.edges(data=True):
            if (u, v) in self.edges_to_ignore:
                # Making sure the error of the edges to ignore gets set to 0
                self.solver.add_constraint(
                    self.edge_error_vars[(u, v)] == 0,
                    name=f"edge_error_u={u}_v={v}",
                )
                continue
            
            # If the edge is not in the edges_to_ignore list, we need to check if it has a flow attribute
            if self.flow_attr not in data:
                utils.logger.error(f"{__name__}: Flow attribute '{self.flow_attr}' not found in edge data for edge {str((u, v))}, and this edge is not in the edges_to_ignore list.")
                raise ValueError(f"Flow attribute '{self.flow_attr}' not found in edge data for edge {str((u, v))}, and this edge is not in the edges_to_ignore list.")
            
            # Getting the flow value of the edge            
            f_u_v = data[self.flow_attr]
            
            # Encoding the error on the edge (u, v) as the difference between 
            # the flow value of the edge and the sum of the weights of the paths that go through it (pi variables)
            # If we minimize the sum of edge_error_vars, then we are minimizing the sum of the absolute errors.
            self.solver.add_constraint(
                f_u_v - self.edge_vars[(u, v)] <= self.edge_error_vars[(u, v)],
                name=f"edge_error_u={u}_v={v}",
            )

            self.solver.add_constraint(
                self.edge_vars[(u, v)] - f_u_v <= self.edge_error_vars[(u, v)],
                name=f"edge_error_u={u}_v={v}",
            )