import numpy as np
import flowpaths as fp
import networkx as nx

from src.cell_type_tree import CellTypeTree

class kFlowCellTypeDecomp(fp.kFlowDecomp):

    optimize_with_greedy = True
    optimize_with_flow_safe_paths = True

    def __init__(
            self,
            G: nx.DiGraph,
            cell_tree: CellTypeTree,
            flow_attr: str,
            cell_flow_attr: str,
            k: int,
            flow_attr_origin: str = "edge",
            weight_type: type = float,
            subpath_constraints: list = [],
            subpath_constraints_coverage: float = 1.0,
            subpath_constraints_coverage_length: float = None,
            length_attr: str = None,
            elements_to_ignore: list = [],
            optimization_options: dict = {},
            solver_options: dict = {},
        ):

        self.cell_tree = cell_tree
        self.cell_types = cell_tree.get_leaf_types()
        self.cell_flow_attr = cell_flow_attr

        super().__init__(
            G = G,
            flow_attr = flow_attr,
            k = k,
            flow_attr_origin = flow_attr_origin,
            weight_type = weight_type,
            subpath_constraints = subpath_constraints,
            subpath_constraints_coverage = subpath_constraints_coverage,
            subpath_constraints_coverage_length = subpath_constraints_coverage_length,
            length_attr = length_attr,
            elements_to_ignore = elements_to_ignore,
            optimization_options = optimization_options,
            solver_options = solver_options,
        )

        
    def _encode_flow_decomposition(self):

        cell_types = self.cell_tree.get_leaf_types()
        normal = False
        
        if self.is_solved():
            return

        self.pi_vars = self.solver.add_variables(
            [(u, v, i, ct) for u, v, i in self.edge_indexes for ct in cell_types],
            name_prefix = "pc",
            lb = 0,
            ub = self.w_max,
            var_type="integer" if self.weight_type == int else "continuous",
        )

        self.path_weights_vars = self.solver.add_variables(
            [(i, ct) for i in self.path_indexes for ct in cell_types],
            name_prefix="v",
            lb=0,
            ub=self.w_max,
            var_type="integer"# if self.weight_type == int else "continuous",
        )

        # We encode that for each edge (u,v), the sum of the weights of the paths going through the edge is equal to the flow value of the edge.
        for u, v, data in self.G.edges(data=True):
            
            if (u, v) in self.edges_to_ignore:
                continue
            
            #f_u_v = data[self.flow_attr]
            ct_u_v = data[self.cell_flow_attr]

            # We encode that edge_vars[(u,v,i)] * self.path_weights_vars[(i)] = self.pi_vars[(u,v,i)],
            # assuming self.w_max is a bound for self.path_weights_vars[(i)]
            for i in range(self.k):
                
                # ct specific variables
                for ct in cell_types:
                    self.solver.add_binary_continuous_product_constraint(
                        binary_var = self.edge_vars[(u, v, i)],
                        continuous_var = self.path_weights_vars[(i, ct)],
                        product_var = self.pi_vars[(u, v, i, ct)],
                        lb = 0,
                        ub = self.w_max,
                        name=f"10ct_u={u}_v={v}_i={i}_ct={ct}",
                    )

            # Total flow of cell in a edge
            for cell_type in cell_types:
                self.solver.add_constraint(
                    self.solver.quicksum(self.pi_vars[(u, v, i, cell_type)] for i in range(self.k)) == int(ct_u_v[self.cell_tree.get_cell_type_index(cell_type)]),
                    name=f"tcf_u={u}_v={v}_ct={cell_type}"
            )
            
    def get_solution(self, remove_empty_paths=False):

        if self._solution is None:            

            self.check_is_solved()
            #weights_sol_dict = self.solver.get_variable_values("w", [int])
            weights_sol_ct_dict = self.solver.get_variable_values("v", [int, str])
            #print(weights_sol_ct_dict)
            #print(weights_sol_ct_dict)
            #self.path_weights_sol = [
            #    (
            #        round(weights_sol_dict[i])
            #        if self.weight_type == self.weight_type#int
            #        else float(weights_sol_dict[i])
            #    )
            #    for i in range(self.k)
            #]
            self.path_weights_ct_sol = self.cell_tree.transform_count_to_arrays(weights_sol_ct_dict)
            #print(self.path_weights_ct_sol)
            #for i, ct in weights_sol_ct_dict:
            #    pass #"self.path_weights_ct_sol[i][self.cell_tree.get_cell_type_index(ct)] = weights_sol_ct_dict[(i, ct)]
    
            
            self.path_weights_sol = self.path_weights_ct_sol.sum(axis = 1)

            if self.flow_attr_origin == "edge":
                self._solution = {
                    "paths": self.get_solution_paths(),
                    "weights": self.path_weights_sol,
                    "ct_weights": self.path_weights_ct_sol,
                }
            elif self.flow_attr_origin == "node":
                self._solution = {
                    "_paths_internal": self.get_solution_paths(),
                    "paths": self.G_internal.get_condensed_paths(self.get_solution_paths()),
                    "weights": self.path_weights_sol,
                    "ct_weights": self.path_weights_ct_sol,
                }

        return self._remove_empty_paths(self._solution) if remove_empty_paths else self._solution