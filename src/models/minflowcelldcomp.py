import time
import copy

import numpy as np
import flowpaths as fp
import networkx as nx

import flowpaths.utils.solverwrapper as sw
import flowpaths.utils as utils

from src.cell_type_tree import CellTypeTree
from .kflowcellcomp import kFlowCellTypeDecomp


class MinFlowCellDecomp(fp.MinFlowDecomp):

    subgraph_lowerbound_size = 20
    subgraph_lowerbound_shift = 18
    optimize_with_given_weights_num_free_paths = 0
    use_min_gen_set_lowerbound = False
    use_min_gen_set_lowerbound_partition_constraints = False
    use_min_gen_set_lowerbound_partition_constraints_min_constraint_len = 2
    use_min_gen_set_lowerbound_partition_constraints_limit_num_constraints = 3
    min_gen_set_remove_sums_of_two = True
    optimize_with_given_weights = False
    use_subgraph_scanning_lowerbound = False
    use_subgraph_scanning_weights_in_given_weights_optimization = True

    def __init__(
            self,
            G: nx.DiGraph,
            cell_tree: CellTypeTree,
            flow_attr: str,
            cell_flow_attr: str,
            flow_attr_origin: str = "edge",
            weight_type: type = float,
            subpath_constraints: list = [],
            subpath_constraints_coverage: float = 1.0,
            subpath_constraints_coverage_length: float = None,
            length_attr: str = None,
            elements_to_ignore: list = [],
            additional_starts: list = [],
            additional_ends: list = [],
            optimization_options: dict = {},
            solver_options: dict = {},
            ):
        
        super().__init__(
            G = G,
            flow_attr = flow_attr,
            flow_attr_origin = flow_attr_origin,
            weight_type = weight_type,
            subpath_constraints = subpath_constraints,
            subpath_constraints_coverage = subpath_constraints_coverage,
            subpath_constraints_coverage_length = subpath_constraints_coverage_length,
            length_attr = length_attr,
            elements_to_ignore = elements_to_ignore,
            additional_starts = additional_starts,
            additional_ends = additional_ends,
            optimization_options = optimization_options,
            solver_options = solver_options,
            )
        
        self.cell_tree = cell_tree
        self.cell_flow_attr = cell_flow_attr
        self.cell_types = cell_tree.get_leaf_types()

    def solve(self) -> bool:

        self.solve_time_start = time.perf_counter()

        if self.optimization_options.get("optimize_with_given_weights", MinFlowCellDecomp.optimize_with_given_weights):            
            self._solve_with_given_weights()

        for i in range(self.get_lowerbound_k(), self.G.number_of_edges()):
            utils.logger.info(f"{__name__}: iteration with k = {i}")
            fd_model = None
            # Checking if we have already found a solution with the same number of paths
            # via the min gen set and given weights approach
            if self._given_weights_model is not None and self._given_weights_model.is_solved():
                if len(self._given_weights_model.get_solution(remove_empty_paths=True)["paths"]) == i:
                    fd_model = self._given_weights_model

            fd_solver_options = copy.deepcopy(self.solver_options)
            if "time_limit" in fd_solver_options:
                fd_solver_options["time_limit"] = self.time_limit - self.solve_time_elapsed

            if fd_model is None:
                fd_model = kFlowCellTypeDecomp(
                    G=self.G,
                    cell_tree=self.cell_tree,
                    flow_attr=self.flow_attr,
                    cell_flow_attr=self.cell_flow_attr,
                    k=i,
                    weight_type=self.weight_type,
                    subpath_constraints=self.subpath_constraints,
                    subpath_constraints_coverage=self.subpath_constraints_coverage,
                    subpath_constraints_coverage_length=self.subpath_constraints_coverage_length,
                    length_attr=self.length_attr,
                    elements_to_ignore=self.edges_to_ignore,
                    optimization_options=self.optimization_options,
                    solver_options=fd_solver_options,
                )
                fd_model.solve()

            if fd_model.is_solved():
                self._solution = fd_model.get_solution(remove_empty_paths=False)
                #print(self._solution)
                if self.flow_attr_origin == "node":
                    # If the flow_attr_origin is "node", we need to convert the solution paths from the expanded graph to paths in the original graph.
                    self._solution["_paths_internal"] = self._solution["paths"]
                    self._solution["paths"] = self.G_internal.get_condensed_paths(self._solution["paths"])
                    #self._solution["ct_weights"] = self
                self.set_solved()
                self.solve_statistics = fd_model.solve_statistics
                self.solve_statistics["mfd_solve_time"] = time.perf_counter() - self.solve_time_start
                self.fd_model = fd_model
                return True
            elif fd_model.solver.get_model_status() != sw.SolverWrapper.infeasible_status:
                # If the model is not solved and the status is not infeasible,
                # it means that the solver stopped because of an unexpected termination,
                # thus we cannot conclude that the model is infeasible.
                # In this case, we stop the search.
                return False

        return False