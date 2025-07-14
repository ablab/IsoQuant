import time
import networkx as nx
import flowpaths.stdigraph as stdigraph
import flowpaths.utils.graphutils as gu
import flowpaths.abstractpathmodeldag as pathmodel
import flowpaths.utils.safetyflowdecomp as sfd
import flowpaths.utils as utils
import flowpaths.nodeexpandeddigraph as nedg

class kFlowDecompWithSB(pathmodel.AbstractPathModelDAG):
    """
    Class to decompose a flow into a given number of weighted paths.
    """
    # storing some defaults
    optimize_with_greedy = True
    optimize_with_flow_safe_paths = True

    def __init__(
        self,
        G: nx.DiGraph,
        flow_attr: str,
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
        """
        Initialize the Flow Decomposition model for a given number of paths `k`.

        Parameters
        ----------
        - `G : nx.DiGraph`
            
            The input directed acyclic graph, as networkx DiGraph.

        - `flow_attr : str`
            
            The attribute name from where to get the flow values on the edges.

        - `k: int`
            
            The number of paths to decompose in.

        - `flow_attr_origin : str`, optional

            The origin of the flow attribute. Default is `"edge"`. Options:
            
            - `"edge"`: the flow attribute is assumed to be on the edges of the graph.
            - `"node"`: the flow attribute is assumed to be on the nodes of the graph. See [the documentation](node-expanded-digraph.md) on how node-weighted graphs are handled.

        - `weight_type : type`, optional
            
            The type of weights (`int` or `float`). Default is `float`.

        - `subpath_constraints : list`, optional
            
            List of subpath constraints. Default is an empty list. 
            Each subpath constraint is a list of edges that must be covered by some solution path, according 
            to the `subpath_constraints_coverage` or `subpath_constraints_coverage_length` parameters (see below).

        - `subpath_constraints_coverage : float`, optional
            
            Coverage fraction of the subpath constraints that must be covered by some solution paths. 
            
            Defaults to `1.0`, meaning that 100% of the edges (or nodes, if `flow_attr_origin` is `"node"`) of 
            the constraint need to be covered by some solution path). 
            See [subpath constraints documentation](subpath-constraints.md#3-relaxing-the-constraint-coverage)

        - `subpath_constraints_coverage_length : float`, optional
            
            Coverage length of the subpath constraints. Default is `None`. If set, this overrides `subpath_constraints_coverage`, 
            and the coverage constraint is expressed in terms of the subpath constraint length. 
            `subpath_constraints_coverage_length` is then the fraction of the total length of the constraint (specified via `length_attr`) needs to appear in some solution path.
            See [subpath constraints documentation](subpath-constraints.md#3-relaxing-the-constraint-coverage)

        - `length_attr : str`, optional
            
            The attribute name from where to get the edge lengths (or node length, if `flow_attr_origin` is `"node"`). Defaults to `None`.
            
            - If set, then the subpath lengths (above) are in terms of the edge/node lengths specified in the `length_attr` field of each edge/node.
            - If set, and an edge/node has a missing edge length, then it gets length 1.

        - `elements_to_ignore : list`, optional

            List of edges (or nodes, if `flow_attr_origin` is `"node"`) to ignore when adding constrains on flow explanation by the weighted paths. Default is an empty list. See [ignoring edges documentation](ignoring-edges.md)

        - `optimization_options : dict`, optional
            
            Dictionary with the optimization options. Default is `None`. See [optimization options documentation](solver-options-optimizations.md).
            This class also supports the optimization `"optimize_with_greedy": True` (this is the default value). This
            will use a greedy algorithm to solve the problem, and if the number of paths returned by it equals a lowerbound on the solution size,
            then we know the greedy solution is optimum, and it will use that. The lowerbound used currently is the edge-width of the graph,
            meaning the minimum number of paths needed to cover all edges. This is a correct lowerbound because any flow decomposition must cover all edges, 
            as they have non-zero flow.

        - `solver_options : dict`, optional
            
            Dictionary with the solver options. Default is `None`. See [solver options documentation](solver-options-optimizations.md).


        Raises
        ----------
        - ValueError: If `weight_type` is not int or float.
        - ValueError: If some edge does not have the flow attribute specified as `flow_attr`.
        - ValueError: If the graph does not satisfy flow conservation on nodes different from source or sink.
        - ValueError: If the graph contains edges with negative (<0) flow values.
        - ValueError: If `flow_attr_origin` is not "node" or "edge".
        """

        # Handling node-weighted graphs
        self.flow_attr_origin = flow_attr_origin
        if self.flow_attr_origin == "node":
            if G.number_of_nodes() == 0:
                utils.logger.error(f"{__name__}: The input graph G has no nodes. Please provide a graph with at least one node.")
                raise ValueError(f"The input graph G has no nodes. Please provide a graph with at least one node.")
            self.G_internal = nedg.NodeExpandedDiGraph(G, node_flow_attr=flow_attr, node_length_attr=length_attr)
            subpath_constraints_internal = self.G_internal.get_expanded_subpath_constraints(subpath_constraints)
            
            edges_to_ignore_internal = self.G_internal.edges_to_ignore
            if not all(isinstance(element_to_ignore, str) for element_to_ignore in elements_to_ignore):
                utils.logger.error(f"elements_to_ignore must be a list of nodes (i.e strings), not {elements_to_ignore}")
                raise ValueError(f"elements_to_ignore must be a list of nodes (i.e strings), not {elements_to_ignore}")
            edges_to_ignore_internal += [self.G_internal.get_expanded_edge(node) for node in elements_to_ignore]
            edges_to_ignore_internal = list(set(edges_to_ignore_internal))
                
        elif self.flow_attr_origin == "edge":
            if G.number_of_edges() == 0:
                utils.logger.error(f"{__name__}: The input graph G has no edges. Please provide a graph with at least one edge.")
                raise ValueError(f"The input graph G has no edges. Please provide a graph with at least one edge.")
            self.G_internal = G
            subpath_constraints_internal = subpath_constraints
            if not all(isinstance(edge, tuple) and len(edge) == 2 for edge in elements_to_ignore):
                utils.logger.error(f"elements_to_ignore must be a list of edges (i.e. tuples of nodes), not {elements_to_ignore}")
                raise ValueError(f"elements_to_ignore must be a list of edges (i.e. tuples of nodes), not {elements_to_ignore}")
            edges_to_ignore_internal = elements_to_ignore
        else:
            utils.logger.error(f"flow_attr_origin must be either 'node' or 'edge', not {self.flow_attr_origin}")
            raise ValueError(f"flow_attr_origin must be either 'node' or 'edge', not {self.flow_attr_origin}")

        self.G = stdigraph.stDiGraph(self.G_internal)
        self.subpath_constraints = subpath_constraints_internal
        self.edges_to_ignore = self.G.source_sink_edges.union(edges_to_ignore_internal)

        if weight_type not in [int, float]:
            utils.logger.error(f"weight_type must be either int or float, not {weight_type}")
            raise ValueError(f"weight_type must be either int or float, not {weight_type}")
        self.weight_type = weight_type

        # Check requirements on input graph:
        # Check flow conservation only if there are no edges to ignore
        satisfies_flow_conservation = gu.check_flow_conservation(G, flow_attr)
        if len(edges_to_ignore_internal) == 0 and not satisfies_flow_conservation:
            utils.logger.error(f"{__name__}: The graph G does not satisfy flow conservation or some edges have missing `flow_attr`. This is an error, unless you passed `edges_to_ignore` to include at least those edges with missing `flow_attr`.")
            raise ValueError("The graph G does not satisfy flow conservation or some edges have missing `flow_attr`. This is an error, unless you passed `edges_to_ignore` to include at least those edges with missing `flow_attr`.")

        # Check that the flow is positive and get max flow value
        self.flow_attr = flow_attr
        self.w_max = self.weight_type(
            self.G.get_max_flow_value_and_check_non_negative_flow(
                flow_attr=self.flow_attr, edges_to_ignore=self.edges_to_ignore
            )
        )

        self.k = k
        
        self.subpath_constraints_coverage = subpath_constraints_coverage
        self.subpath_constraints_coverage_length = subpath_constraints_coverage_length
        self.length_attr = length_attr

        self.pi_vars = {}
        self.path_weights_vars = {}

        self.path_weights_sol = None
        self._solution = None
        self._lowerbound_k = None
        
        self.solve_statistics = {}
        self.optimization_options = optimization_options.copy() or {}

        greedy_solution_paths = None
        self.optimize_with_greedy = self.optimization_options.get("optimize_with_greedy", kFlowDecompWithSB.optimize_with_greedy)
        self.optimize_with_flow_safe_paths = self.optimization_options.get("optimize_with_flow_safe_paths", kFlowDecompWithSB.optimize_with_flow_safe_paths)
        
        # We can apply the greedy algorithm only if 
        # - there are no edges to ignore (in the original input graph), and 
        # - the graph satisfies flow conservation
        if self.optimize_with_greedy and len(edges_to_ignore_internal) == 0 and satisfies_flow_conservation:
            if self._get_solution_with_greedy():
                greedy_solution_paths = self._solution["paths"]
                self.optimization_options["external_solution_paths"] = greedy_solution_paths
        
        if self.optimize_with_flow_safe_paths and satisfies_flow_conservation:
            start_time = time.perf_counter()
            self.optimization_options["external_safe_paths"] = sfd.compute_flow_decomp_safe_paths(G=G, flow_attr=self.flow_attr)
            self.solve_statistics["flow_safe_paths_time"] = time.perf_counter() - start_time
            # If we optimize with flow safe paths, we need to disable optimizing with safe paths and sequences
            if self.optimization_options.get("optimize_with_safe_paths", False):
                utils.logger.error(f"{__name__}: Cannot optimize with both flow safe paths and safe paths")
                raise ValueError("Cannot optimize with both flow safe paths and safe paths")
            if self.optimization_options.get("optimize_with_safe_sequences", False):
                utils.logger.error(f"{__name__}: Cannot optimize with both flow safe paths and safe sequences")
                raise ValueError("Cannot optimize with both flow safe paths and safe sequences")
        
        self.optimization_options["trusted_edges_for_safety"] = self.G.get_non_zero_flow_edges(flow_attr=self.flow_attr, edges_to_ignore=self.edges_to_ignore)

        # Call the constructor of the parent class AbstractPathModelDAG
        super().__init__(
            G=self.G, 
            k=self.k,
            subpath_constraints=self.subpath_constraints, 
            subpath_constraints_coverage=self.subpath_constraints_coverage, 
            subpath_constraints_coverage_length=self.subpath_constraints_coverage_length,
            length_attr=self.length_attr, 
            optimization_options=self.optimization_options,
            solver_options=solver_options,
            solve_statistics=self.solve_statistics,
        )

        # If already solved with a previous method, we don't create solver, not add paths
        if self.is_solved():
            return

        # This method is called from the super class AbstractPathModelDAG
        self.create_solver_and_paths()

        # This method is called from the current class to encode the flow decomposition
        self._encode_flow_decomposition()

        # The given weights optimization
        self._encode_given_weights()

        utils.logger.info(f"{__name__}: initialized with graph id = {utils.fpid(G)}, k = {self.k}")

    def _encode_flow_decomposition(self):
        
        # Encodes the flow decomposition constraints for the given graph.
        # This method sets up the path weight variables and the edge variables encoding
        # the sum of the weights of the paths going through the edge.

        # If already solved, no need to encode further
        if self.is_solved():
            return

        # pi vars from https://arxiv.org/pdf/2201.10923 page 14
        self.pi_vars = self.solver.add_variables(
            self.edge_indexes,
            name_prefix="pi",
            lb=0,
            ub=self.w_max,
            var_type="integer" if self.weight_type == int else "continuous",
        )
        self.path_weights_vars = self.solver.add_variables(
            self.path_indexes,
            name_prefix="w",
            lb=0,
            ub=self.w_max,
            var_type="integer" if self.weight_type == int else "continuous",
        )

        self.symmetry_vars = self.solver.add_variables(
            self.path_indexes,
            name_prefix = "u",
            lb = 0,
            ub = 1,
            var_type = "integer",
        )

        # Symmetry breaking constraints by Reima Kuosmanen
        for i in range(len(self.safe_lists), self.k):

            edge_list = self.G.edges()
            m = len(edge_list)

            u, v = edge_list[0]

            self.solver.add_constraint(
                self.edge_vars[(u, v, i)] - self.edge_vars[(u, v, i - 1)] == self.symmetry_vars[(u, v, i)]
            )

            for p in range(1, m):
                
                u, v = edge_list[p]

                # Has increased
                self.solver.add_constraint(
                    self.edge_vars[(u, v, i)] - self.edge_vars[(u, v, i - 1)] >= self.symmetry_vars[(u, v, i)]
                )

                # increase means the index can lower
                self.solver.add_constraint(
                    self.edge_vars[(u, v, i)] - self.edge_vars[(u, v, i - 1)] >= 
                    self.solver.quicksum(
                        -self.symmetry_vars[(u0, v0, i)]
                        for u0, v0 in edge_list[0:p]
                    )
                )


        # We encode that for each edge (u,v), the sum of the weights of the paths going through the edge is equal to the flow value of the edge.
        for u, v, data in self.G.edges(data=True):
            if (u, v) in self.edges_to_ignore:
                continue
            f_u_v = data[self.flow_attr]

            # We encode that edge_vars[(u,v,i)] * self.path_weights_vars[(i)] = self.pi_vars[(u,v,i)],
            # assuming self.w_max is a bound for self.path_weights_vars[(i)]
            for i in range(self.k):
                self.solver.add_binary_continuous_product_constraint(
                    binary_var=self.edge_vars[(u, v, i)],
                    continuous_var=self.path_weights_vars[(i)],
                    product_var=self.pi_vars[(u, v, i)],
                    lb=0,
                    ub=self.w_max,
                    name=f"10_u={u}_v={v}_i={i}",
                )

            self.solver.add_constraint(
                self.solver.quicksum(self.pi_vars[(u, v, i)] for i in range(self.k)) == f_u_v,
                name=f"10d_u={u}_v={v}_i={i}",
            )

    def _encode_given_weights(self):

        weights = self.optimization_options.get("given_weights", None)
        if weights is None:
            return
        
        if self.optimization_options.get("optimize_with_safe_paths", False):
            utils.logger.error(f"{__name__}: Cannot optimize with both given weights and safe paths")
            raise ValueError("Cannot optimize with both given weights and safe paths")
        if self.optimization_options.get("optimize_with_safe_sequences", False):
            utils.logger.error(f"{__name__}: Cannot optimize with both given weights and safe sequences")
            raise ValueError("Cannot optimize with both given weights and safe sequences")
        if self.optimization_options.get("optimize_with_safe_zero_edges", False):
            utils.logger.error(f"{__name__}: Cannot optimize with both given weights and safe zero edges")
            raise ValueError("Cannot optimize with both given weights and safe zero edges")
        if self.optimization_options.get("optimize_with_flow_safe_paths", False):
            utils.logger.error(f"{__name__}: Cannot optimize with both given weights and flow safe paths")
            raise ValueError("Cannot optimize with both given weights and flow safe paths")
            
        if len(weights) > self.k:
            utils.logger.error(f"Length of given weights ({len(weights)}) is greater than k ({self.k})")
            raise ValueError(f"Length of given weights ({len(weights)}) is greater than k ({self.k})")

        for i, weight in enumerate(weights):
            self.solver.add_constraint(
                self.path_weights_vars[i] == weight,
                name=f"given_weight_{i}",
            )

        self.solver.set_objective(
            self.solver.quicksum(self.edge_vars[(u, v, i)] for u, v in self.G.edges() for i in range(self.k)),
            sense="minimize",
        )

    def _get_solution_with_greedy(self):
        
        # Attempts to find a solution using a greedy algorithm.
        # This method first decomposes the problem using the maximum bottleneck approach.
        # If the number of paths obtained is less than or equal to the specified limit `k`,
        # it sets the solution and marks the problem as solved. It also records the time
        # taken to solve the problem using the greedy approach.

        # Returns
        # -------
        # - bool: True if a solution is found using the greedy algorithm, False otherwise.
        

        start_time = time.perf_counter()
        (paths, weights) = self.G.decompose_using_max_bottleneck(self.flow_attr)

        # Check if the greedy decomposition satisfies the subpath constraints
        if self.subpath_constraints:
            for subpath in self.subpath_constraints:
                if self.subpath_constraints_coverage_length is None:
                    # By default, the length of the constraints is its number of edges 
                    constraint_length = len(subpath)
                    # And the fraction of edges that we need to cover is self.subpath_constraints_coverage
                    coverage_fraction = self.subpath_constraints_coverage
                else:
                    constraint_length = sum(self.G[u][v].get(self.length_attr, 1) for (u,v) in subpath)
                    coverage_fraction = self.subpath_constraints_coverage_length
                # If the subpath is not covered enough by the greedy decomposition, we return False
                if gu.max_occurrence(subpath, paths, edge_lengths={(u,v): self.G[u][v].get(self.length_attr, 1) for (u,v) in subpath}) < constraint_length * coverage_fraction:
                    return False
        
        if len(paths) <= self.k:
            # If paths contains strictly less than self.k paths, 
            # then we add arbitrary paths (i.e. we repeat the first path) with 0 weights to reach self.k paths.
            paths += [paths[0] for _ in range(self.k - len(paths))]
            weights += [0 for _ in range(self.k - len(weights))]
            # self._solution = {
            #     "paths": paths,
            #     "weights": weights,
            # }
            if self.flow_attr_origin == "edge":
                self._solution = {
                    "paths": paths,
                    "weights": weights,
                }
            elif self.flow_attr_origin == "node":
                self._solution = {
                    "_paths_internal": paths,
                    "paths": self.G_internal.get_condensed_paths(paths),
                    "weights": self.path_weights_sol,
                }
            self.set_solved()
            self.solve_statistics = {}
            self.solve_statistics["greedy_solve_time"] = time.perf_counter() - start_time
            return True

        return False

    def _remove_empty_paths(self, solution):
        """
        Removes empty paths from the solution. Empty paths are those with 0 or 1 nodes.

        Parameters
        ----------
        - `solution: dict`
            
            The solution dictionary containing paths and weights.

        Returns
        -------
        - `solution: dict`
            
            The solution dictionary with empty paths removed.

        """
        non_empty_paths = []
        non_empty_weights = []
        for path, weight in zip(solution["paths"], solution["weights"]):
            if len(path) > 1:
                non_empty_paths.append(path)
                non_empty_weights.append(weight)
        return {"paths": non_empty_paths, "weights": non_empty_weights}
    


    def get_solution(self, remove_empty_paths=False):
        """
        Retrieves the solution for the flow decomposition problem.

        If the solution has already been computed and cached as `self.solution`, it returns the cached solution.
        Otherwise, it checks if the problem has been solved, computes the solution paths and weights,
        and caches the solution.

        Parameters
        ----------

        - `remove_empty_paths: bool`, optional

            If `True`, removes empty paths from the solution. Default is `False`. These can happen only if passed the optimization option `"allow_empty_paths" : True`.

        Returns
        -------
        - `solution: dict`
        
            A dictionary containing the solution paths (key `"paths"`) and their corresponding weights (key `"weights"`).

        Raises
        ------
        - `exception` If model is not solved.
        """

        if self._solution is None:            

            self.check_is_solved()
            weights_sol_dict = self.solver.get_variable_values("w", [int])
            self.path_weights_sol = [
                (
                    round(weights_sol_dict[i])
                    if self.weight_type == int
                    else float(weights_sol_dict[i])
                )
                for i in range(self.k)
            ]

            if self.flow_attr_origin == "edge":
                self._solution = {
                    "paths": self.get_solution_paths(),
                    "weights": self.path_weights_sol,
                }
            elif self.flow_attr_origin == "node":
                self._solution = {
                    "_paths_internal": self.get_solution_paths(),
                    "paths": self.G_internal.get_condensed_paths(self.get_solution_paths()),
                    "weights": self.path_weights_sol,
                }

        return self._remove_empty_paths(self._solution) if remove_empty_paths else self._solution

    def is_valid_solution(self, tolerance=0.001):
        """
        Checks if the solution is valid by comparing the flow from paths with the flow attribute in the graph edges.

        Raises
        ------
        - ValueError: If the solution is not available (i.e., self.solution is None).

        Returns
        -------
        - bool: True if the solution is valid, False otherwise.

        Notes
        -------
        - get_solution() must be called before this method.
        - The solution is considered valid if the flow from paths is equal
            (up to `TOLERANCE * num_paths_on_edges[(u, v)]`) to the flow value of the graph edges.
        """

        if self._solution is None:
            utils.logger.error(f"{__name__}: Solution is not available. Call get_solution() first.")
            raise ValueError("Solution is not available. Call get_solution() first.")

        solution_paths = self._solution.get("_paths_internal", self._solution["paths"])
        solution_weights = self._solution["weights"]
        solution_paths_of_edges = [
            [(path[i], path[i + 1]) for i in range(len(path) - 1)]
            for path in solution_paths
        ]

        flow_from_paths = {(u, v): 0 for (u, v) in self.G.edges()}
        num_paths_on_edges = {e: 0 for e in self.G.edges()}
        for weight, path in zip(solution_weights, solution_paths_of_edges):
            for e in path:
                flow_from_paths[e] += weight
                num_paths_on_edges[e] += 1

        for u, v, data in self.G.edges(data=True):
            if self.flow_attr in data and (u,v) not in self.edges_to_ignore:
                if (
                    abs(flow_from_paths[(u, v)] - data[self.flow_attr])
                    > tolerance * num_paths_on_edges[(u, v)]
                ):
                    utils.logger.debug(f"Flow validation failed for edge ({u}, v): expected {data[self.flow_attr]}, got {flow_from_paths[(u, v)]}")
                    return False

        return True
    
    def get_objective_value(self):
        
        self.check_is_solved()

        if self._solution is None:
            self.get_solution()

        return self.k
    
    def get_lowerbound_k(self):

        if self._lowerbound_k != None:
            return self._lowerbound_k

        self._lowerbound_k = self.G.get_width(edges_to_ignore=self.edges_to_ignore)

        return self._lowerbound_k
