import time
import networkx as nx
import flowpaths.stdigraph as stdigraph
from .kflowdecompSB import kFlowDecompWithSB
import flowpaths.abstractpathmodeldag as pathmodel
import flowpaths.utils.solverwrapper as sw
import flowpaths.utils.graphutils as gu
import flowpaths.mingenset as mgs
import flowpaths.utils as utils
import flowpaths.nodeexpandeddigraph as nedg
import copy
import math

class MinFlowDecompSB(pathmodel.AbstractPathModelDAG): # Note that we inherit from AbstractPathModelDAG to be able to use this class to also compute safe paths. 
    """
    A class to decompose a network flow if a directed acyclic graph into a minimum number of weighted paths.
    """

    # Default optimization parameters
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
        flow_attr: str,
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
        """
        Initialize the Minimum Flow Decomposition model, minimizing the number of paths.

        Parameters
        ----------
        - `G : nx.DiGraph`
            
            The input directed acyclic graph, as networkx DiGraph.

        - `flow_attr : str`
            
            The attribute name from where to get the flow values on the edges.

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

        - `subpath_constraints_coverage: float`, optional
            
            Coverage fraction of the subpath constraints that must be covered by some solution paths. 
            
            Defaults to `1.0`, meaning that 100% of the edges (or nodes, if `flow_attr_origin` is `"node"`) of 
            the constraint need to be covered by some solution path). 
            See [subpath constraints documentation](subpath-constraints.md#3-relaxing-the-constraint-coverage)

        - `subpath_constraints_coverage_length : float`, optional
            
            Coverage length of the subpath constraints. Default is `None`. If set, this overrides `subpath_constraints_coverage`, 
            and the coverage constraint is expressed in terms of the subpath constraint length. 
            `subpath_constraints_coverage_length` is then the fraction of the total length of the constraint (specified via `length_attr`) needs to appear in some solution path.
            See [subpath constraints documentation](subpath-constraints.md#3-relaxing-the-constraint-coverage)

        - `length_attr: str`, optional
            
            The attribute name from where to get the edge lengths (or node length, if `flow_attr_origin` is `"node"`). Defaults to `None`.
            
            - If set, then the subpath lengths (above) are in terms of the edge/node lengths specified in the `length_attr` field of each edge/node.
            - If set, and an edge/node has a missing edge length, then it gets length 1.

        - `elements_to_ignore : list`, optional

            List of edges (or nodes, if `flow_attr_origin` is `"node"`) to ignore when adding constrains on flow explanation by the weighted paths. 
            Default is an empty list. See [ignoring edges documentation](ignoring-edges.md)

        - `additional_starts: list`, optional
            
            List of additional start nodes of the paths. Default is an empty list. See [additional start/end nodes documentation](additional-start-end-nodes.md). **You can set this only if `flow_attr_origin` is `"node"`**.

        - `additional_ends: list`, optional
            
            List of additional end nodes of the paths. Default is an empty list. See [additional start/end nodes documentation](additional-start-end-nodes.md). **You can set this only if `flow_attr_origin` is `"node"`**.

        - `optimization_options : dict`, optional
            
            Dictionary with the optimization options. Default is an empty dict. See [optimization options documentation](solver-options-optimizations.md).
            This class also supports the optimization `"optimize_with_greedy": True` (this is the default value). This
            will use a greedy algorithm to solve the problem, and if the number of paths returned by it equals a lowerbound on the solution size,
            then we know the greedy solution is optimum, and it will use that. The lowerbound used currently is the edge-width of the graph,
            meaning the minimum number of paths needed to cover all edges. This is a correct lowerbound because any flow decomposition must cover all edges, 
            as they have non-zero flow.

        - `solver_options : dict`, optional
            
            Dictionary with the solver options. Default is `{}`. See [solver options documentation](solver-options-optimizations.md).

        Raises
        ------
        `ValueError`

        - If `weight_type` is not `int` or `float`.
        - If some edge does not have the flow attribute specified as `flow_attr`.
        - If the graph does not satisfy flow conservation on nodes different from source or sink.
        - If the graph contains edges with negative (<0) flow values.
        - If the graph is not acyclic.
        - If `flow_attr_origin` is not "node" or "edge".
        """

        # Handling node-weighted graphs
        self.flow_attr_origin = flow_attr_origin
        if self.flow_attr_origin == "node":
            if G.number_of_nodes() == 0:
                utils.logger.error(f"{__name__}: The input graph G has no nodes. Please provide a graph with at least one node.")
                raise ValueError(f"The input graph G has no nodes. Please provide a graph with at least one node.")
            if len(additional_starts) + len(additional_ends) == 0:
                self.G_internal = nedg.NodeExpandedDiGraph(
                    G=G, 
                    node_flow_attr=flow_attr, 
                    node_length_attr=length_attr)
            else:
                self.G_internal = nedg.NodeExpandedDiGraph(
                    G=G, 
                    node_flow_attr=flow_attr, 
                    node_length_attr=length_attr,
                    try_filling_in_missing_flow_attr=True,
                    additional_starts=additional_starts,
                    additional_ends=additional_ends,
                )
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
            if len(additional_starts) + len(additional_ends) > 0:
                utils.logger.error(f"additional_starts and additional_ends are not supported when flow_attr_origin is 'edge'.")
                raise ValueError(f"additional_starts and additional_ends are not supported when flow_attr_origin is 'edge'.")
            self.G_internal = G
            subpath_constraints_internal = subpath_constraints
            if not all(isinstance(edge, tuple) and len(edge) == 2 for edge in elements_to_ignore):
                utils.logger.error(f"elements_to_ignore must be a list of edges (i.e. tuples of nodes), not {elements_to_ignore}")
                raise ValueError(f"elements_to_ignore must be a list of edges (i.e. tuples of nodes), not {elements_to_ignore}")
            edges_to_ignore_internal = elements_to_ignore

        else:
            utils.logger.error(f"flow_attr_origin must be either 'node' or 'edge', not {self.flow_attr_origin}")
            raise ValueError(f"flow_attr_origin must be either 'node' or 'edge', not {self.flow_attr_origin}")

        self.G = self.G_internal
        self.subpath_constraints = subpath_constraints_internal
        self.edges_to_ignore = edges_to_ignore_internal
        
        self.flow_attr = flow_attr
        self.weight_type = weight_type
        self.subpath_constraints_coverage = subpath_constraints_coverage
        self.subpath_constraints_coverage_length = subpath_constraints_coverage_length
        self.length_attr = length_attr
        self.optimization_options = optimization_options
        self.solver_options = solver_options
        self.time_limit = self.solver_options.get("time_limit", sw.SolverWrapper.time_limit)
        self.solve_time_start = None

        self.solve_statistics = {}
        self._solution = None
        self._lowerbound_k = None
        self._is_solved = None

        # Internal variables
        self._generating_set = None
        self._all_subgraph_weights = None
        self._given_weights_model = None
        self._source_flow = None

        utils.logger.info(f"{__name__}: initialized with graph id = {utils.fpid(G)}")

    def solve(self) -> bool:
        """
        Attempts to solve the flow decomposition problem using a model with varying number of paths.

        This method iterates over a range of possible path numbers, creating and solving a flow decomposition model for each count.
        If a solution is found, it stores the solution and relevant statistics, and returns True. If no solution is found after
        iterating through all possible path counts, it returns False.

        Returns:
            bool: True if a solution is found, False otherwise.

        Note:
            This overloads the `solve()` method from `AbstractPathModelDAG` class.
        """
        self.solve_time_start = time.perf_counter()

        if self.optimization_options.get("optimize_with_given_weights", MinFlowDecompSB.optimize_with_given_weights):            
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
                fd_model = kFlowDecompWithSB(
                    G=self.G,
                    flow_attr=self.flow_attr,
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
                self._solution = fd_model.get_solution(remove_empty_paths=True)
                if self.flow_attr_origin == "node":
                    # If the flow_attr_origin is "node", we need to convert the solution paths from the expanded graph to paths in the original graph.
                    self._solution["_paths_internal"] = self._solution["paths"]
                    self._solution["paths"] = self.G_internal.get_condensed_paths(self._solution["paths"])
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

    def _solve_with_given_weights(self) -> bool:

        all_weights = set({self.G.edges[e][self.flow_attr] for e in self.G.edges() if self.flow_attr in self.G.edges[e]})
        all_weights_list = list(all_weights)
        
        # We call this so that the generating set is computed and stored in the class, if this optimization is activated
        _ = self.get_lowerbound_k()

        if self._generating_set is not None:
            all_weights.update(self._generating_set)
            all_weights_list = list(all_weights)

        # print("all_weights_list", sorted(all_weights_list))

        if self.optimization_options.get("use_subgraph_scanning_weights_in_given_weights_optimization", MinFlowDecompSB.use_subgraph_scanning_weights_in_given_weights_optimization):
            if self._all_subgraph_weights is not None:
                all_weights.update(self._all_subgraph_weights)
                all_weights_list = list(all_weights)
        
        # print("all_weights_list", sorted(all_weights_list))

        given_weights_optimization_options = copy.deepcopy(self.optimization_options)
        given_weights_optimization_options["optimize_with_greedy"] = False
        given_weights_optimization_options["optimize_with_safe_paths"] = False
        given_weights_optimization_options["optimize_with_safe_sequences"] = False
        given_weights_optimization_options["optimize_with_zero_safe_edges"] = False
        given_weights_optimization_options["optimize_with_flow_safe_paths"] = False
        given_weights_optimization_options["allow_empty_paths"] = True
        given_weights_optimization_options["given_weights"] = all_weights_list
        utils.logger.info(f"{__name__}: Solving with given weights = {given_weights_optimization_options['given_weights']}")

        given_weights_kfd_solver_options = copy.deepcopy(self.solver_options)
        if "time_limit" in given_weights_kfd_solver_options:
            given_weights_kfd_solver_options["time_limit"] = self.time_limit - self.solve_time_elapsed

        given_weights_kfd_solver = kFlowDecompWithSB(
            G=self.G,
            k = len(given_weights_optimization_options["given_weights"]) + self.optimization_options.get("optimize_with_given_weights_num_free_paths", MinFlowDecompSB.optimize_with_given_weights_num_free_paths),
            flow_attr=self.flow_attr,
            weight_type=self.weight_type,
            subpath_constraints=self.subpath_constraints,
            subpath_constraints_coverage=self.subpath_constraints_coverage,
            subpath_constraints_coverage_length=self.subpath_constraints_coverage_length,
            length_attr=self.length_attr,
            elements_to_ignore=self.edges_to_ignore,
            optimization_options=given_weights_optimization_options,
            solver_options=given_weights_kfd_solver_options,
            )
        given_weights_kfd_solver.solve()

        if given_weights_kfd_solver.is_solved():
            self._given_weights_model = given_weights_kfd_solver
            sol = self._given_weights_model.get_solution(remove_empty_paths=True)
            utils.logger.info(f"{__name__}: found an MFD solution with given weights in {len(sol['paths'])} paths weights {sol['weights']}")
        else:
            utils.logger.info(f"{__name__}: did NOT found an MFD solution with given weights")

    def _get_source_flow(self):
        if self._source_flow is None:
            self._source_flow = 0
            for v in self.G.nodes():
                if self.G.in_degree(v) == 0:
                    for _, _, data in self.G.out_edges(v, data=True):
                        if self.flow_attr in data:
                            self._source_flow += data[self.flow_attr]
            utils.logger.debug(f"{__name__}: source_flow = {self._source_flow}")
            return self._source_flow
        else:
            return self._source_flow

    def _get_partition_constraints_for_min_gen_set(
            self, 
            min_constraint_len: int = 1, 
            limit_num_constraints: int = None) -> list:
        """
        Get the partition constraints for the min gen set problem. 

        Returns
        -------
        - `partition_constraints: list`
        
            A list of partition constraints, where each constraint is a list of numbers, 
            such that numbers of each constraint sum up to the flow value.
        """

        partition_constraints = set()

        # Sources get level 0, nodes with in-neighbors only sources get level 1, and so on
        level = dict()
        max_level = 0
        for node in nx.topological_sort(self.G):
            if self.G.in_degree(node) == 0:
                level[node] = 0
            else:
                level[node] = max([level[n] for n in self.G.predecessors(node)]) + 1
                max_level = max(max_level, level[node])

        level_edges = [[] for _ in range(max_level)]

        # An edge (u,v) is at level i if level[u] <= i and level[v] >= i+1
        for u, v in self.G.edges():
            for i in range(level[u], level[v]):
                level_edges[i].append((u, v))

        source_flow = self._get_source_flow()

        # Now we create the partition constraints
        for i in range(max_level):
            
            level_flow_sum = 0
            level_flow_parts = []
            for u, v in level_edges[i]:
                if (u, v) in self.edges_to_ignore or self.flow_attr not in self.G.edges[u, v]:
                    continue
                level_flow_sum += self.G.edges[u, v][self.flow_attr]
                level_flow_parts.append(self.G.edges[u, v][self.flow_attr])
            
            if level_flow_sum == source_flow and len(level_flow_parts) >= min_constraint_len:
                # We add the constraint for this level
                partition_constraints.add(tuple(sorted(level_flow_parts)))

        partition_constraints_list = [list(constraint) for constraint in partition_constraints]
        
        partition_constraints_list = sorted(partition_constraints_list, key=lambda x: min(x), reverse=False)

        # We remove the constraints whose min value equals the min value of the previous constraint
        partition_constraints_list = [constraint for i, constraint in enumerate(partition_constraints_list) if i == 0 or min(constraint) != min(partition_constraints_list[i-1])]

        if limit_num_constraints is not None:
            partition_constraints_list = partition_constraints_list[:limit_num_constraints]

        utils.logger.debug(f"{__name__}: partition_constraints = {partition_constraints_list}")

        return partition_constraints_list
    
    def _get_lowerbound_with_min_gen_set(self) -> int:

        min_gen_set_start_time = time.perf_counter()
        all_weights = list(set({self.G.edges[e][self.flow_attr] for e in self.G.edges() if self.flow_attr in self.G.edges[e]}))
        # Get the source_flow as the sum of the flow values on all the edges exiting the source nodes
        # (i.e., nodes with in-degree 0)
        source_flow = self._get_source_flow()
        # source_flow = sum(self.G.nodes[n].get(self.flow_attr, 0) for n in self.G.nodes() if self.G.in_degree(n) == 0)
        current_lowerbound_k = self._lowerbound_k if self._lowerbound_k is not None else 1
        min_gen_set_lowerbound = None

        partition_constraints = None
        if self.optimization_options.get("use_min_gen_set_lowerbound_partition_constraints", MinFlowDecompSB.use_min_gen_set_lowerbound_partition_constraints):
            partition_constraints = self._get_partition_constraints_for_min_gen_set(
                min_constraint_len = self.optimization_options.get("use_min_gen_set_lowerbound_partition_constraints_min_constraint_len", MinFlowDecompSB.use_min_gen_set_lowerbound_partition_constraints_min_constraint_len),
                limit_num_constraints = self.optimization_options.get("use_min_gen_set_lowerbound_partition_constraints_limit_num_constraints", MinFlowDecompSB.use_min_gen_set_lowerbound_partition_constraints_limit_num_constraints),
                )
        mingenset_solver_options = copy.deepcopy(self.solver_options)
        if "time_limit" in mingenset_solver_options:
            mingenset_solver_options["time_limit"] = self.time_limit - self.solve_time_elapsed

        mingenset_model = mgs.MinGenSet(
            numbers = all_weights, 
            total = source_flow, 
            weight_type = self.weight_type,
            lowerbound = current_lowerbound_k,
            partition_constraints=partition_constraints,
            remove_sums_of_two = self.optimization_options.get("min_gen_set_remove_sums_of_two", MinFlowDecompSB.min_gen_set_remove_sums_of_two),
            solver_options = mingenset_solver_options,
            )
        mingenset_model.solve()
    
        # If we solved the min gen set problem, we store it and the model, and return the number of elements in the generating set
        if mingenset_model.is_solved():        
            self._generating_set = mingenset_model.get_solution()
            min_gen_set_lowerbound = len(self._generating_set)
            utils.logger.info(f"{__name__}: found a min gen set solution with {min_gen_set_lowerbound} elements ({self._generating_set})")
        else:
            utils.logger.info(f"{__name__}: did NOT find a min gen set solution")
            exit(0)
        
        self.solve_statistics["min_gen_set_solve_time"] = time.perf_counter() - min_gen_set_start_time
        
        return min_gen_set_lowerbound
    
    def _get_lowerbound_with_subgraph_scanning(self) -> int:

        start_time = time.perf_counter()
        topo_order = list(nx.topological_sort(self.G))
        right_node_index = 0
        subgraph_subpath_constraints = []    
        current_lowerbound_k = self._lowerbound_k if self._lowerbound_k is not None else 1    
        subgraph_scanning_lowerbound = 0
        all_subgraph_weights = set()

        right_node_index = MinFlowDecompSB.subgraph_lowerbound_size

        while right_node_index < self.G.number_of_nodes() - 1:
            
            # print("right_node_index", right_node_index)
            subgraph = gu.get_subgraph_between_topological_nodes(
                self.G, 
                topo_order=topo_order, 
                left=right_node_index - MinFlowDecompSB.subgraph_lowerbound_size, 
                right=right_node_index)

            subgraph_subpath_constraints = [c for c in self.subpath_constraints if all(n in subgraph.nodes() for n in c)]
            subgraph_edges_to_ignore = [e for e in self.edges_to_ignore if all(n in subgraph.nodes() for n in e)]
            
            subgraph_optimization_options = copy.deepcopy(self.optimization_options)
            subgraph_optimization_options["use_subgraph_scanning_lowerbound"] = False
            subgraph_optimization_options["lowerbound_k"] = current_lowerbound_k

            subgraph_solver_options = copy.deepcopy(self.solver_options)
            if "time_limit" in subgraph_solver_options:
                subgraph_solver_options["time_limit"] = self.time_limit - self.solve_time_elapsed

            subgraph_mfd_solver = MinFlowDecomp(
                    G=subgraph,
                    flow_attr=self.flow_attr,
                    weight_type=self.weight_type,
                    subpath_constraints=subgraph_subpath_constraints,
                    subpath_constraints_coverage=self.subpath_constraints_coverage,
                    subpath_constraints_coverage_length=self.subpath_constraints_coverage_length,
                    length_attr=self.length_attr,
                    elements_to_ignore=subgraph_edges_to_ignore,
                    optimization_options=subgraph_optimization_options,
                    solver_options=subgraph_solver_options,
                )
            subgraph_mfd_solver.solve()
            if subgraph_mfd_solver.is_solved():
                subgraph_mfd_solution = subgraph_mfd_solver.get_solution()
                subgraph_scanning_lowerbound = max(subgraph_scanning_lowerbound, len(subgraph_mfd_solution["weights"]))
                all_subgraph_weights.update(subgraph_mfd_solution["weights"])

            right_node_index = min(right_node_index + MinFlowDecompSB.subgraph_lowerbound_shift, self.G.number_of_nodes() - 1)

        # Removing zero weight
        if 0 in all_subgraph_weights:
            all_subgraph_weights.remove(0)

        self._all_subgraph_weights = all_subgraph_weights

        self.solve_statistics["subgraph_scanning_time"] = time.perf_counter() - start_time
        
        return subgraph_scanning_lowerbound if subgraph_scanning_lowerbound > 0 else None

    @property
    def solve_time_elapsed(self):
        """
        Returns the elapsed time since the start of the solve process.

        Returns
        -------
        - `float`
        
            The elapsed time in seconds.
        """
        return time.perf_counter() - self.solve_time_start if self.solve_time_start is not None else None

    def get_solution(self):
        """
        Retrieves the solution for the flow decomposition problem.

        Returns
        -------
        - `solution: dict`
        
            A dictionary containing the solution paths (key `"paths"`) and their corresponding weights (key `"weights"`).

        Raises
        -------
        - `exception` If model is not solved.
        """
        self.check_is_solved()
        return self._solution
    
    def get_objective_value(self):

        self.check_is_solved()

        # Number of paths
        return len(self._solution["paths"])

    def is_valid_solution(self) -> bool:
        return self.fd_model.is_valid_solution()
    
    def get_lowerbound_k(self):

        if self._lowerbound_k != None:
            return self._lowerbound_k
        
        stG = stdigraph.stDiGraph(self.G)

        self._lowerbound_k = self.optimization_options.get("lowerbound_k", 1)

        all_weights = set({int(self.G.edges[e][self.flow_attr]) for e in self.G.edges() if self.flow_attr in self.G.edges[e]})
        
        self._lowerbound_k = max(self._lowerbound_k, math.ceil(math.log2(len(all_weights))))

        self._lowerbound_k = max(self._lowerbound_k, stG.get_width(edges_to_ignore=self.edges_to_ignore))

        if self.optimization_options.get("use_min_gen_set_lowerbound", MinFlowDecompSB.use_min_gen_set_lowerbound):  
            mingenset_lowerbound = self._get_lowerbound_with_min_gen_set()
            if mingenset_lowerbound is not None:
                self._lowerbound_k = max(self._lowerbound_k, mingenset_lowerbound)

        if self.optimization_options.get("use_subgraph_scanning_lowerbound", MinFlowDecompSB.use_subgraph_scanning_lowerbound):
            subgraph_scanning_lowerbound = self._get_lowerbound_with_subgraph_scanning()
            if subgraph_scanning_lowerbound is not None:
                self._lowerbound_k = max(self._lowerbound_k, subgraph_scanning_lowerbound)
        
        return self._lowerbound_k
