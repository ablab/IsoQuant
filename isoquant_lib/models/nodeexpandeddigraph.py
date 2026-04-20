import networkx as nx
from flowpaths.utils import graphutils as gu
import flowpaths.utils as utils
from copy import deepcopy

class NodeExpandedDiGraph(nx.DiGraph):
    
    def __init__(
            self,
            G: nx.DiGraph,
            node_flow_attr: str,
            node_cell_flow_attr: str = None,
            try_filling_in_missing_flow_attr: bool = False,
            node_length_attr: str = None,
            additional_starts: list = [],
            additional_ends: list = [],
            ):
        """
        This class is a subclass of the networkx DiGraph class. It is used to represent a directed graph
        where all nodes `v` have been "expanded" or "subdivided" into an edge `(v.0, v.1)`. This is useful for representing
        graphs where the flow values, or weights, are associated with the nodes, rather than the edges. 
        These expanded edges are then added to the `edges_to_ignore` list, available as a property of this class.

        !!! info "Using this class"

            - Create a `NodeExpandedDiGraph` object by passing a directed graph `G` and the attribute name `node_flow_attr` from where to get the flow values / weights on the nodes.
            - Pass the edges from the `edges_to_ignore` attribute of this class to the decomposition models, in order to ignore all original edges of the graph,
              and thus consider in the constraints only the new edges added in the expanded graph (which have flow values).
            - Solve the decomposition model on the expanded graph.
            - Use the `get_condensed_paths` method to condense the solution paths (which are in the expanded graph) to paths in the original graph.

        Parameters
        ----------
        - `G : nx.DiGraph`
            
            The input directed graph, as networkx DiGraph.

        - `node_flow_attr : str`

            The attribute name from where to get the flow values / weights on the nodes. 
            This attribute for each `v` is then set to the edge `(v.0, v.1)` connecting the new expanded nodes.
            This attribute can be missing from some nodes of the graph, in which case 
            the edge `(v.0, v.1)` is added to `edges_to_ignore`.

        - `try_filling_in_missing_flow_attr : bool`, optional

            If `True`, try filling in missing flow values in the expanded graph (i.e., if some original edge does not have the attribute specified by `node_flow_attr`), 
            by setting them to the flow values of a flow between the sources and the sinks of the original graph, such that for every edge where `node_flow_attr` is present,
            the demand and capacity are equal to the value specified by `node_flow_attr`.
            Default is `False`.

        - `node_length_attr : str`, optional

            The attribute name from where to get the length values on the nodes. Default is `None`. 
            If you specify this attribute, it may be missing from some nodes of the graph.

        !!! example "Example"

            ```python
            import flowpaths as fp
            import networkx as nx

            graph = nx.DiGraph()
            graph.add_node("s", flow=13)
            graph.add_node("a", flow=6)
            graph.add_node("b", flow=9)
            graph.add_node("c", flow=13)
            graph.add_node("d", flow=6)
            graph.add_node("t", flow=13)

            # Adding edges
            graph.add_edges_from([("s", "a"), ("s", "b"), ("a", "b"), ("a", "c"), ("b", "c"), ("c", "d"), ("c", "t"), ("d", "t")])

            # Expand the graph
            ne_graph = fp.NodeExpandedDiGraph(graph, node_flow_attr="flow")

            # Solve the problem on the expanded graph
            mfd_model = fp.MinFlowDecomp(
                ne_graph, 
                flow_attr="flow",
                edges_to_ignore=ne_graph.edges_to_ignore,
                )
            mfd_model.solve()

            if mfd_model.is_solved():
                # Getting the solution in the expanded graph
                solution = mfd_model.get_solution()
                # Condensing the paths in the expanded graph to paths in the the original graph
                original_paths = ne_graph.get_condensed_paths(solution["paths"])
                print("Original paths:", original_paths)
                print("Weights:", solution["weights"])
            ```
        """
        super().__init__()

        if G.number_of_nodes() == 0:
                utils.logger.error(f"{__name__}: The input graph G has no nodes. Please provide a graph with at least one node.")
                raise ValueError(f"The input graph G has no nodes. Please provide a graph with at least one node.")
        if not all(isinstance(node, str) for node in G.nodes()):
            utils.logger.error(f"{__name__}: Graph id {utils.fpid(G)}: every node of the graph must be a string.")
            raise ValueError("Every node of the graph must be a string.")

        self.original_G = deepcopy(G)

        if "id" in G.graph:
            self.graph["id"] = G.graph["id"]
            
        self.node_flow_attr = node_flow_attr
        self.node_cell_flow_attr = node_cell_flow_attr
        self.node_length_attr = node_length_attr

        self.global_source_id = 'source' + str(id(self))
        self.global_sink_id = 'sink' + str(id(self))

        self._edges_to_ignore = []

        for node in G.nodes:
            node0 = node + '.0'
            node1 = node + '.1'
            self.add_node(node0, **G.nodes[node])
            self.add_node(node1, **G.nodes[node])
            self.add_edge(node0, node1, **G.nodes[node])
            if self.node_flow_attr in G.nodes[node]:
                self[node0][node1][self.node_flow_attr] = G.nodes[node][self.node_flow_attr]
            else:
                self._edges_to_ignore.append((node0, node1))
            if self.node_length_attr is not None:
                if self.node_length_attr in G.nodes[node]:
                    self[node0][node1][self.node_length_attr] = G.nodes[node][self.node_length_attr]

            # Adding in-coming edges
            for pred in G.predecessors(node):
                pred1 = pred + '.1'
                self.add_edge(pred1, node0, **G.edges[pred, node])
                self._edges_to_ignore.append((pred1, node0))
                
                # If the edge (pred,node) does not have the length attribute, set it to 0
                if self.node_length_attr is not None:
                    if self.node_length_attr not in G.edges[pred, node]:
                        self[pred1][node0][self.node_length_attr] = 0

            # Adding out-going edges
            for succ in G.successors(node):
                succ0 = succ + '.0'
                self.add_edge(node1, succ0, **G.edges[node, succ])
                # This is not necessary, as the edge (node1, succ0) has already been added above, for succ
                # self._edges_to_ignore.append((node1, succ0))

        if (len(additional_starts) + len(additional_ends) > 0) and not try_filling_in_missing_flow_attr:
            utils.logger.error(f"{__name__}: If `additional_starts` or `additional_ends` are specified, `try_filling_in_missing_flow_attr` must be set to True.")
            raise ValueError("If `additional_starts` or `additional_ends` are specified, `try_filling_in_missing_flow_attr` must be set to True.")

        if additional_starts != []:
            self.add_node(self.global_source_id + '.0')
            self.add_node(self.global_source_id + '.1')
            new_edge = (self.global_source_id + '.0', self.global_source_id + '.1')
            self.add_edges_from([new_edge])
            self._edges_to_ignore.append(new_edge)

            for node in additional_starts:
                if node not in G.nodes:
                    utils.logger.error(f"{__name__}: Node {node} not in the original graph.")
                    raise ValueError(f"Node {node} not in the original graph.")
                new_edge = (self.global_source_id + '.1', node + '.0', )
                self.add_edges_from([new_edge])
                self._edges_to_ignore.append(new_edge)

        if additional_ends != []:
            self.add_node(self.global_sink_id + '.0')
            self.add_node(self.global_sink_id + '.1')
            new_edge = (self.global_sink_id + '.0', self.global_sink_id + '.1')
            self.add_edges_from([new_edge])
            self._edges_to_ignore.append(new_edge)
            for node in additional_ends or []:
                if node not in G.nodes:
                    utils.logger.error(f"{__name__}: Node {node} not in the original graph.")
                    raise ValueError(f"Node {node} not in the original graph.")
                new_edge = (node + '.1', self.global_sink_id + '.0')
                self.add_edges_from([new_edge])
                self._edges_to_ignore.append(new_edge)

        if try_filling_in_missing_flow_attr:
            self._try_filling_in_missing_flow_values()

        nx.freeze(self)    

    def _try_filling_in_missing_flow_values(self):
        
        # Fills in missing flow values in the expanded graph, by setting them to the flow values of a maximum flow
        # between the source and the sink of the original graph, with capacity equal to the flow values of the nodes.
        
        # Create a max flow instance where the edges with missing flow have infinite capacity
        # and the edges with flow have capacity equal to the flow value
        # The source are the sources of the original graph, and the sink are the sinks of the original graph
        # We solve this with the networkx package
        
        min_cost_flow_network = nx.DiGraph(self)

        source_node = f'source{id(min_cost_flow_network)}'
        sink_node = f'sink{id(min_cost_flow_network)}'
        demands_attr = f'demand{id(min_cost_flow_network)}'
        capacities_attr = f'capacities{id(min_cost_flow_network)}'
        costs_attr = f'costs{id(min_cost_flow_network)}'
        
        min_cost_flow_network.add_node(source_node)
        min_cost_flow_network.add_node(sink_node)

        # For every node in the graph, add an edge from the source to the node.0 with infinite capacity
        for node in self.original_G.nodes:
            if self.original_G.in_degree(node) == 0:
                min_cost_flow_network.add_edge(source_node, node + '.0')
                min_cost_flow_network[source_node][node + '.0'].update({demands_attr: 0, capacities_attr: float('inf'), costs_attr: 0})
            if self.original_G.out_degree(node) == 0:
                min_cost_flow_network.add_edge(node + '.1', sink_node)
                min_cost_flow_network[node + '.1'][sink_node].update({demands_attr: 0, capacities_attr: float('inf'), costs_attr: 0})

        for u, v, data in self.edges(data=True):
            if self.node_flow_attr not in data:
                min_cost_flow_network[u][v].update({demands_attr: 0, capacities_attr: float('inf'), costs_attr: 0})
            else:
                min_cost_flow_network[u][v].update({demands_attr: data[self.node_flow_attr], capacities_attr: data[self.node_flow_attr], costs_attr: 0})
        
        flow_value, flow_dict = gu.min_cost_flow(
            G = min_cost_flow_network, 
            s = source_node, 
            t = sink_node, 
            demands_attr = demands_attr,
            capacities_attr = capacities_attr,
            costs_attr = costs_attr,
            )

        if flow_dict is not None:
            for edge in self.edges:
                if self.node_flow_attr not in self.edges[edge]:
                    self.edges[edge][self.node_flow_attr] = flow_dict[edge[0]][edge[1]]

    @property
    def edges_to_ignore(self):
        """
        List of edges to ignore when solving the decomposition model on the expanded graph. 

        These are the edges of the original graph, since only the new edges that have been introduced 
        for every node must considered in the decomposition model, with flow value from the node attribute `node_flow_attr`.
        """
        return self._edges_to_ignore
    
    def get_expanded_additional_starts(self, additional_starts):
        
        return [self.get_expanded_edge(node)[0] for node in additional_starts]
    
    def get_expanded_additional_ends(self, additional_ends):
        
        return [self.get_expanded_edge(node)[1] for node in additional_ends]
    
    def get_expanded_subpath_constraints(self, subpath_constraints):
        """
        Expand a list of subpath constraints from the original graph (where every constraint is a list **nodes** or **edges**
        in the original graph) to a list of subpath constraints in the expanded graph (where every constraint 
        is a list of edges in the expanded graph). 
        
        - If the constraints are lists of nodes, then the corresponding edges are of type `('node.0', 'node.1')` 
        - If the constraints are lists of edges of the form `(u,v)`, then the corresponding edges are of type `('u.1', 'v.0')`

        Parameters
        ----------
        - `subpath_constraints : list`
            
            List of subpath constraints in the original graph.

        Returns
        -------
        - `expanded_constraints : list`
            
            List of subpath constraints in the expanded graph.
        """
        
        # Check if the subpath constraints are lists of node or lists of edges,
        # and expand them accordingly, using the two functions already implemented.
        
        if not isinstance(subpath_constraints, list):
            utils.logger.error(f"{__name__}: Subpath constraints must be a list.")
            raise ValueError("Subpath constraints must be a list.")
        if not all(isinstance(constraint, list) for constraint in subpath_constraints):
            utils.logger.error(f"{__name__}: Subpath constraints must be a list of lists.")
            raise ValueError("Subpath constraints must be a list of lists.")
        
        if len(subpath_constraints) == 0:
            return []
        
        if isinstance(subpath_constraints[0][0], str):
            return self._get_expanded_subpath_constraints_nodes(subpath_constraints)
        elif isinstance(subpath_constraints[0][0], tuple):
            return self._get_expanded_subpath_constraints_edges(subpath_constraints)
        else:
            utils.logger.error(f"{__name__}: Subpath constraints must be a list of lists of nodes or edges.")
            raise ValueError("Subpath constraints must be a list of lists of nodes or edges.")

    def _get_expanded_subpath_constraints_nodes(self, subpath_constraints):
        """
        Expand a list of subpath constraints from the original graph (where every constraint is a list **nodes**
        in the original graph) to a list of subpath constraints in the expanded graph (where every constraint 
        is a list of edges (of type `('node.0', 'node.1`)` in the expanded graph).

        Parameters
        ----------
        - `subpath_constraints : list`
            
            List of subpath constraints in the original graph (as lists of **nodes**).

        Returns
        -------
        - `expanded_constraints : list`
            
            List of subpath constraints in the expanded graph (as lists of edges of type `('node.0', 'node.1`)`).
        """

        expanded_constraints = []

        for constraint in subpath_constraints:
            expanded_constraint = []
            for node in constraint:
                if node not in self.original_G.nodes:
                    utils.logger.error(f"{__name__}: Node {node} not in the original graph.")
                    raise ValueError(f"Node {node} not in the original graph.")
                expanded_constraint.append((node + '.0', node + '.1'))
            expanded_constraints.append(expanded_constraint)

        return expanded_constraints
    
    def _get_expanded_subpath_constraints_edges(self, subpath_constraints):
        """
        Expand a list of subpath constraints from the original graph (where every constraint is a list **edges** 
        in the original graph) to a list of subpath constraints in the expanded graph where every constraint 
        is a list of the corresponding edges in the expanded graph). 

        The subpath constraints also contain the expanded nodes, so for every edge `(u,v)` in the subpath constraint,
        we add the expanded node `u` (as `(u.0, u.1)`) and the edge `(u.1, v.0)` in the expanded graph.
        For the last edge `(u,v)` in the constraint, we also add the expanded node `v` (as `(v.0, v.1)`).

        Parameters
        ----------
        - `subpath_constraints : list`
            
            List of subpath constraints in the original graph (as lists of **edges**).

        Returns
        -------
        - `expanded_constraints : list`
            
            List of subpath constraints in the expanded graph (as lists of corresponding edges).
        """

        expanded_constraints = []

        for constraint in subpath_constraints:
            expanded_constraint = []
            for i, edge in enumerate(constraint):
                if edge not in self.original_G.edges:
                    utils.logger.error(f"{__name__}: Edge {edge} not in the original graph.")
                    raise ValueError(f"Edge {edge} not in the original graph.")
                # For every edge (u,v) in the subpath constraint, we add the expanded node u (as (u.0, u.1)) and the edge (u.1, v.0)
                expanded_constraint.append(self.get_expanded_edge(edge[0]))
                expanded_constraint.append((edge[0] + '.1', edge[1] + '.0'))
                
                # For the last edge (u,v) in the constraint, we also add the expanded node v (as (v.0, v.1))
                if i == len(constraint) - 1:
                    expanded_constraint.append(self.get_expanded_edge(edge[1]))

            expanded_constraints.append(expanded_constraint)

        return expanded_constraints

    def get_expanded_edge(self, graph_element):
        """
        Given a node or an edge in the original graph, return the corresponding expanded edge in the expanded graph.
        
        - If `graph_element` is some `node`, it returns the edge `('node.0', 'node.1')`.
        - If `graph_element` is some `edge`, it returns the edge `('edge[0].1', 'edge[1].0')`.
        
        This is useful for converting nodes and edges in the original graph to the corresponding edges in the expanded graph.
        """

        if isinstance(graph_element, str):   
            node = graph_element         
            if node not in self.original_G.nodes:
                utils.logger.error(f"{__name__}: Node {node} not in the original graph.")
                raise ValueError(f"Node {node} not in the original graph.")
            return (node + '.0', node + '.1')
        
        elif isinstance(graph_element, tuple):
            edge = graph_element
            if edge not in self.original_G.edges:
                utils.logger.error(f"{__name__}: Edge {edge} not in the original graph.")
                raise ValueError(f"Edge {edge} not in the original graph.")
            return (edge[0] + '.1', edge[1] + '.0')
        else:
            utils.logger.error(f"{__name__}: Graph element must be a node or an edge, not {graph_element}.")
            raise ValueError(f"Graph element must be a node or an edge, not {graph_element}.")

    def get_condensed_paths(self, paths):
        """
        Condense a list of paths from the expanded graph to the original graph. 
        
        This assumes that:

        - The nodes in the expanded graph are named as 'node.0' and 'node.1', where 'node' is the name of the node in the
        original graph. 
        - The paths are lists of nodes in the expanded graph, where the nodes are ordered as 'nodeA.0', 'nodeA.1', 'nodeB.0', 'nodeB.1', etc.
        Meaning that we always have two nodes from the same original node in sequence.

        Parameters
        ----------
        - `paths : list`
            
            List of paths in the expanded graph.

        Returns
        -------
        - `condensed_paths: list`
            
            List of paths in the original graph.

        Raises
        ------
        - `ValueError`
            
            - If the node names in the expanded_path on even positions (starting from 0) do not end with `.0`.
            - If these node names (with the suffix `.0` removed) are not in the original graph.
        """

        condensed_paths = []
        for path in paths:
            condensed_path = []
            for i in range(0, len(path) - 1, 2):
                # Raise an error if the last two symbols of path[i] are not '.0'
                if path[i][-2:] != '.0':
                    utils.logger.error(f"{__name__}: Invalid node name in path: {path[i]}")
                    raise ValueError(f"Invalid node name in path: {path[i]}")
                
                node = path[i][:-2]
                # Raise an error if the node is not in the original graph and is not the global source or sink
                if node not in self.original_G.nodes and node not in [self.global_source_id, self.global_sink_id]:
                    utils.logger.error(f"{__name__}: Node {node} not in the original graph.")
                    raise ValueError(f"Node {node} not in the original graph.")
                
                # Append the node to the condensed path only if it is not the global source or sink
                if node not in [self.global_source_id, self.global_sink_id]:
                    condensed_path.append(node)
            condensed_paths.append(condensed_path)
        return condensed_paths
    
    def get_condensed_graph(self):
        """
        Return a condensed version of the expanded graph, by mapping the node attributes `self.node_flow_attr` and `self.node_length_attr` 
        from expanded edges back to the corresponding nodes of the original graph.

        This is useful when creating an expanded graph by some means, or correcting the flow values of an expanded graph, and then condensing this graph.
        """

        condensed_graph = deepcopy(self.original_G)
        for node in condensed_graph.nodes:
            expanded_edge = self.get_expanded_edge(node)
            if self.node_flow_attr in self.edges[expanded_edge]:
                condensed_graph.nodes[node][self.node_flow_attr] = self.edges[expanded_edge][self.node_flow_attr]
            if self.node_cell_flow_attr is not None:
                if self.node_cell_flow_attr in self.edges[expanded_edge]:
                    condensed_graph.nodes[node][self.node_cell_flow_attr] = self.edges[expanded_edge][self.node_cell_flow_attr]
            if self.node_length_attr is not None:
                if self.node_length_attr in self.edges[expanded_edge]:
                    condensed_graph.nodes[node][self.node_length_attr] = self.edges[expanded_edge][self.node_length_attr]
                
        return condensed_graph
    