import math

from ortools.linear_solver import pywraplp
from collections import defaultdict
from src.intron_graph import VERTEX_polya, VERTEX_polyt, VERTEX_read_end, VERTEX_read_start
from src.visualization import visualize

def flatten(l):       return [item for sublist in l for item in sublist]
def head(x):          return int(x.name().split('_')[1])
def tail(x):          return int(x.name().split('_')[2])
def exists_edge(e,p): return len(list(filter(lambda x: head(e)==x[0] and tail(e)==x[1], p))) >= 1
def h(e):             return e[0]
def t(e):             return e[1]

class Enc:

    def __init__(self,n,E,F={},R=[],solver_id='GLOP_LINEAR_PROGRAMMING',eps=1000):
        self.n = n
        self.m = len(E)
        self.source = 0
        self.target = self.n-1
        self.E = E
        self.F = F
        self.R = R
        self.f = 0
        self.k = 0
        self.w_max = 0
        self.constraints = []

        self.epsilon = eps #we could have a function that computes epsilon based on features of the data

        self.solver = pywraplp.Solver.CreateSolver(solver_id)
        if not self.solver:
            exit(0) #TODO: display error message and catch exception
        
        for e in self.E:
            u,_=e
            if u==self.source:
                self.f += self.F[e] #sum of out-going flow from the source
                self.k += 1         #k = |{v in V : f(s,v)>0}| trivial lower bound, think of funnels. (we assume that every edge has positive flow)
            self.w_max = max(self.w_max,self.F[e])
        
        self.edge_vars  = []
        self.phi_vars   = []
        self.gam_vars   = []
        self.weights    = []
        self.slack_vars = []
        self.path_vars  = []

    def add_constraint(self, constraint):
        self.constraints.append(str(constraint))

    def show_constraints(self):
        #TODO
        return self.constraints

    def clear(self):
        self.constraints = []
        self.edge_vars   = []
        self.phi_vars    = []
        self.gam_vars    = []
        self.slack_vars  = []
        self.weights     = []
        self.path_vars   = []
        self.solver.Clear()

    def display_stats(self):
        return

    def encode(self):

        for i in range(self.k): #DANGER the order of the variables, e.g. order of edge_varibles and self.phi_vars must be the same
            self.edge_vars  += [ [self.solver.BoolVar('x_{}_{}_{}'.format(h(e),t(e),i))              for e in self.E ] ]
            self.phi_vars   += [ [self.solver.IntVar(0, self.w_max,'f_{}_{}_{}'.format(h(e),t(e),i)) for e in self.E ] ]
            self.gam_vars   += [ [self.solver.IntVar(0, self.w_max,'g_{}_{}_{}'.format(h(e),t(e),i)) for e in self.E ] ]
            self.weights    += [  self.solver.IntVar(1, self.w_max,'w_{}'.format(i)) ]
            self.slack_vars += [  self.solver.IntVar(0, self.w_max,'s_{}'.format(i)) ]
            self.path_vars  += [ [self.solver.BoolVar('r_{}_{}'.format(i,j)) for j in range(len(self.R)) ]]

        print()
        print(self.edge_vars)
        print(self.phi_vars)
        print(self.gam_vars)
        print(self.weights)
        print(self.slack_vars)
        print(self.path_vars)
        print()

        #The identifiers of the constraints come from https://www.biorxiv.org/content/10.1101/2023.03.20.533019v1.full.pdf

        #14a, 14b
        for i in range(self.k):
            self.solver.Add( sum( filter(lambda edge: head(edge)==self.source, self.edge_vars[i]) ) == 1 )
            self.solver.Add( sum( filter(lambda edge: tail(edge)==self.target, self.edge_vars[i]) ) == 1 )
            
        #14c
        for i in range(self.k):
            for v in range(1,self.n-1): #find all wedges u->v->w for a fixed v (excluding {s,t})
                self.solver.Add( sum( filter(lambda edge: tail(edge)==v, self.edge_vars[i]) ) - sum( filter(lambda edge: head(edge)==v, self.edge_vars[i]) ) == 0 )

        #14d, 14e
        for e in range(self.m):
            self.solver.Add( self.F[self.E[e]] - sum(self.phi_vars[i][e] for i in range(self.k)) <=  sum(self.gam_vars[i][e] for i in range(self.k)) )
            self.solver.Add( self.F[self.E[e]] - sum(self.phi_vars[i][e] for i in range(self.k)) >= -sum(self.gam_vars[i][e] for i in range(self.k)) )

        #14f, 14i
        for i in range(self.k):
            for e in range(self.m): # for phi,x in tuple(zip(self.phi_vars[i],self.edge_vars[i])): self.solver.Add( phi <= self.w_max * x )
                self.solver.Add( self.phi_vars[i][e] <= self.w_max * self.edge_vars[i][e] )
                self.solver.Add( self.gam_vars[i][e] <= self.w_max * self.edge_vars[i][e] )

        #14g, 14j
        for i in range(self.k):
            for phi in self.phi_vars[i]:
                self.solver.Add( phi <= self.weights[i] )
            for gam in self.gam_vars[i]:
                self.solver.Add( gam <= self.slack_vars[i] )

        #14h, 14k
        for i in range(self.k):
            for e in range(self.m):
                self.solver.Add( self.phi_vars[i][e] >= self.weights[i]    - (1 - self.edge_vars[i][e]) * self.w_max )
                self.solver.Add( self.gam_vars[i][e] >= self.slack_vars[i] - (1 - self.edge_vars[i][e]) * self.w_max )

        #Example of a subpath constraint: R=[ [(1,3),(3,5)], [(0,1)] ], means that we have 2 paths to cover, the first one is 1-3-5. the second path is just a single edge 0-1
        def EncodeSubpathConstraints():
            #7a
            for i in range(self.k):
                for j in range(len(self.R)):
                    self.solver.Add( sum( filter( lambda e: exists_edge(e,self.R[j]), self.edge_vars[i] )) >= len(self.R[j])*self.path_vars[i][j] )
            #7b
            for r in range(len(self.R)):
                self.solver.Add( sum( self.path_vars[i][r] for i in range(len(self.path_vars)) ) >= 1 )

        if self.R!=[]:
            EncodeSubpathConstraints()

        #Add objective function: the goal is to minimize the total slack of the paths
        self.solver.Minimize( sum( slack for slack in self.slack_vars ) )

    def print_solution(self,solution):
        opt, slack, paths = solution
        print("Solution has size ", opt, "with slack", slack, "and with the following weight-slack-path decomposition")
        for p in paths:
            print(p)
    
    def build_solution(self):
        opt   = self.k
        slack = self.solver.Objective().Value()
        paths = []
        for i in range(self.k):
            path = list(map(lambda e : tail(e), list(filter(lambda e : e.solution_value()==1, self.edge_vars[i]))))
            paths.append( (int(self.weights[i].solution_value()), int(self.slack_vars[i].solution_value()), path[:-1]) )
        return (opt, slack, paths)

    def linear_search(self):
        print("Feasibility")
        #Feasibility: Find first feasible solution (there always exists one) and clear the solver for next stage
        while True:
            self.encode()
            status = self.solver.Solve()
            if (status == pywraplp.Solver.OPTIMAL):
                previous_slack = self.solver.Objective().Value()
                solution       = self.build_solution()
                break
            self.clear()
            self.k += 1

        print("Optimality")
        #Optimality: Find the k for which the difference in slacks of two consecutive iterations becomes sufficiently small
        while self.k <= min(self.m,self.f):

            self.encode()
            status = self.solver.Solve()
            current_slack = self.solver.Objective().Value()
            
            if previous_slack-current_slack < self.epsilon:
                print(previous_slack,current_slack)
                self.print_solution(solution)
                exit(0)
                #break

            solution       = self.build_solution()
            previous_slack = current_slack
            self.clear()
            self.k += 1


'''
# Transform intron_graph into a flow matrix F
- add super source S and an edge (S,v) to every node v with 0 in-degree (leftmost guys in each layer); add super target T and an edge (v,T) from every node v with 0-outdegree (rightmost guys in each layer)
- define f(S,v)=out-going flow of v and f(v,T)=incoming-flow of v
- DAG may contain multiple connected components but thats fine, we can still add from super source to every 0-indegree guy... in the future think of possible opts (threads running in different components?)

Input:  intron_graph
Output: abstracted DAG in the form (number of nodes, edge list, flow dictionary)
''' 
def intron_to_matrix(intron_graph):
    intron2vertex = dict()
    vertex2intron = dict()
    vertex_id = 1

    # add all intron vertices
    for intron in intron_graph.intron_collector.clustered_introns.keys():
        intron2vertex[intron] = vertex_id
        vertex2intron[vertex_id] = intron
        vertex_id += 1

    # add all starting vertices (avoid repeating incoming from the same starting vertex)
    for intron in intron_graph.incoming_edges.keys():
        for preceeding_intron in intron_graph.incoming_edges[intron]:
            if preceeding_intron[0] in [VERTEX_polyt, VERTEX_read_start] and preceeding_intron not in intron2vertex:
                intron2vertex[preceeding_intron] = vertex_id
                vertex2intron[vertex_id] = preceeding_intron
                vertex_id += 1

    # add all terminal vertices (avoid repeating outgoing from the same terminal vertex)
    for intron in intron_graph.outgoing_edges.keys():
        for subsequent_intron in intron_graph.outgoing_edges[intron]:
            if subsequent_intron[0] in [VERTEX_polya, VERTEX_read_end] and subsequent_intron not in intron2vertex:
                intron2vertex[subsequent_intron] = vertex_id
                vertex2intron[vertex_id] = subsequent_intron
                vertex_id += 1

    source = 0
    target = vertex_id

    # create edges
    edge_list = []
    edge_set = set()
    starting_introns = defaultdict(int)
    for intron in intron_graph.incoming_edges.keys():
        for preceeding_intron in intron_graph.incoming_edges[intron]:
            edge_list.append((intron2vertex[preceeding_intron], intron2vertex[intron]))
            edge_set.add((intron2vertex[preceeding_intron], intron2vertex[intron]))
            if preceeding_intron[0] in [VERTEX_polyt, VERTEX_read_start]:
                starting_introns[preceeding_intron] += intron_graph.edge_weights[(preceeding_intron, intron)]

    terminal_introns = defaultdict(int)
    for intron in intron_graph.outgoing_edges.keys():
        for subsequent_intron in intron_graph.outgoing_edges[intron]:
            if subsequent_intron[0] in [VERTEX_polya, VERTEX_read_end]:
                edge_list.append((intron2vertex[intron], intron2vertex[subsequent_intron]))
                edge_set.add((intron2vertex[intron], intron2vertex[subsequent_intron]))
                terminal_introns[subsequent_intron] += intron_graph.edge_weights[(intron, subsequent_intron)]

    flow_dict = defaultdict(int)
    for intron in intron_graph.incoming_edges.keys():
        for preceeding_intron in intron_graph.incoming_edges[intron]:
            u = intron2vertex[preceeding_intron]
            v = intron2vertex[intron]
            flow_dict[(u, v)] = intron_graph.edge_weights[(preceeding_intron, intron)]

    for intron in intron_graph.outgoing_edges.keys():
        for subsequent_intron in intron_graph.outgoing_edges[intron]:
            if subsequent_intron[0] in [VERTEX_polya, VERTEX_read_end]:
                u = intron2vertex[intron]
                v = intron2vertex[subsequent_intron]
                flow_dict[(u, v)] = intron_graph.edge_weights[(intron, subsequent_intron)]

    # add connection to super source and total weight
    for starting_intron in starting_introns.keys():
        starting_vertex = intron2vertex[starting_intron]
        edge_list.append((source, starting_vertex))
        edge_set.add((source, starting_vertex))
        flow_dict[(source, starting_vertex)] = starting_introns[starting_intron]

    for terminal_intron in terminal_introns.keys():
        terminal_vertex = intron2vertex[terminal_intron]
        edge_list.append((terminal_vertex, target))
        edge_set.add((terminal_vertex, target))
        flow_dict[(terminal_vertex, target)] = terminal_introns[terminal_intron]

    #FIXME edge_set is not necessary
    assert len(edge_list) == len(edge_set)
    assert len(edge_list) == len(flow_dict)

    print("Returning from graph transformation")

    print("N:",vertex_id+1)
    print("Edges:",edge_list)
    print("Flows:",flow_dict)

    return vertex_id+1, edge_list, flow_dict


def Encode_ILP(intron_graph):

    n,E,F = intron_to_matrix(intron_graph)

    visualize((E,F))

    e = Enc(n,E,F)
    e.linear_search()


#https://developers.google.com/optimization/reference/python/linear_solver/pywraplp