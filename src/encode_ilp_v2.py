import math

from ortools.linear_solver import pywraplp
from collections import defaultdict
from src.intron_graph import VERTEX_polya, VERTEX_polyt, VERTEX_read_end, VERTEX_read_start

def flatten(l):       return [item for sublist in l for item in sublist]
def head(x):          return int(x.name().split('_')[1])
def tail(x):          return int(x.name().split('_')[2])
def exists_edge(e,p): return len(list(filter(lambda x: head(e)==x[0] and tail(e)==x[1], p))) >= 1

class Enc:

    def __init__(self,n,E,F={},F_low={},F_high={},R=[],solver_id='SCIP'):
        self.n = n
        self.m = len(E)
        self.source = 0
        self.target = self.n-1
        self.E = E
        self.F = F
        if not F:
            self.F = F_high
        self.F_low  = F_low
        self.F_high = F_high
        self.R = R
        self.f = 0
        self.k = 0
        self.w_max = 0
        self.constraints = []

        self.solver = pywraplp.Solver.CreateSolver(solver_id)
        if not self.solver:
            exit(0) #TODO: display error message and catch exception
        
        for e in self.E:
            u,_=e
            if u==self.source:
                self.f += self.F[e] #sum of out-going flow from the source
                self.k += 1         #k = |{v in V : f(s,v)>0}| trivial lower bound, think of funnels. (we assume that every edge has positive flow)
            self.w_max = max(self.w_max,self.F[e])
        
        self.edge_vars = []
        self.pi_vars   = []
        self.weights   = []
        self.path_vars = []

        assert(self.F == self.F_high)
        print(self.F_low)
        print(self.F_high)

    def add_constraint(self, constraint):
        self.constraints.append(str(constraint))

    def show_constraints(self):
        #TODO
        return self.constraints

    def clear(self):
        self.constraints = []
        self.edge_vars   = []
        self.pi_vars     = []
        self.weights     = []
        self.path_vars   = []
        self.solver.Clear()

    def display_stats(self):
        return

    def encode(self):
        for i in range(self.k): #DANGER the order of the variables in edge_varibles and self.pi_vars must be the same
            self.edge_vars += [ [self.solver.BoolVar('x_{}_{}_{}'.format(e[0],e[1],i))              for e in self.E ] ]
            self.pi_vars   += [ [self.solver.IntVar(0, self.w_max,'p_{}_{}_{}'.format(e[0],e[1],i)) for e in self.E ] ]
            self.weights   += [  self.solver.IntVar(1, self.w_max,'w_{}'.format(i)) ]
            self.path_vars += [ [self.solver.BoolVar('r_{}_{}'.format(i,j)) for j in range(len(self.R)) ]]

        print(self.edge_vars)
        print(self.pi_vars)
        print(self.weights)
        print(self.path_vars)

        #3a, 3b
        for i in range(self.k):
            self.solver.Add( sum( filter(lambda edge: head(edge)==self.source, self.edge_vars[i]) ) == 1 )
            self.solver.Add( sum( filter(lambda edge: tail(edge)==self.target, self.edge_vars[i]) ) == 1 )
            
        #3c
        for i in range(self.k):
            for v in range(1,self.n-1): #find all wedges u->v->w for a fixed v (excluding {s,t})
                self.solver.Add( sum( filter(lambda edge: tail(edge)==v, self.edge_vars[i]) ) - sum( filter(lambda edge: head(edge)==v, self.edge_vars[i]) ) == 0 )

        def EncodeExactFlow():
            #5a (to be exchanged with constraint 9a)
            for e in range(self.m):
                self.solver.Add( sum( self.pi_vars[i][e] for i in range(self.k)) == self.F[self.E[e]] )

        #5b
        for i in range(self.k):
            for pi,x in tuple(zip(self.pi_vars[i],self.edge_vars[i])):
                self.solver.Add( pi <= self.w_max * x )
        #5c
        for i in range(self.k):
            for pi in self.pi_vars[i]:
                self.solver.Add( pi <= self.weights[i] )
        #5d
        for i in range(self.k):
            for pi,x in tuple(zip(self.pi_vars[i],self.edge_vars[i])):
                self.solver.Add( pi >= self.weights[i] - (1 - x) * self.w_max )
        
        #Example of a subpath constraint: R=[ [(1,3),(3,5)], [(0,1)] ], means that we have 2 paths to cover, the first one is 1-3-5. the second path is just a single edge 0-1
        def EncodeSubpathConstraints():
            #7a
            for i in range(self.k):
                for j in range(len(self.R)):
                    self.solver.Add( sum( filter( lambda e: exists_edge(e,self.R[j]), self.edge_vars[i] )) >= len(self.R[j])*self.path_vars[i][j] )
            #7b
            for r in range(len(self.R)):
                self.solver.Add( sum( self.path_vars[i][r] for i in range(len(self.path_vars)) ) >= 1 )
            
        def EncodeInexactFlow():
            #9a (to be exchanged with constraint 5a)
            for e in range(self.m):
                self.solver.Add( sum( self.pi_vars[i][e] for i in range(self.k) ) <= self.F_high[self.E[e]] )
                self.solver.Add( sum( self.pi_vars[i][e] for i in range(self.k) ) >= self.F_low [self.E[e]] )

        def EncodePathError():
            #9a (to be exchanged with constraint 5a)
            for e in range(self.m):
                self.solver.Add( sum( self.pi_vars[i][e] for i in range(self.k) ) <= self.F_high[self.E[e]] )
                self.solver.Add( sum( self.pi_vars[i][e] for i in range(self.k) ) >= self.F_low [self.E[e]] )

        if self.R!=[]:
            EncodeSubpathConstraints()
        if len(self.F_low) != 0 and len(self.F_high) != 0:
            EncodeInexactFlow()
        else:
            EncodeExactFlow()

        #Add a trivial objective function - we are just interested in deciding whether or not the constraints form a feasible region
        self.solver.Minimize(1)

    def print_solution(self):
        print("Solution has size ", self.k, " with the following weight-path decomposition")
        for i in range(self.k):
            path = list(map(lambda e : tail(e), list(filter(lambda e : e.solution_value()==1, self.edge_vars[i]))))
            print(int(self.weights[i].solution_value()), path[:-1])

    def linear_search(self):
        print(self.m,self.f)
        while self.k <= min(self.m,self.f):

            self.encode()
            if self.solver.Solve() == pywraplp.Solver.OPTIMAL:
                self.print_solution()
                break

            self.clear()
            self.k += 1


'''
# Transform intron_graph into a flow matrix F
- add super source S and an edge (S,v) to every node v with 0 in-degree (leftmost guys in each layer); add super target T and an edge (v,T) from every node v with 0-outdegree (rightmost guys in each layer)
- define f(S,v)=out-going flow of v and f(v,T)=incoming-flow of v
- DAG may contain multiple connected components but thats fine, we can still add from super source to every 0-indegree guy... in the future think of possible opts (threads running in different components?)
'''

def intron_to_matrix(intron_graph):
    edge_count, edge_list, flow_dict = intron_graph.to_adjacency_list()

    alfa = 0.5
    f_low   = {}
    f_high = {}
    for edge in flow_dict:
        f_low[edge] = math.floor((1-alfa)*flow_dict[edge])
        f_high[edge] = math.ceil((1+alfa)*flow_dict[edge])

    return edge_count, edge_list, f_low, f_high


def Encode_ILP(intron_graph):

    n,E,Fl,Fh = intron_to_matrix(intron_graph)

    R = [[(1,3),(3,5)]]

    e = Enc(n,E,F_low=Fl,F_high=Fh)
    e.linear_search()


#https://developers.google.com/optimization/reference/python/linear_solver/pywraplp