from ortools.linear_solver import pywraplp

def flatten(l):       return [item for sublist in l for item in sublist]
def head(x):          return int(x.name().split('_')[1])
def tail(x):          return int(x.name().split('_')[2])
def exists_edge(e,p): return len(list(filter(lambda x: head(e)==x[0] and tail(e)==x[1], p))) >= 1

class Enc:

    def __init__(self,F,R=[],F_low=[],F_high=[],solver_id='SCIP'):
        self.n = len(F)
        self.m = len(list(filter(lambda f_uv : f_uv!=0, flatten(F)))) #every edge has positive flow
        self.source = 0
        self.target = self.n-1
        self.F = F
        self.R = R
        self.f = sum(F[self.source]) #sum of out-going flow from the source
        self.w_max = max(map(max,F))
        self.constraints = []

        self.solver   = pywraplp.Solver.CreateSolver(solver_id)
        if not self.solver:
            exit(0) #TODO: display error message and catch exception
        
        self.k = len(list(filter(lambda f_sv: f_sv!=0, F[self.source]))) #k = |{v in V : f(s,v)>0}| trivial lower bound, think of funnels

        self.edge_vars = []
        self.pi_vars   = []
        self.weights   = []
        self.path_vars = []

    def add_constraint(self, constraint):
        self.constraints.append(str(constraint))

    def clear(self):
        self.constraints = []
        self.edge_vars   = []
        self.pi_vars     = []
        self.weights     = []
        self.path_vars   = []
        self.solver.Clear()

    def display_stats(self):
        print("TODO display_stats")
        return

    def encode(self):
        for i in range(self.k): #DANGER the order of the variables in edge_varibles and self.pi_vars must be the same
            self.edge_vars += [ [self.solver.BoolVar('x_{}_{}_{}'.format(u,v,i))              for u in range(self.n) for v in range(self.n) if self.F[u][v]!=0] ]
            self.pi_vars   += [ [self.solver.IntVar(0, self.w_max,'p_{}_{}_{}'.format(u,v,i)) for u in range(self.n) for v in range(self.n) if self.F[u][v]!=0] ]
            self.weights   += [  self.solver.IntVar(1, self.w_max,'w_{}'.format(i)) ]
            self.path_vars += [ [self.solver.BoolVar('r_{}_{}'.format(i,j)) for j in range(len(self.R)) ]] #R=[ [(1,3),(3,5)], [(0,1)] ], 2 paths to cover: 1-3-5 and a single edge 0-1

        #3a, 3b
        for i in range(self.k):
            self.solver.Add( sum( filter(lambda edge: head(edge)==self.source, self.edge_vars[i]) ) == 1 )
            self.solver.Add( sum( filter(lambda edge: tail(edge)==self.target, self.edge_vars[i]) ) == 1 )
            
        #3c
        for i in range(self.k):
            for v in range(1,self.n-1): #find all wedges u->v->w for a fixed v excluding {s,t}
                self.solver.Add( sum( filter(lambda edge: tail(edge)==v, self.edge_vars[i]) ) - sum( filter(lambda edge: head(edge)==v, self.edge_vars[i]) ) == 0 )

        #5a
        for e in range(self.m):
            self.solver.Add( sum( self.pi_vars[i][e] for i in range(self.k)) == self.F[head(self.pi_vars[0][e])][tail(self.pi_vars[0][e])] ) #[0] bcause we just want the edge (u,v), the path is irrelevant

        #5b
        for i in range(self.k):
            for pi,edge in tuple(zip(self.pi_vars[i],self.edge_vars[i])):
                self.solver.Add( pi <= self.w_max * edge )
        #5c
        for i in range(self.k):
            for pi in self.pi_vars[i]:
                self.solver.Add( pi <= self.weights[i] )
        
        #5d
        for i in range(self.k):
            for pi,edge in tuple(zip(self.pi_vars[i],self.edge_vars[i])):
                self.solver.Add( pi >= self.weights[i] - (1 - edge) * self.w_max )
        
        #Example of a subpath constraint: R=[ [(1,3),(3,5)], [(0,1)] ], means that we have 2 paths to cover, the first one is 1-3-5. the second path is just a single edge 0-1
        def EncodeSubpathConstraints():
            #7a
            for i in range(self.k):
                for j in range(len(self.R)):
                    self.solver.Add( sum( filter( lambda e: exists_edge(e,self.R[j]), self.edge_vars[i] )) >= len(self.R[j])*self.path_vars[i][j] )
            #7b
            for r in range(len(self.R)):
                self.solver.Add( sum( self.path_vars[i][r] for i in range(self.k) ) >= 1 )
            
        def EncodeInexactFlow():
            #9a (to be exchanged with constraint 5a)
            for e in range(self.m):#[0] bcause we just want the edge (u,v), the path is irrelevant
                self.solver.Add( sum( self.pi_vars[i][e] for i in range(self.k) ) <= self.F_high[head(self.pi_vars[0][e])][self.tail(self.pi_vars[0][e])] )
                self.solver.Add( sum( self.pi_vars[i][e] for i in range(self.k) ) >= self.F_low [head(self.pi_vars[0][e])][self.tail(self.pi_vars[0][e])] )

        EncodeSubpathConstraints()

        #Add trivial objective function (we are just interested in deciding whether or not these constraints form a feasible region)
        self.solver.Minimize(1)

    def print_solution(self):
        print("Solution has size ", self.k, " with the following weight-path decomposition")
        for i in range(self.k):
            path = list(map(lambda e : tail(e), list(filter(lambda e : e.solution_value()==1, self.edge_vars[i]))))
            print(int(self.weights[i].solution_value()), path[:-1])

    def linear_search(self):
        while self.k <= min(self.m,self.f):

            self.encode()

            if self.solver.Solve() == pywraplp.Solver.OPTIMAL:
                self.print_solution()
                break

            self.clear()
            self.k += 1

def intron_to_matrix(intron_graph):
    F = [[0,6,7,0,0,0],
         [0,0,2,4,0,0],
         [0,0,0,9,0,0],
         [0,0,0,0,6,7],
         [0,0,0,0,0,6],
         [0,0,0,0,0,0]]
    #TODO
    #dictionary mapping nodes from intron_graph to nodes of ILP graph, i.e. a dict: NxN->N
    #same thing the other way around, i.e. a dict N->NxN
    return F


#Transform intron_graph into flow matrix F
# dont forget to add super source to every 0 in-degree vertex (leftmost guys in each layer) and a super target from every 0-outdegree vertex (rightmost guys in each layer) 
# F = intron_to_matrix(intron_graph)

# DAG may contain multiple connected components but thats fine, we can still add from super source to every guy... future think of possible opts (threads running in different components?) 


def Encode_ILP(intron_graph):

    F = intron_to_matrix(intron_graph)

    R = [[(1,3),(3,5)]]

    e = Enc(F,R=R)
    e.linear_search()

if __name__ == "__main__":
    Encode_ILP("")
    print('\nExiting ILP module now...')