from src.intron_graph import VERTEX_polya, VERTEX_polyt, VERTEX_read_end, VERTEX_read_start
from src.visualization import visualize
from collections import defaultdict
from gurobipy import GRB
import gurobipy as gp
import logging

logger = logging.getLogger('IsoQuant')

TOLERANCE = 0.1  #tolerance allowed for Gurobi numerical values

def tail(grb_edge): return int(grb_edge.split("[")[1].split(",")[0])
def head(grb_edge): return int(grb_edge.split("[")[1].split(",")[1])
def path(grb_edge): return int(grb_edge.split("[")[1].split(",")[2])

class Encode:

    def __init__(self,n,E,F,epsilon,timeout,threads,R=[]):
        self.n          = n
        self.m          = len(E)
        self.source     = 0
        self.target     = self.n-1
        self.E          = E
        self.F          = F
        self.R          = R
        self.epsilon    = epsilon #we could have a function that finds a good epsilon based on features of the data
        self.timeout    = timeout
        self.threads    = threads
        self.w_max      = 0
        self.edge_vars  = {}
        self.phi_vars   = {}
        self.gam_vars   = {}
        self.weights    = {}
        self.slacks     = {}
        self.spc_vars   = {}
        self.model      = self.create_solver()
        
        source_outdegree_nonzero = 0
        target_indegree_nonzero  = 0
        for (u,v) in self.E:
            if u==self.source and F[(u,v)]!=0:
                source_outdegree_nonzero += 1
            if v==self.target and F[(u,v)]!=0:
                target_indegree_nonzero += 1
            self.w_max = max(self.w_max,self.F[(u,v)])
        self.k = max(source_outdegree_nonzero, target_indegree_nonzero) #a lower bound for the largest edge antichain of the graph

    def create_solver(self):
        env = gp.Env(empty=True) #https://docs.gurobi.com/projects/optimizer/en/current/reference/python/env.html#Env.Env
        env.setParam('OutputFlag'   ,              0)
        env.setParam('LogToConsole' ,              0)
        env.setParam('TimeLimit'    ,   self.timeout)
        env.setParam('Threads'      ,   self.threads)
        env.start()
        model = gp.Model("MFD_Robust", env=env)
        if not model:
            logger.error("FATAL, could not create Gurobi model")
            exit(0) #TODO: use try except
        return model

    def clear(self):
        #self.env.close()
        #self.model.close()
        self.model       = self.create_solver()
        self.edge_vars   = {}
        self.phi_vars    = {}
        self.gam_vars    = {}
        self.weights     = {}
        self.slacks      = {}
        self.spc_vars    = {}

    def display_stats(self):
        print("###################################\nSTATS####################################\n")
        self.model.display()
        self.model.printStats()
        print(self.model.numVars)
        for v in self.model.getVars():
            print(v,v.VarName)

        print("#####################################\n######################################\n")
        return

    def solve(self):
        self.model.optimize()

    #returns a value in ]0,1] that depends on the position of the edge in the graph (because of non uniform coverages, etc), 1 is maximum confidence
    def confidence(self, edge):
        return 1

    def encode(self):

        # Create variables
        edge_indexes    = [ (u,v,i) for i in range(self.k) for (u, v) in self.E        ]
        path_indexes    = [ (    i) for i in range(self.k)                             ]
        subpath_indexes = [ (i,j  ) for i in range(self.k) for j in range(len(self.R)) ]

        #try to replace all the CONTINUOUS variables with INTEGER variables
        self.edge_vars = self.model.addVars(   edge_indexes, vtype=GRB.BINARY    ,  name='e'                     )
        self.spc_vars  = self.model.addVars(subpath_indexes, vtype=GRB.BINARY    ,  name='r'                     )
        self.phi_vars  = self.model.addVars(   edge_indexes, vtype=GRB.CONTINUOUS,  name='p', lb=0, ub=self.w_max)
        self.gam_vars  = self.model.addVars(   edge_indexes, vtype=GRB.CONTINUOUS,  name='g', lb=0, ub=self.w_max)
        self.weights   = self.model.addVars(   path_indexes, vtype=GRB.CONTINUOUS,  name='w', lb=1, ub=self.w_max) #might be worth to test lower bound to zero and assess the consequences on ObjVal
        self.slacks    = self.model.addVars(   path_indexes, vtype=GRB.CONTINUOUS,  name='s', lb=0, ub=self.w_max)

        #The identifiers of the constraints come from https://www.biorxiv.org/content/10.1101/2023.03.20.533019v1.full.pdf page 13

        for i in range(self.k):
            self.model.addConstr( self.edge_vars.sum(self.source,'*',i) == 1, "14a_i={}".format(i) )
            self.model.addConstr( self.edge_vars.sum('*',self.target,i) == 1, "14b_i={}".format(i) )

        for i in range(self.k):
            for v in range(1,self.n-1): #find all wedges u->v->w for v in V\{s,t}
                self.model.addConstr( self.edge_vars.sum('*',v,i) - self.edge_vars.sum(v,'*',i) == 0, "14c_v={}_i={}".format(v,i) )

        for (u,v) in self.E:
            f_uv = self.F[(u,v)]
            phi_sum = self.phi_vars.sum(u,v,'*')
            gam_sum = self.gam_vars.sum(u,v,'*')
            self.model.addConstr( f_uv - phi_sum <=  gam_sum / self.confidence((u,v)), "14d_u={}_v={}".format(u,v) )
            self.model.addConstr( f_uv - phi_sum >= -gam_sum / self.confidence((u,v)), "14e_u={}_v={}".format(u,v) )
            for i in range(self.k):
                self.model.addConstr( self.phi_vars[u,v,i] <= self.w_max * self.edge_vars[u,v,i]                        , "14f_u={}_v={}_i={}".format(u,v,i) ) 
                self.model.addConstr( self.gam_vars[u,v,i] <= self.w_max * self.edge_vars[u,v,i]                        , "14i_u={}_v={}_i={}".format(u,v,i) ) #FIXME w_max vs M
                self.model.addConstr( self.phi_vars[u,v,i] <= self.weights[i]                                           , "14g_u={}_v={}_i={}".format(u,v,i) )
                self.model.addConstr( self.gam_vars[u,v,i] <= self.slacks [i]                                           , "14j_u={}_v={}_i={}".format(u,v,i) )
                self.model.addConstr( self.phi_vars[u,v,i] >= self.weights[i] - (1 - self.edge_vars[u,v,i]) * self.w_max, "14h_u={}_v={}_i={}".format(u,v,i) )
                self.model.addConstr( self.gam_vars[u,v,i] >= self.slacks [i] - (1 - self.edge_vars[u,v,i]) * self.w_max, "14k_u={}_v={}_i={}".format(u,v,i) ) #FIXME w_max vs M

        #Example of a subpath constraint: R=[ [(1,3),(3,5)], [(0,1)] ], means that we have 2 paths to cover, the first one is 1-3-5. the second path is just a single edge 0-1
        def EncodeSubpathConstraints():
            for i in range(self.k):
                for j in range(len(self.R)):
                    edgevars_on_subpath = list(map(lambda e: self.edge_vars[e[0],e[1],i], self.R[j]))
                    self.model.addConstr( sum(edgevars_on_subpath) >= len(self.R[j]) * self.spc_vars[i,j] )
            for j in range(len(self.R)):
                self.model.addConstr( self.spc_vars.sum('*',j) >= 1 )

        if self.R!=[]:
            EncodeSubpathConstraints()

        self.model.setObjective( self.slacks.sum(), GRB.MINIMIZE )


    def print_solution(self,solution):
        opt, slack, paths = solution
        print("\n#####SOLUTION#####\n","> FD size   :",   opt,"\n> Slack sum :", slack,"\n> Weight-Slack-Path decomposition:")
        for p in paths:
            print(*p)
        visualize((self.E,self.F),paths)

    ''' About the "build_solution" function:
    paths in the solution are of the form [ (weight_0, slack_0, p_0), ... (weight_i, slack_i, p_i), ...]
    where p_i is of the form [v_1,v_2,...,v_k], where source and target are not included (see build_solution)
    '''    
    def build_solution(self):
        paths = []
        for i in range(self.k):
            path = []
            u    = self.source
            while u != self.target:
                edges = list(filter(lambda grb_var : grb_var.X>1-TOLERANCE , self.edge_vars.select(u,'*',i) ))
                assert(len(edges)==1)
                v = head(edges[0].VarName)
                path.append(v)
                u = v
            paths.append( (self.weights[i].X, self.slacks[i].X, path[:-1]) )

        return (self.k, self.model.ObjVal, paths)

    def get_paths_from_solution(self, solution):
        _,_,p = solution
        self.print_solution(solution)
        paths,weights = list(map(lambda x : x[2], p)), list(map(lambda x : x[0], p))
        return (paths,weights)
    
    def linear_search(self):
        print("############# Feasibility #############")
        logger.info(">>> Feasibility, starting with " + str(self.k))
        #Feasibility: Find first feasible solution (there always exists one) and clear the solver for next stage
        while True:
            self.encode()
            self.solve()
            if self.model.status == GRB.OPTIMAL:
                previous_slack = self.model.ObjVal
                solution       = self.build_solution()
                self.clear()
                self.k += 1
                break
            self.clear()
            self.k += 1

        if previous_slack == 0:
            logger.info(">>> No optimization required, found optimal solution in feasibility step with width " + str(self.k-1))
            return self.get_paths_from_solution(solution)

        print("#############Optimality#############")
        logger.info(">>> Optimality, starting with " + str(self.k))
        #Optimality: Find the k for which the difference in slacks of two consecutive iterations becomes sufficiently small
        while self.k <= self.m:
            self.encode()
            self.solve()

            if self.model.status == GRB.TIME_LIMIT: #TODO: use try except
                raise GRB_TimeOut("ilp.Robust.optimize_linear while optimizing")

            ''' NOTE: there is an interplay between the lower bound of the weight variables and this block of code
            if self.model.ObjVal >= previous_slack: #if we do not improve by allowing more paths we stop
                self.final_k = self.k-1
                break
            '''

            solution = self.build_solution()

            if self.model.ObjVal==0: #if we found a perfect solution we stop
                return self.get_paths_from_solution(solution)
            
            assert(previous_slack != 0)

            #The stopping condition can also be (previous_slack - self.model.ObjVal < self.epsilon), or perhaps even another thing. Can play with this
            if 1-self.model.ObjVal/previous_slack < self.epsilon: # if the relative improvement is small we also stop
                return self.get_paths_from_solution(solution)
            
            previous_slack = self.model.ObjVal
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
class Intron2Graph:
    
    def __init__(self, intron_graph):
        
        self.intron2vertex = dict()
        self.vertex2intron = dict()
        self.vertex_id = 1

        # add all intron vertices
        for intron in intron_graph.intron_collector.clustered_introns.keys():
            self.intron2vertex[intron] = self.vertex_id
            self.vertex2intron[self.vertex_id] = intron
            self.vertex_id += 1

        source = 0
        target = self.vertex_id

        # create edges
        self.edge_list = []
        edge_set = set()
        starting_introns = defaultdict(int)
        for intron in intron_graph.incoming_edges.keys():
            for preceeding_intron in intron_graph.incoming_edges[intron]:
                self.edge_list.append((self.intron2vertex[preceeding_intron], self.intron2vertex[intron]))
                edge_set.add((self.intron2vertex[preceeding_intron], self.intron2vertex[intron]))
                if preceeding_intron[0] in [VERTEX_polyt, VERTEX_read_start]:
                    starting_introns[preceeding_intron] += intron_graph.edge_weights[(preceeding_intron, intron)]

        terminal_introns = defaultdict(int)
        for intron in intron_graph.outgoing_edges.keys():
            for subsequent_intron in intron_graph.outgoing_edges[intron]:
                if subsequent_intron[0] in [VERTEX_polya, VERTEX_read_end]:
                    terminal_introns[subsequent_intron] += intron_graph.edge_weights[(intron, subsequent_intron)]

        self.flow_dict = defaultdict(int)
        for intron in intron_graph.incoming_edges.keys():
            for preceeding_intron in intron_graph.incoming_edges[intron]:
                u = self.intron2vertex[preceeding_intron]
                v = self.intron2vertex[intron]
                self.flow_dict[(u, v)] = intron_graph.edge_weights[(preceeding_intron, intron)]

        # add connection to super source and total weight
        for starting_intron in starting_introns.keys():
            starting_vertex = self.intron2vertex[starting_intron]
            self.edge_list.append((source, starting_vertex))
            edge_set.add((source, starting_vertex))
            self.flow_dict[(source, starting_vertex)] = starting_introns[starting_intron]

        for terminal_intron in terminal_introns.keys():
            terminal_vertex = self.intron2vertex[terminal_intron]
            self.edge_list.append((terminal_vertex, target))
            edge_set.add((terminal_vertex, target))
            self.flow_dict[(terminal_vertex, target)] = terminal_introns[terminal_intron]

        #FIXME edge_set is not necessary
        assert len(self.edge_list) == len(edge_set)
        assert len(self.edge_list) == len(self.flow_dict)

        self.vertex_id += 1
        logger.debug(self.vertex2intron)
        logger.debug(self.flow_dict)
        #return self.vertex_id+1, self.edge_list, self.flow_dict

    def transcript_to_path(self, transcript):
        vertex_list = list(map(lambda t: self.intron2vertex[t], transcript))
        return list(zip(vertex_list[:-1], vertex_list[1:]))

    def transcripts_to_paths(self, transcripts):
        return list(map(self.transcript_to_path, transcripts))

    def path_to_transcript(self,path):
        return list(map(lambda v : self.vertex2intron[v], path))
    
    def paths_to_transcripts(self, paths):
        return list(map(self.path_to_transcript, paths))


class GRB_TimeOut(Exception):
    def __init__(self, message:str):
        super(GRB_TimeOut, self).__init__('Gurobi TimeOut: ' + message)


def Encode_ILP(intron_graph, transcripts_constraints=[], epsilon=0.25, timeout=300, threads=5):

    #check if graph is empty
    if len(intron_graph.outgoing_edges)==0:
        logger.warning("Empty edge list in graph.")
        return []

    g = Intron2Graph(intron_graph)

    path_constraints = g.transcripts_to_paths(transcripts_constraints)

    n,E,F = g.vertex_id, g.edge_list, g.flow_dict

    visualize((E,F))
    e = Encode(n, E, F, epsilon, timeout, threads, path_constraints)
    e.encode()
    paths,weights = e.linear_search()
    
    transcripts = g.paths_to_transcripts(paths)
    weighted_transcripts = list(zip(transcripts,weights))
    
    return weighted_transcripts
