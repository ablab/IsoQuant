#To recall how the Sugiyama algorithm works https://blog.disy.net/sugiyama-method/

#There is a python package who does this for us...
#https://github.com/bdcht/grandalf

#A raw implementation
#https://github.com/gml4gtk/pysugiyama/blob/main/sugiyama.py

#R package that  implements many visualization algorithms
#https://igraph.org/r/doc/layout_with_sugiyama.html

from graphviz import Digraph

'''
TODO label the graph
def network2dot(G):
    dot = Digraph(format='pdf')
    dot.graph_attr['rankdir'] = 'LR' # Display the graph in landscape mode
    dot.node_attr['shape'] = 'rectangle' # Rectangle nodes
    
    for (u,v) in G.edges():
        att = G[u][v]
        dot.edge(str(u),str(v),label=f'{att["weight"]}')

    return dot
'''

def network2dot(G):
    dot = Digraph(format='pdf')
    dot.graph_attr['rankdir'] = 'LR' # Display the graph in landscape mode
    dot.node_attr['shape'] = 'rectangle' # Rectangle nodes

    (E,F) = G
    
    for (u,v) in E:
        dot.edge(str(u),str(v),label=str(F[(u,v)]))

    dot.render(directory='.', view=True)
    

def visualize(G):
    #should allow more arguments. e.g., a set of paths that we want to highlight in the drawing
    network2dot(G)