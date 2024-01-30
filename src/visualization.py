from graphviz import Digraph

def network2dot(G, paths=[]):
    dot = Digraph(format='pdf')
    dot.graph_attr['rankdir'] = 'LR'        # Display the graph in landscape mode
    dot.node_attr['shape']    = 'rectangle' # Rectangle nodes

    E,F    = G
    colors = ['red','blue','green','purple','brown','cyan','yellow','pink','grey']

    for (u,v) in E:
        dot.edge(str(u),str(v),label=str(F[(u,v)]))

    for index, (weight,slack,path) in enumerate(paths):
        pathColor = colors[index % len(colors)]
        for i in range(len(path)-1):
            dot.edge(str(path[i]), str(path[i+1]), label=str(weight) + ", " + str(slack), fontcolor=pathColor, color=pathColor, penwidth='2.0')
        if len(path) == 1:
            dot.node(str(path[0]), color=pathColor, penwidth='2.0')
           
    dot.render(directory='.', view=True)
   

def visualize(G, paths=[]):
    network2dot(G, paths)