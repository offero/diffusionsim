# -*- coding: utf-8 -*-
"""
Module to generate graph structures.

:Author: Christopher Kirkos

Implementation
--------------
"""

from __future__ import division

import random
from random import sample
import networkx as nx
from itertools import combinations, chain
import pylab
from pylab import plt
from warnings import filterwarnings
from exceptions import RuntimeWarning

class DINetworkGenerator(object):
    """Abstract base class for DISim network generators. These objects
    are iterators that return graph objects on each iteration."""
    
    def __init__(self, n):
        """Construct the network generator object.
        
        :param integer n: The number of nodes in the network
        """
        self.n = n

    def __iter__(self):
        return self
    
    def next(self):
        raise NotImplementedError("The next() method should be overridden by"
                                  "a subclass.")


class DICorePeriphNxGenerator(DINetworkGenerator):
    """Generates a core-periphery network like the one discussed in [AR1997]_ 
    using NetworkX [HSS2008]_.
    
    `"First set of simulations ... core-perphiphery networks with fully-linked
    cores ... network density was varied by varying the number of network ties 
    beyond the core"` ([AR1997]_ p. 297).
    
    The core of the network is completely connected and the edges with/between
    the periphery are generated randomly.
    """
    
    def __init__(self, numCoreNodes, numPeriphNodes, pties, seed=None, 
                 *args, **kwargs):
        """Construct the network generator object.
        
        :param integer numCoreNodes: The number of nodes in the Core (>0).
        :param integer numPeriphNodes: The number of nodes in the Periphery (>0).
        :param integer pties: The number of additional ties to generate in 
                              the periphery.
        :param integer seed: A number to seed the random number generator.
                             (Optional)
        """
        assert(numCoreNodes>=0 and numPeriphNodes>=0)        
        n = numCoreNodes + numPeriphNodes
        super(DICorePeriphNxGenerator, self).__init__(n, *args, **kwargs)
        self.numCoreNodes = numCoreNodes
        self.numPeriphNodes = numPeriphNodes
        self.pties = pties
        
        # customize random number generator seed value if desired
        if not seed is None:
            random.seed(seed)
        
    def next(self):
        return generateARCorePeriph(self.numCoreNodes, self.numPeriphNodes, 
                                    self.pties)
        

def dissimilarProduct(A,B):
    """Generator for all combinations of 2 lists where the items are not the
    same."""
    return ((x,y) for x in A for y in B if x!=y)      

def setDefaultNodeAttrs(G):
    """Helper function to set default attributes on a new network. This 
    funciton modifies the graph itself.
    
    :param networkx.Graph G: A networkx Graph object.
    """
    for a in G.nodes():
        G.node[a]['adopted'] = False
        G.node[a]['influence'] = []


def generateARCorePeriph(numCoreNodes, numPeriphNodes, pties, show=False):
    """Generates a core-periphery network like the one discussed in [AR1997]_ 
    using NetworkX [HSS2008]_. 
    
    `"First set of simulations ... core-perphiphery networks with fully-linked
    cores ... network density was varied by varying the number of network ties 
    beyond the core"` ([AR1997]_ p. 297).
    
    :param integer numCoreNodes: The number of nodes in the Core (>0).
    :param integer numPeriphNodes: The number of nodes in the Periphery (>0).
    :param integer pties: The number of additional ties to generate in 
                          the periphery.
    """
    assert(numCoreNodes>=0 and numPeriphNodes>=0)
    
    # total number of nodes (n) in graph
    n = numCoreNodes + numPeriphNodes
    
    # Determine the number of core nodes, rounded to nearest int
    #core = int(round(n*cpratio))
    coreNodes = range(numCoreNodes)
    periphNodes = range(numCoreNodes, n)
    
    # Construct initial core network 
    #G = generators.complete_graph(core)
    # manually generate complete multigraph
    G = nx.empty_graph(numCoreNodes, create_using=nx.MultiGraph())
    G.add_edges_from( combinations(coreNodes,2) )
    for a in G.nodes():
        G.node[a]['core']=True

    G.add_nodes_from([(pn,{'core':False}) for pn in periphNodes])
    G.name="random core-periphery(%s)"%(n)

    # iterator for all periph to core edges
    pcedges = dissimilarProduct(periphNodes, coreNodes)
    # iterator for all periph to periph edges
    ppedges = dissimilarProduct(periphNodes, periphNodes)
    
    allPotentialEdges = chain(pcedges, ppedges)
        #sampEdges = sample(tuple(pcedges), pcties)
        #sampEdges.extend(sample(tuple(ppedges), ppties))
    
    # random sampling of edges to add
    sampEdges = sample(tuple(allPotentialEdges), pties)
    
    # add extra non-core (peripheral) edges to the network
    G.add_edges_from(sampEdges)
    
    if show:
        nx.draw(G)
        pylab.show()

    return G

def drawAdoptionNetworkGV(G, writeFile=None, writePng=None):
    """Generates the GraphViz adoption network. Optionally writes the output
    to DOT and/or PNG files.
    
    This function expected the node attributes 'adopted' and 'influence' to be 
    pre-populated. If the 'adopted' attribute is True, the node is given a 
    different color, showing adoption visually. If the 'influence' attribute 
    contains a list of nodes that played a role in a given node's adoption,
    then the edges from those nodes to the given target node is highlighted.
    
    :param networkx.MultiGraph G: The networkx MultiGraph object. This object 
                                  should be pre-populated with node attributes.
                                  A MultiGraph is needed because this function
                                  augments the graph with additional edges to
                                  represent information/influence flow visually.
    :param string writeDot: The filename/path to which to save the DOT file.
    :param String writePng: The filename/path to which to save the PNG file.
    :return: pygraphviz.AGraph object augmented with style attributes.
    """
    
    # pygraphviz calls the graphviz subprocess and issues a warning from its
    # output about node size. This output kills the console.
    filterwarnings(action="ignore", category=RuntimeWarning, 
                            module="agraph")
    
    colorAdopted = "dodgerblue"
    colorNonAdopted = "firebrick1"
    
    # global default graph attributes
    gvGraph = nx.to_agraph(G)
    gvGraph.graph_attr['splines']="true"
    #gvGraph.graph_attr['size']="6,6!"
    
    # global default node attributes
    gvGraph.node_attr['shape']='circle'
    gvGraph.node_attr['fixedsize']='true'
    gvGraph.node_attr['width']='0.25'
    gvGraph.node_attr['style']='filled'
    gvGraph.node_attr['fontcolor']='white'
    gvGraph.node_attr['fontsize']='8'
    
    # global default edge attributes
    gvGraph.edge_attr['dir']="none"
    gvGraph.edge_attr['weight']="1"
    gvGraph.edge_attr['color']="gray35"
    
    coreNodes = []
    
    # set custom colors for nodes and edges
    for node in gvGraph.nodes():    
        # set adopted node color
        node.attr['fillcolor'] = colorAdopted \
                                if node.attr['adopted'] == "True" \
                                else colorNonAdopted
        # add to list of core nodes
        if node.attr['core'] == "True":
            coreNodes.append(node)
            
        # set edge color of node influences
        # ni is a string from a list, could be '[]', or '[1,2,3]'
        # strip first and last characters, the brackets, from the string 
        influenceNodes = node.attr['influence'][1:-1]
        if len(influenceNodes) > 0:
            # turn string to list of individual node ids, iterate through list
            # TODO: find a better way, maybe use the original networkx.Graph
            influenceNodes = influenceNodes.split(', ')
            for ni in influenceNodes:
                # node influenced by ni
                #edge = gvGraph.get_edge(ni, node)
                #edge.attr['dir']='forward'
                #edge.attr['color']="#1E90FF8F"
                #edge.attr['penwidth']="4"
                gvGraph.add_edge(ni, node, 
                              weight="1",
                              dir="forward", 
                              color="#1E90FFAF", 
                              penwidth="4")
    
    # creating a cluster subgraph will cause graphviz to group them visually
    coreGraph = gvGraph.add_subgraph(coreNodes, "clusterCoreNodes")
    # color nodes of the core a different border
    #coreGraph.node_attr['color']="yellowgreen"
    for node in coreGraph.nodes():
        node.attr['color']="yellowgreen"
    
    for edge in coreGraph.edges():
        edge.attr['len']='1.5'
    
    if writeFile != None:
        gvGraph.write(writeFile)
    if writePng != None:
        gvGraph.draw(writePng, 'png', 'neato')
    
    return gvGraph

def drawAdoptionNetworkMPL(G, fnum=1, show=False, writeFile=None):
    """Draws the network to matplotlib, coloring the nodes based on adoption. 
    Looks for the node attribute 'adopted'. If the attribute is True, colors 
    the node a different color, showing adoption visually. This function assumes
    that the node attributes have been pre-populated.
    
    :param networkx.Graph G: Any NetworkX Graph object.
    :param int fnum: The matplotlib figure number. Defaults to 1.
    :param bool show: 
    :param string writeFile: A filename/path to save the figure image. If not
                             specified, no output file is written.
    """
    Gclean = G.subgraph([n for n in G.nodes() if n not in nx.isolates(G)])
    plt.figure(num=fnum, figsize=(6,6))
    # clear figure
    plt.clf()
    
    # Blue ('b') node color for adopters, red ('r') for non-adopters. 
    nodecolors = ['b' if Gclean.node[n]['adopted'] else 'r' \
                  for n in Gclean.nodes()]
    layout = nx.spring_layout(Gclean)
    nx.draw_networkx_nodes(Gclean, layout, node_size=80, 
                           nodelist=Gclean.nodes(), 
                           node_color=nodecolors)
    nx.draw_networkx_edges(Gclean, layout, alpha=0.5) # width=4
    
    # TODO: Draw labels of Ii values. Maybe vary size of node.
    # TODO: Color edges blue based on influences from neighbors
    
    influenceEdges = []
    for a in Gclean.nodes():
        for n in Gclean.node[a]['influence']:
            influenceEdges.append((a,n))
    
    if len(influenceEdges)>0:
        nx.draw_networkx_edges(Gclean, layout, alpha=0.5, width=5,
                               edgelist=influenceEdges,
                               edge_color=['b']*len(influenceEdges))
    
    #some extra space around figure
    plt.xlim(-0.05,1.05)
    plt.ylim(-0.05,1.05)
    plt.axis('off')
    
    if writeFile != None:
        plt.savefig(writeFile)
    
    if show:
        plt.show()


