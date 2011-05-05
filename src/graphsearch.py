'''
Created on May 4, 2011

:author: Christopher Kirkos

Theory
------

Boundary Pressure Points
++++++++++++++++++++++++

*"In sum, we define a boundary pressure point as a concentration of social 
ties linking potential adopters of an innovation in one segment of a 
network to a potential adopter in another segment of that network."
([AR1997]_ p. 300).* 

*"... in the second set of simulations, we operationalized boundary 
pressure points by counting each non-focal potential adopter that 
communicates with at least half of the focal potential adopters. We also 
tried proportions other than one half, and the results did not differ 
substantially." ([AR1997]_ p. 300).* 

Boundary Weaknesses
+++++++++++++++++++

*"...we operationalized boundary weaknesses by counting each non-focal 
potential adopter that satisfied two conditions: the potential adopter
had to communicate with a focal potential adopter, and it had to have
assessed profits high enough such that one adoption would create enough
impetus for this potential adopter to adopt" ([AR1997]_ p. 300).*

*"A boundary weakness occurs ... when potential adopter F both has ties
bridging two sides of a boundary and has a low adoption threshold. A single
adoption can cause such a weakly linked potential adopter to adopt..." 
([AR1997]_ p. 300).* 

*"In sum, we define a boundary weakness as a social tie linking a potential
adopter of an innovation in one segment of a network to a potential adopter,
in another segment of that network, who is highly predisposed to adopting
this innovation." ([AR1997]_ p. 300).* 


Implementation
--------------

'''

WPP_CACHE = {}
def findWeaknessesAndPressurePoints(G, A_i, proportion=1/2, 
                                    targetSegment='periphery',
                                    addGraphAttrs=True,
                                    ignoreCache=False, 
                                    clearCache=False):
    """Searches a graph for nodes that match the conditions for boundary
    weaknesses and boundary pressure points as given by [AR1997]_.
    
    :param networkx.Graph G: The networkx Graph object representing the network
    :param int A_i: The ambiguity level for this simulation to calculate 
                    potential bandwagon pressure. For boundary weakness 
                    calculation.
    :param float proportion: The proportion of nodes from an alternate segment 
                             B that are required to neighbor a given target 
                             node from segment A (where A <union> B = {}).
                             For pressure point calculation.
    :param str targetSegment: The attribute of the Graph `G` that identifies
                              the segment from which to determine weaknesses
                              and pressure points. Specify `periphery` for 
                              trickle down diffusion and `core` for trickle
                              up diffusion.
    :param bool addGraphAttrs: Augments the Graph `G` with node attributes
                               representing weaknesses and pressure points.
    :param bool ignoreCache: Cause function to recalculate for the given
                             Graph `G` and update the cache accordingly.
    :param bool clearCache: Reset the cache, deleting all existing entries.
    :return tuple: A tuple of 2 lists, the first list contains the node IDs 
                   that were identified as being boundary weaknesses, the
                   second contains node ID's of pressure points.
    
    .. todo:: 
        Generalize segments; expand from just 'core':T/F attributes
        Create a generalized "segment" attribute of which "coreA", "periphA",
        etc. can be values.
    
    .. todo::
        See if there's a way to automatically determine graph 'dirtyness' for
        cache invalidation.
    
    """
    # cache the result for multiple calls
    # WARNING: The results become invalid if the graph's edge configuration
    # changes!
    global WPP_CACHE
    if clearCache:
        WPP_CACHE.clear()
    
    cacheKey = (G,A_i,proportion,targetSegment)
    if not ignoreCache and cacheKey in WPP_CACHE:
        return WPP_CACHE[cacheKey]
    
    weakNodes=[]
    pressurePointNodes=[]
    # Segment 'A' nodes (targetSegment)
    A = [n for n in G.nodes() if targetSegment in G.node[n]['segments']]
    n_a = len(A) # number of nodes in A
    n_b = G.number_of_nodes() - n_a # number of nodes not in A
    
    for a_i in A:
        # The set of neighbors of a_i not in the same segment
        B_a_i = [n for n in G.neighbors(a_i) if targetSegment not in \
                                                    G.node[n]['segments']]        
        # Detect boundary weakness
        if len(B_a_i) > 0:            
            Bc_ik = G.node[a_i]['I'] + (A_i * (1/G.number_of_nodes()))
            if Bc_ik > 0: 
                weakNodes.append(a_i)
                if addGraphAttrs:
                    G.node[a_i]['weak']=True
        # Detect pressure point
        if len(B_a_i) >= n_b * proportion:
            pressurePointNodes.append(a_i)
            if addGraphAttrs:
                    G.node[a_i]['ppoint']=True
    
    WPP_CACHE[cacheKey] = (weakNodes, pressurePointNodes)
    
    return (weakNodes, pressurePointNodes) 

