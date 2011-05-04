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


def findWeaknessesAndPressurePoints(G, A_i, proportion=1/2, 
                                    targetSegment='core'):
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
                              and pressure points. 
    :return tuple: A tuple of 2 values, the first is the number of boundary
                   weaknesses and the second is the number of boundary 
                   pressure points.
    
    .. todo:: 
        Generalize segments; expand from just 'core':T/F attributes
        Create a generalized "segment" attribute of which "coreA", "periphA",
        etc. can be values.
    
    """
    countWeaknesses = 0
    countPPoints = 0
    # Segment 'A' nodes (targetSegment)
    A = [n for n in G.nodes() if G.node[n][targetSegment]]
    n_a = len(A) # number of nodes in A
    n_b = G.number_of_nodes() - n_a # number of nodes not in A
    
    for a_i in A:
        # The set of neighbors of a_i not in the same segment
        B_a_i = [n for n in G.neighbors(a_i) if not G.node[n][targetSegment]]        
        # Detect boundary weakness
        if len(B_a_i) > 0:            
            Bc_ik = G.node[a_i]['I'] + (A_i * (1/G.number_of_nodes()))
            if Bc_ik > 0: 
                countWeaknesses += 1
        # Detect pressure point
        if len(B_a_i) >= n_b*proportion:
            countPPoints += 1
    
    return (countWeaknesses, countPPoints)
