'''
Created on May 4, 2011

:author: Christopher Kirkos
'''

def findBoundaryPressurePoints(G, proportion=1/2):
    """
    
    *"In sum, we define a boundary pressure point as a concentration of social 
    ties linking potential adopters of an innovation in one segment of a 
    network to a potential adopter in another segment of that network."
    ([AR1997]_ p. 300).* 

    *"... in the second set of simulations, we operationalized boundary 
    pressure points by counting each non-focal potential adopter that 
    communicates with at least half of the focal potential adopters. We also 
    tried proportions other than one half, and the results did not differ 
    substantially." ([AR1997]_ p. 300).* 
    
    n_a, n_b are the number of nodes in segment a and b, respectively
    proportion = 1/2
    count = 0
    pressurePoints = []
    for each node a_i in segment A
        B_a_i is the subset of nodes in segment B that are neighbors of a_i
        if |B_a_i| >= n_b*proportion # comm. with at least half... 
            count += 1
            pressurePoints.append(a_i)
    
    TODO: Generalize segments; expand from just 'core':T/F attributes
          Create a generalized "segment" attribute of which "coreA", "periphA"
          etc. can be values.
    
    """
    
    count = 0
    # Segment 'A' nodes
    A = [n for n in G.nodes() if G.node[n]['core']]
    n_a = len(A) # number of nodes in A
    n_b = G.number_of_nodes() - n_a # number of nodes not in A
    
    for a_i in A:
        # The set of neighbors of a_i not in the same segment
        B_a_i = [n for n in G.neighbors(a_i) if not G.node[n]['core']]
        if len(B_a_i) >= n_b*proportion:
            count += 1
    
    return count

def findBoundaryWeaknesses(G):
    """Searches a graph for nodes that match the condition of being a boundary
    weakness as given by [AR1997]_.
    
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
    
    Pseudocode
    Directed, flow from segment A to segment B
    
    count = 0
    weaknesses = [] # list of nodes that are boundary weaknesses
    for each node a_i in segment A
        B_a_i is the subset of nodes in segment B that are neighbors of a_i 
        if exists at least 1 node b_i in B_a_i:
            calculate constrained B_i,k -> Bc_i,k
            Bc_i,k = I_i + (A_i * 1/n)
            if Bc_i,k > 0:
                weaknesses.append(a_i)
                
    
    len(B_a_i)
    
    """
    
    
