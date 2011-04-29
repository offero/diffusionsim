.. background

Background Information on DISim
*******************************

DISim is a program that simulates various Diffusion of Innovation models

1997 Adoption Threshold Model
=============================

Model with "Bandwagon Assessment" [AR1997]_

Original Threshold Model
------------------------

:math:`B_i,k = I_i + (A_i * P_k-1)`

:math:`B_i,k` is potential adopter i's "bandwagon assessment" of the innovation in cycle k

:math:`I_i` is potential adopter *i*'s individual assessment of the innovation's profitability
(given that it knows about the innovation?)

:math:`A_i * P_{k-1}` is the bandwagon pressure

:math:`P_{k-1}` is the information that creates bandwagon pressure after k-1 cycles

In fad models, :math:`P_{k-1}` is the proportion of adopters in a collectivity
 after cycle k-1 (from 1993 model). Num adopters / total population.

In learning models, :math:`P_{k-1}` is information about profitability; the average of
information across all adopters.

:math:`A_i` is the weight given to the information in :math:`P_{k-1}`. How much bandwagon pressure
affects agent *i*. Inversly related to trust in self-assessment. The higher the
ambiguity of assessment, the more the agent is influenced by bandwagon pressure
(relies on it).

If :math:`B_{i,k}` is evaluated >0, then adopt

1st, Model assumes global pressure (p. 296). Then, they modify to consider
social network structure (core-preiphery) (p. 297).

Modified Network Model
----------------------

:math:`P_{k-1}` replaced with :math:`P_{i,k-1}`

:math:`P_{i,k-1}` is the number of neighbors who have adopted / total population
Max bandwagon pressure becomes # neighbors / total population

"We assume that assessed returns are normally distributed..." (p. 297).

Trickle-down diffusion (focal to non-focal)
  Initiate by choosing 1 adopter in the focal stratum (core)
Trickle-up diffusion (non-focal to focal)
  1 adopter in the non-focal stratum (periphery)

"Other than the seed, for any potential adopter to adopt, it had to find out
information about the innovation through the network (i.e., communicate with
an adopter), and it had to find the innovation adoptable as a result of
finding out this information (i.e., positive Bi, k)." ([AR1997]_ pp 297)

1999 Network Model Adaptation
=============================

Refer to [RA1999]_

"Initial assessed profits and achieved profits are independently drawn from
the same normal distribution" ([RA1999]_ p. 368)

.. math::
	:label: 1999eq

	B_i,k = I_i * (A * \sum(r_i * D_{i,k-1})) + \left((1-A) *
	        \frac{\sum(p_i * D_{i,k-L-1})}{n_k-L-1}\right)

	0 \leq A \leq 1

:math:`p_i` = actual profits achieved by org i
*L* = # cycles required to transmit profitability information (lag)

