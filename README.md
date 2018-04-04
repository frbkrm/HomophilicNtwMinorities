# HomophilicNtwMinorities

The code generates networks with two groups that can vary in size. The network is generated using combination of preferential attachment and homophily. The homophily is a parameter that indicates the attractiveness between two nodes i and j, denoted as $h_{ij}$.

In this model, nodes have only two attributes, so they are assigned to two groups with different sizes. The homophily parameter alters between 0 to 1. To generate the network, you give these parameters as input:


N: total number of nodes in the network

m: minimum number of edges that an arrival node has

f_a : Minority fraction varies between 0 to 1.

$h_{ab}$: the probability of connection between group a and b

$h_{ba}$: the probability of connection between group b to a

The homophily parameters are ranging between 0 to 1. If homophily is symmetric: $h_{ab}$ = $h_{ba}$

For more information please check the paper https://arxiv.org/abs/1702.00150.

For visualization of this model check this repository: https://github.com/neuronphysics/Homophilic-Ntw-Viz

