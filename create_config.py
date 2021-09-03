"""
This python script can be used to create a configuration file with
which the C++ model 'games_on_graphs' can be parametrized.

Once the configuration file (by default config.txt) is created, the
model can be started as following (type in terminal):

./games_on_graphs -c config.txt
"""

import numpy
import igraph
import json


def getCommunityMeasures (E, c_i, c_j):
    """
    Given the community matrix, the function returns the number
    of internal and external links for the two modules i and j.
    :param E:
    :param i:
    :param j:
    :return:
    """

    L_j_in = E[c_j, c_j]
    L_i_in = E[c_i, c_i]
    L_j_out = numpy.sum(E[c_j,:]) - L_j_in
    L_i_out = numpy.sum(E[c_i,:]) - L_i_in

    return L_i_in, L_i_out, L_j_in, L_j_out

def modularityFromCommunityMatrix (E):
    """
    Computes the modularity from the community matrix. The matrix
    is assumed to be an upper triangular matrix with diagonal elements
    being the number of links within modules and non-diagonal elements
    being number of links between modules.

    The functions is implemented according to:
    [Generating graphs that approach a prescribed modularity,
     Computer Communications, 2013]

    :param E:
    :return:
    """

    E = numpy.triu(E)               # make upper triangular matrix
    c = float(numpy.shape(E)[0])    # number of modules
    L = float(numpy.sum(E))         # total number of links
    E_ = numpy.copy(E)
    numpy.fill_diagonal(E_, 0.0)
    L_inter = float(numpy.sum(E_))  # total number of inter module links

    sum = 0.0
    for i in range(0, int(c)):
        for j in range(i, int(c)):
            L_i_in, L_i_out, L_j_in, L_j_out = getCommunityMeasures(E, i, j)
            D_c_i = 2.0 * L_i_in + L_i_out  # sum of degrees in module i
            D_c_j = 2.0 * L_j_in + L_j_out  # sum of degrees in module j
            sum += numpy.power((D_c_i-D_c_j) / (2.0*L), 2.0)

    Q = 1.0 - 1.0/c - L_inter/L - 1.0/(2.0*c) * sum
    return Q

def constructModularGraphEdgeList(c, L_out, N_c, L_in=None):
    """
    Creates a graph with N_c nodes per module L_in edges within
    each module and L_out/(c-1) edges between each pair of modules.
    If L_in is None, then each module is fully connected.

    If the likelihood of a disconnected module is too high given L_in,
    an Exception is thrown.
    This likelihood is computed based on Erdos Reni model results:
    https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93R%C3%A9nyi_model

    :return:
    """

    c = int(c)
    L_out = int(L_out)
    N_c = int(N_c)
    N = N_c*c

    ############################################
    # Check preconditions
    ############################################
    if L_out % (c - 1) != 0:
        raise Exception("Unequal number of edges between modules.")

    if L_out % (c - 1) > (N_c*N_c):
        raise Exception("L_out is two high for the given N_c.")

    if not L_in is None:
        if L_in > (N_c*N_c-N_c)/2.0:
            raise Exception("L_in is too high.")
        # Compute edge probability (Erdos Reni) to check
        # whether L_in is too low (i.e. disconnected module is expected)
        epsilon = 0.01
        p = L_in / ((N_c*N_c-N_c)/2.0)
        if p < (1.0-epsilon) * numpy.log(1.0*N_c) / (1.0*N_c):
            raise Exception("L_in is too low. Disconnected module very likely.")

    # Create a graph with N nodes
    g = igraph.Graph()
    vertex_indices = range(0, N)
    g.add_vertices(vertex_indices)

    ############################################
    # Connect modules internally (L_in)
    ############################################
    edges_intra_module = []
    for c_i in range(c):

        # Vertex indices in module i. Create all possible node pairs
        # and store in list.
        vertex_indices_c_i = vertex_indices[c_i * N_c:c_i * N_c + N_c]
        edges_intra_module_c_i = []
        for k in range(N_c):
            for l in range(k + 1, N_c):
                edges_intra_module_c_i.append((vertex_indices_c_i[k], vertex_indices_c_i[l]))

        # If L_in is given shuffle edge list and take only L_in first ones
        if not L_in is None:
            numpy.random.shuffle(list(edges_intra_module_c_i))
            edges_intra_module_c_i = edges_intra_module_c_i[:L_in]

        edges_intra_module += edges_intra_module_c_i


    ############################################
    # Create L_out/(c-1) random edges between each
    # pair of modules
    ############################################
    edges_inter_module = []
    for c_i in range(c):
        for c_j in range(c_i + 1, c):

            # List of vertex indices of the two modules i and j
            vertex_indices_c_i = vertex_indices[c_i * N_c:c_i * N_c + N_c]
            vertex_indices_c_j = vertex_indices[c_j * N_c:c_j * N_c + N_c]

            # Each pair of vertex indices i,j can potentially be connected.
            # Hence, we have to distribute L_out / (c - 1) ones randomly in the
            # matrix with columns as j vertex indices and rows as i vertex indices.
            # Here a flat index of this matrix entries is created incrementing
            # row and column wise:
            #
            # 0 1 2
            # 3 4 5
            # 6 7 8
            # We then pick L_out / (c - 1) of these indices from a shuffled list
            # containing all N_c*N_c of them
            flat_ij_indices = range(0, N_c*N_c)
            numpy.random.shuffle(list(flat_ij_indices))

            for k in range(0, int(L_out / (c - 1))):
                flat_index = flat_ij_indices[k]
                i_index = int(flat_index / N_c)
                j_index = int(flat_index % N_c)
                edges_inter_module.append((vertex_indices_c_i[i_index], vertex_indices_c_j[j_index]))

    # Add the edge lists to the graph
    g.add_edges(edges_intra_module)
    g.add_edges(edges_inter_module)

    return g

def createNetwork(N, num_modules, L_in, L_out):
    E = numpy.zeros(shape=(num_modules, num_modules))
    E[:] = L_out / (num_modules - 1)
    numpy.fill_diagonal(E, L_in)
    Q = modularityFromCommunityMatrix(E)
    Q = Q

    N_c = N / num_modules
    print(E)
    g = constructModularGraphEdgeList(num_modules, L_out, N_c, L_in)
    clusters = g.clusters()
    if len(clusters) > 1:
        print("Graph is not connected")
        exit(1)

    return g

##########################################################
# Creates a configuration file for the C++ implementation
# of the mode (games_on_graphs)
##########################################################
config_file = 'config.txt'
N = 40
num_modules = 4
L_in = 40
L_out = 9
S = ['C'] * N   # Initial strategies of the N individuals (C=Cooperator, D=Defector)
S[0] = 'D'
T = 100
b = 1.5
sampleTime = 10
randomSeed = 100

# Create a network with N nodes, num_module modules,
# with each module having L_in edges, and L_out edges
# running between each pair of module.
graph = createNetwork(N, num_modules, L_in, L_out)

# Create a GML file of the graph
tempFile = "gml_graph.conf"
graph.write_gml(tempFile)
file = open(tempFile, 'r')
graphGML = file.read()
file.close()

simulationConfiguration = {
    "NETWORK": graphGML,
    "STRATEGIES": "".join(S),
    "STRATEGIES_ENFORCED" : "",
    "B": b,
    "T": T,
    "SAMPLE_TIME": sampleTime,
    "RANDOM_SEED": randomSeed
}

with open(config_file, 'w') as outfile:
    json.dump(simulationConfiguration, outfile)












