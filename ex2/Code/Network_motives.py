import numpy as np
import plotly.graph_objects as go
from scipy.stats import zscore, norm


def create_common_neighbours_matrix(matrix):
    """
    Create a matrix where each cell is the number of common neighbours between the i'th and j'th gene.
    :param matrix: The binarized matrix of interactions
    :return: The matrix of common neighbours
    """
    # Remove the diagonal, as we don't want to count the gene itself as a common neighbour
    np.fill_diagonal(matrix, 0)

    # Create the common neighbours matrix, by multiplying the matrix by itself
    common_neighbours_matrix = matrix @ matrix

    return common_neighbours_matrix


def create_any_interaction_matrix(matrix):
    """
    Create a matrix where each cell is the total number of interactions between the i'th and j'th gene, no matter the
    sign of the interaction, and who is the regulator and who is the regulated.
    :param matrix: The matrix of interactions
    :return: The matrix of any interactions
    """
    abs_matrix = np.abs(matrix)
    total_interactions = abs_matrix + abs_matrix.T

    # Binarize the matrix
    binarized_matrix = np.where(total_interactions > 0, 1, 0)
    return binarized_matrix


def check_neighbours_rule(matrix, domain):
    # Create the binarized matrix, which tells if there is any interaction between the i'th and j'th gene.
    binarized_matrix = create_any_interaction_matrix(matrix)
    common_neighbours_matrix = create_common_neighbours_matrix(binarized_matrix)
    # get the upper triangle of the interactions matrix, without the diagonal, to an array
    upper_triangle_interactions = binarized_matrix[np.triu_indices(binarized_matrix.shape[0], k=1)]
    # get the upper triangle of the neighbours matrix, without the diagonal, to an array

    upper_triangle_neighbours = common_neighbours_matrix[np.triu_indices(common_neighbours_matrix.shape[0], k=1)]
    # plot a histogram of the number of interaction genes, as a function of number of common neighbours, using plotly
    fig = go.Figure()
    fig.add_trace(go.Histogram(x=upper_triangle_neighbours, y=upper_triangle_interactions, histfunc="avg",
                               marker=dict(line=dict(color="black", width=2))))

    fig.update_layout(title=f"The Probability of Interactions between two {domain} as a Function of Common Neighbours",
                      xaxis_title="Number of common neighbours",
                      yaxis_title=f"Probability of interaction between two {domain}",
                      bargap=0.1)
    fig.write_html(f"../Figs/{domain}.html")
    fig.show()
    return common_neighbours_matrix


def count_motifs(matrix, verbose=False):
    abs_matrix = np.abs(matrix)
    total_interactions = abs_matrix + abs_matrix.T
    np.fill_diagonal(total_interactions, 0)
    mutual_interactions = np.sum(total_interactions == 2)
    if verbose:
        print("Number of mutual interactions in Matrix: ", mutual_interactions)

    # Extract coordinates of mutual interactions
    mutual_interactions_coordinates = np.argwhere(total_interactions == 2)
    # For each mutual interaction, iterate over all other nueurons, and check if it's regulated by the first neuron,
    # and regulates the second neuron

    motif_in_network = []
    for interaction in mutual_interactions_coordinates:
        Y = interaction[0]
        W = interaction[1]
        for X in range(celeg_matrix.shape[0]):
            if ((celeg_matrix[Y, X] != 0 and celeg_matrix[X, Y] == 0 and celeg_matrix[X, W] != 0)
                    and celeg_matrix[W, X] == 0) and X != W and X != Y:
                for Z in range(celeg_matrix.shape[0]):
                    if celeg_matrix[W, Z] == 0 and celeg_matrix[Z, W] == 0 and celeg_matrix[Z, Y] == 0 and\
                            celeg_matrix[Y, Z] != 0 and celeg_matrix[X, Z] == 0 and celeg_matrix[Z, X] == 0 and\
                            Z != W and Z != Y and Z != X:
                        motif_in_network.append([W, X, Y, Z])

    num_motifs = len(motif_in_network)
    if verbose:
        print("Number of Motifs in network: ", num_motifs)

    return num_motifs


if __name__ == '__main__':
    # load the coliAdj.txt to a numpy array
    ecoly_matrix = np.genfromtxt("../Resources/coliAdj.txt", delimiter=',', invalid_raise=False)
    celeg_matrix = np.genfromtxt("../Resources/elegansAdj.txt", delimiter=' ', invalid_raise=False, dtype=np.int32)

    # Q1
    print("Number of Transcription Factors: ", np.sum(np.count_nonzero(ecoly_matrix, axis=1) >= 1))

    # Q2
    # Count number of inhibitory interactions, and number of activatory interactions
    inhibitory = np.sum(ecoly_matrix < 0)
    activatory = np.sum(ecoly_matrix > 0)
    print("Number of inhibitory interactions: ", inhibitory)
    print("Number of activatory interactions: ", activatory)

    # Q3
    # Count number of non regulated genes, where the sum of it's column is 0
    non_regulated = np.sum(np.count_nonzero(ecoly_matrix, axis=0) == 0)
    print("Number of non regulated genes: ", non_regulated)

    # Q4a
    # check_neighbours_rule(
    #     ecoly_matrix,
    #     domain="E. coli genes"
    # )

    # Q5a
    abs_matrix = np.abs(ecoly_matrix)
    total_interactions = abs_matrix + abs_matrix.T
    np.fill_diagonal(total_interactions, 0)
    mutual_interactions = np.sum(total_interactions == 2)
    print("Number of mutual interactions in E. coly: ", mutual_interactions)

    # Q4b
    # check_neighbours_rule(
    #     celeg_matrix,
    #     domain="C. elegans neurons"
    # )

    #Q5b
    # count the number of mutual interactions in the celeg matrix
    num_orig_motifs = count_motifs(celeg_matrix, verbose=True)

    # randomize the celeg matrix 1000 times, and count the number of motifs in each randomization
    num_motifs = []
    for i in range(1000):
        print(f"Randomization number {i}")
        shuffled = np.random.permutation(celeg_matrix.flatten()).reshape(celeg_matrix.shape)
        num_motifs.append(count_motifs(shuffled))
        print(num_motifs[i])

    # Calculate the Z-score and p-value of the number of motifs in the original celeg matrix
    z_score = (num_orig_motifs - np.mean(num_motifs)) / np.std(num_motifs)
    p_value = 1 - norm.cdf(z_score)
    print("Z-score: ", z_score)
    print("p-value: ", p_value)













