# Libraries
import numpy as np

def cost_function(visibility, time_visibility, distance, sunlight, alpha:float, beta:float, gamma:float, total_satellites:int):
    """
    Calculates the cost matrix (the weight of each edge in an undirected graph representing a satellite network where edges are potential ISLs and nodes are satellites).

    :param visibility:
    :param time_visibility:
    :param distance:
    :param sunlight:
    :param alpha:
    :param beta:
    :param gamma:
    :param total_satellites:
    :return:
    """
    # Where satellites are not visible to one another, set the cost as infinity (represented by -1), otherwise 0
    cost_matrix = np.where(visibility == 1, np.zeros((total_satellites, total_satellites)), -1)

    # Min-Max Scale/Normalise distances to ensure equal consideration of distance and other metrics included in cost
    min_dist = np.nanmin(np.where(distance > 0, distance, np.nan))
    max_dist = np.max(distance)

    if np.isnan(min_dist):
        raise ValueError("Minimum distance between satellites must be greater than 0.")
    elif max_dist == -1:
        raise ValueError("Maximum distance between satellites must be greater than 0.")

    # Below ensures no divide by zero error (unlikely all distances are the same, but possible)
    if min_dist == max_dist:
        distance /= max_dist
    else:
        # Min-Max Scale
        distance = (distance - min_dist) / (max_dist - min_dist)

    # Calculate costs/weights according to cost function (included in paper) - added  1 to time_visibility to ensure no
    # divide by 0 error. Gamma is probability of satellite failure due to solar flares - 0 if in Earth's shadow,
    # otherwise gamma (gamma could be found via deep learning image classification of pictures of sun)
    cost_matrix = np.where(visibility == 0, cost_matrix, (alpha * (1 / (time_visibility + 1))) + (beta * distance) + (gamma * sunlight))

    return cost_matrix

# References
# Adding 1 to Numpy Array - https://www.reddit.com/r/learnpython/comments/12hsf8k/trying_to_add_1_to_the_element_at_a_certain_index/
