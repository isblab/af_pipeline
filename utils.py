import numpy as np
from scipy.spatial import distance_matrix


def get_distance_map( coords1, coords2 ):
	"""
	Create an all-v-all distance map.
	"""
	distance_map = distance_matrix( coords1, coords2 )

	return distance_map


def get_contact_map( distance_map, contact_threshold ):
	"""
	Given the distance map, create a binary contact map by thresholding distances.
	"""
	contact_map = np.where( distance_map <= contact_threshold, 1, 0 )

	return contact_map
