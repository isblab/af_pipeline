import numpy as np
from scipy.spatial import distance_matrix

from typing import List


def str_join( prots: List, sep: str ):
	"""
	Given a list of strings, join them by the provided separator (sep).
	"""
	return f"{sep}".join( prots )


def str_split( prot_pair: str, sep: str ):
	"""
	Split the given string by the  provided separator (sep).
	"""
	return prot_pair.split( sep )


def get_distance_map( coords1: np.array, coords2: np.array ):
	"""
	Create an all-v-all distance map.
	"""
	distance_map = distance_matrix( coords1, coords2 )

	return distance_map


def get_contact_map( distance_map: np.array, contact_threshold: float ):
	"""
	Given the distance map, create a binary contact map by thresholding distances.
	"""
	contact_map = np.where( distance_map <= contact_threshold, 1, 0 )

	return contact_map


def get_interaction_map( coords1: np.array, coords2: np.array, contact_threshold: float, map_type: str ):
	"""
	Create an interaction map, given the input coordinates.
	Can create either of:
		Distance map
		Contact map
	"""
	distance_map = get_distance_map( coords1, coords2 )

	if map_type == "distance":
		return distance_map
	
	elif map_type == "contact":
		contact_map = get_contact_map( distance_map, contact_threshold )
		return contact_map
	
	else:
		raise Exception( "Invalid map_type specified..." )

