import numpy as np

from typing import Dict

from parser import AfParser, ResidueSelect
from utils import get_distance_map, get_contact_map

"""
This script will contain modules for performing analysis for the AF2/3 prediction.
Define a class to perform housekeeping jobs like parsing, extracting coords, plddt, pae, etc.
Derive child classes for downstream analysis.
"""

class Initialize( AfParser ):
	def __init__( self, struct_file_path: str, data_file_path: str ):
		super().__init__( struct_file_path, data_file_path )
		# AF2/3 structure file path.
		self.struct_file_path = struct_file_path
		# AF2/3 structure data path.
		self.data_file_path = data_file_path

		self.get_attributes()


	def get_attributes( self ):
		"""
		Extract the following from the input files:
			1. Ca coordinates.
			2. Ca pLDDT.
			3. Average PAE matrix.
		"""
		# Residue positions of all residues for each chain.
		self.res_dict = self.get_residue_positions()
		self.lengths_dict = self.get_chain_lengths( self.res_dict )
		# Ca-coords of all residues for each chain.
		self.coords_dict = self.get_ca_coordinates()
		# Ca-plddt of all residues for each chain.
		self.plddt_dict = self.get_ca_plddt()
		# Average PAE matrix.
		data = self.get_data_dict()
		self.pae = self.get_pae( data )
		# Get minPAE for each residue.
		self.min_pae = self.get_min_pae( avg_pae = self.pae, 
										lengths_dict = self.lengths_dict, 
										mask_intrachain = True,
										return_dict = True )


class SaveConfidentPredictions( Initialize ):
	def __init__( self, struct_file_path: str, data_file_path: str ):
		super().__init__( struct_file_path, data_file_path )

		self.apply_plddt = True
		self.apply_pae = True
		self.plddt_cutoff = 70
		self.pae_cutoff = 5

	def save_confident_regions( self ):
		"""
		Select confident residues based on plddt and min_pae.
		Save the confident residues as a CIF file.
		"""
		confident_residues = {}

		if not apply_plddt and not apply_pae:
			raise Esception( "No confidence filter applied..." )

		for chain in self.res_dict:
			confident_residues[chain] = []
			for res in self.res_dict[chain]:
				select = False
				if self.apply_plddt;
					if self.plddt_dict[chain][res] >= self.plddt_cutoff:
						select = True

				if self.apply_pae;
					if self.min_pae_dict[chain][res] <= self.pae_cutoff:
						select = True

				if select:
					confident_residues[chain].append( res )

		ResidueSelect( confident_residues )


class Interaction( Initialize ):
	# Create a contact map.
	# Obtain confident interactions.
	# Return the required interacting residues.
	def __init__( self, struct_file_path: str, 
						data_file_path: str, 
						interacting_region: Dict ):
		super().__init__( struct_file_path, data_file_path )
		self.interacting_region = interacting_region

		self.contact_threshold = 8
		self.plddt_cutoff = 70
		self.pae_cutoff = 5


	def get_chains_n_indices( self ):
		"""
		Obtain the chain IDs and residues indices 
			for the required interacting region.
		residue_index = residue_position - 1
		"""
		chain1, chain2 = self.interacting_region.keys()
		mol1_res1, mol1_res2 = self.interacting_region[chain1]
		mol1_res1 -= 1
		mol2_res1, mol2_res2 = self.interacting_region[chain2]
		mol2_res1 -= 1

		return [chain1, chain2], [mol1_res1, mol1_res2], [mol2_res1, mol2_res2]


	def get_required_coords( self, chains, mol1_res, mol2_res ):
		"""
		Get the coordinates for the interacting region 
			for which confident interactions are required.
		"""
		chain1, chain2 = chains
		start1, end1 = mol1_res
		start2, end2 = mol2_res
		coords1 = self.coords_dict[chain1][start1:end1,:]
		coords2 = self.coords_dict[chain2][start2:end2,:]

		return coords1, coords2


	def get_required_plddt( self, chains, mol1_res, mol2_res ):
		"""
		Get the plddt for the interacting region 
			for which confident interactions are required.
		"""
		chain1, chain2 = chains
		start1, end1 = mol1_res
		start2, end2 = mol2_res
		plddt1 = self.plddt_dict[chain1][start1:end1]
		plddt2 = self.plddt_dict[chain2][start2:end2]

		return plddt1, plddt2


	def get_required_pae( self, chains, mol1_res, mol2_res ):
		"""
		Get the PAE matrix for the interacting region.
			For this we need the cumulative residue index 
				uptil the required residue position.
		"""
		chain1, chain2 = chains
		start1, end1 = mol1_res
		start2, end2 = mol2_res

		# Count total residues till start1 and start2.
		cum_start1, cum_start2 = 0, 0
		for chain in self.res_dict:
			if chain == chain1:
				cum_start1 += start1
				break
			else:
				cum_start1 += len( self.res_dict[chain] )

		for chain in self.res_dict:
			if chain == chain2:
				cum_start2 += start2
				break
			else:
				cum_start2 += len( self.res_dict[chain] )

		cum_end1 = cum_start1 + ( end1 - start1 )
		cum_end2 = cum_start2 + ( end2 - start2 )

		pae = self.pae[cum_start1:cum_end1, cum_start2:cum_end2]

		return pae


	def get_confident_interactions( self ):
		"""
		For the specified regions in the predicted structure, 
			obtain all confident interacting residue pairs.
		"""
		chains, mol1_res, mol2_res = self.get_chains_n_indices()

		coords1, coords2 = self.get_required_coords( chains, mol1_res, mol2_res )

		contact_map = get_contact_map( 
								get_distance_map( coords1, coords2 ),
								self.contact_threshold
								 )

		plddt1, plddt2 = self.get_required_plddt( chains, mol1_res, mol2_res )
		pae = self.get_required_pae( chains, mol1_res, mol2_res )

		plddt1 = np.where( plddt1 >= self.plddt_cutoff, 1, 0 )
		plddt2 = np.where( plddt2 >= self.plddt_cutoff, 1, 0 )
		plddt_matrix = plddt1 * plddt2.T

		pae = np.where( pae >= self.pae_cutoff, 1, 0 )
		print( contact_map.shape, "  ", plddt_matrix.shape, "  ", pae.shape )
		confident_interactions = contact_map * plddt_matrix * pae

		return confident_interactions

