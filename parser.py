import numpy as np
import os
import json
import pickle as pkl
import Bio
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select

from typing import Dict



class ResidueSelect( Select ):
	def __init__( self, confident_residues ):
		self.confident_residues = confident_residues

	def accept_residue( self, residue ):
		chain = residue.parent.id
		return residue.id[1] in self.confident_residues[chain]


class AfParser():
	def __init__( self, struct_file_path: str, data_file_path: str ):
		# AF2/3 structure file path.
		self.struct_file_path = struct_file_path
		# AF2/3 structure data file path.
		self.data_file_path = data_file_path
		
		# Biopython Structure object.
		self.structure = self.get_structure( 
									self.get_parser()
									 )


	def get_parser( self ):
		"""
		Get the required parser (PDB/CIF) for the input file.
		"""
		ext = os.path.splitext( self.struct_file_path )[1]

		if "pdb" in ext:
			parser = PDBParser()
		elif "cif" in ext:
			parser = MMCIFParser()
		else:
			raise Exception( "Incorrect file format.. Suported .pdb/.cif only." )

		return parser


	def get_structure( self, parser: Bio.PDB.PDBParser ):
		"""
		Return the Biopython Structure object for the input file.
		"""
		basename = os.path.basename( self.struct_file_path )
		structure = parser.get_structure( basename, self.struct_file_path )

		return structure


	def get_residues( self ):
		"""
		Get all residues in the structure.
		"""
		coords = []
		for model in self.structure:
			for chain in model:
				chain_id = chain.id[0]
				for residue in chain:
					yield residue, chain_id


	def save_pdb( self, res_select_obj: Bio.PDB.Select, out_file: str ):
		"""
		Given the ResidueSelect object, save the structure as a PDB file.
		"""
		io = PDBIO()
		io.set_structure( self.structure )
		io.save( out_file, res_select_obj )


	def extract_perresidue_quantity( self, residue, quantity ): 
		"""
		Given the Biopython residue object, return the specified quantity:
			1. residue position
			2. Ca-coordinate
			3. Ca-pLDDT
		"""
		symbol = residue.get_resname()
		# Using representative atoms as specified by AF3.
		if symbol in ["DA", "DG", "DC", "DT", "GLY"]:
			if symbol == "GLY":
				# Use Ca-atom for Glycine.
				rep_atom = "CA"
			elif symbol in ["DA", "DG"]:
				# C4 for purines.
				rep_atom = "C4"
			else:
				# C2 for pyrimidines.
				rep_atom = "C2"
		else:
			# Use Cb-atom for all other amino acids.
			rep_atom = "CB"

		if quantity == "res_pos":
			return residue.id[1]

		elif quantity == "coords":
			coords = residue[rep_atom].coord
			return coords
		
		elif quantity == "plddt":
			plddt = residue[rep_atom].bfactor
			return plddt
		
		else:
			raise Exception( f"Specified quantity: {quantity} does not exist..." )


	def get_residue_positions( self ):
		"""
		Get the residue positions for all residues.
		"""
		res_dict = {}
		for residue, chain_id in self.get_residues():
			res_id = self.extract_perresidue_quantity( residue, "res_pos" )
			if chain_id not in res_dict.keys():
				res_dict[chain_id] = np.array( res_id )
			else:
				res_dict[chain_id] = np.append( res_dict[chain_id], res_id )

		res_dict = {k: v.reshape( -1, 1 ) for k, v in res_dict.items()}

		return res_dict


	def get_chain_lengths( self, res_dict: Dict ):
		"""
		Create a dict containing the length of all chains in the system 
			and the total length of the system.
		"""
		lengths_dict = {}
		lengths_dict["total"] = 0
		for chain in res_dict:
			chain_length = len( res_dict[chain] )
			lengths_dict[chain] = chain_length
			lengths_dict["total"] += chain_length

		return lengths_dict


	def get_ca_coordinates( self ):
		"""
		Get the coordinates for all Ca atoms of all residues.
		"""
		coords_dict = {}
		for residue, chain_id in self.get_residues():
			coords = self.extract_perresidue_quantity( residue, "coords" )
			if chain_id not in coords_dict.keys():
				coords_dict[chain_id] = np.array( coords )
			else:
				coords_dict[chain_id] = np.append( coords_dict[chain_id], coords )

		coords_dict = {k: v.reshape( -1, 3 ) for k, v in coords_dict.items()}

		return coords_dict


	def get_ca_plddt( self ):
		"""
		Get the pLDDT score for all Ca atoms of all residues.
		"""
		plddt_dict = {}
		for residue, chain_id in self.get_residues():
			plddt = self.extract_perresidue_quantity( residue, "plddt" )
			if chain_id not in plddt_dict.keys():
				plddt_dict[chain_id] = np.array( [plddt] )
			else:
				plddt_dict[chain_id] = np.append( plddt_dict[chain_id], plddt )

		plddt_dict = {k: v.reshape( -1, 1 ) for k, v in plddt_dict.items()}

		return plddt_dict


	def get_data_dict( self ):
		"""
		Parse the AF2/3 data file.
			AF2 data file is saved as a .pkl file 
				whereas for AF3 it's stored as .json.
		"""
		ext = os.path.splitext( self.data_file_path )[1]

		if "pkl" in ext:
			with open( self.data_file_path, "rb" ) as f:
				data = pkl.load( f )

		elif "json" in ext:
			with open( self.data_file_path, "r" ) as f:
				data = json.load( f )
		else:
			raise Exception( "Incorrect file format.. Suported .pkl/.json only." )

		return data


	def get_pae( self, data: Dict ):
		"""
		Return the PAE matrix from the data dict.
			AF2/3 PAE matrix is asymmetric.
			Hence, we consider the average PAE: ( PAE + PAE.T )/2.
		"""
		# For AF2.
		if "predicted_aligned_error" in data.keys():
			pae = np.array( data["predicted_aligned_error"] )
		# For AF3.
		elif "pae" in data.keys():
			pae = np.array( data["pae"] )
		else:
			raise Exception( "PAE matrix not found..." )

		avg_pae = ( pae + pae.T )/2

		return avg_pae


	def create_interchain_mask( self, lengths_dict: Dict ):
		"""
		Create a binary 2D mask for selecting only interchain interactions.
		"""
		sys_len = lengths_dict.pop( "total" )
		interchain_mask = np.ones( ( sys_len, sys_len ) )

		prev = 0
		for chain in lengths_dict:
			l = lengths_dict[chain]
			curr = prev + l
			interchain_mask[prev:curr:, prev:curr] = 100
			prev += l

		return interchain_mask


	def get_min_pae( self, avg_pae: np.array, lengths_dict: Dict, 
						mask_intrachain: bool, return_dict: bool ):
		"""
		Given the averaged PAE matrix, obtain min PAe values for all residues.
			Esentially return a vector containing row-wise min PAE values.
		min_pae indicates whether the minimum error in the interaction with some residue.

		If mask_intrachain, only interchain interactions are considered.
		"""
		if mask_intrachain:
			interchain_mask = self.create_interchain_mask( lengths_dict )
			avg_pae = avg_pae * interchain_mask
		min_pae = np.min( avg_pae, axis = 1 )

		# Convert to a dict.
		min_pae_dict = {}
		start = 0
		for chain in lengths_dict:
			if chain != "total":
				end = start + lengths_dict[chain]
				min_pae_dict[chain] = min_pae[start:end]

				start = end

		if return_dict:
			return min_pae_dict
		else:
			return min_pae

