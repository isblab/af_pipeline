import numpy as np
import os
import json
import pickle as pkl
import Bio
from Bio.PDB import PDBParser, MMCIFParser

from typing import Dict

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


	def get_residue_positions( self ):
		"""
		Get the residue positions for all residues.
		"""
		res_dict = {}
		for residue, chain_id in self.get_residues():
			if chain_id not in res_dict.keys():
				res_dict[chain_id] = np.array( [residue.id[1]] )
			else:
				res_dict[chain_id] = np.append( res_dict[chain_id], residue.id[1] )

		res_dict = {k: v.reshape( -1, 1 ) for k, v in res_dict.items()}

		return res_dict


	def extract_quantity( self, residue, quantity ):
		"""
		Given the Biopython residue object, return the specified quantity:
			1. Ca-coordinate
			2. Ca-pLDDT
		"""
		symbol = residue.get_resname()
		if symbol in ["DA", "DG", "DC", "DT"]:
			atom = "P"
		else:
			atom = "CA"

		if quantity == "coords":
			coords = residue[atom].coord
			return coords
		
		elif quantity == "plddt":
			plddt = residue[atom].bfactor
			return plddt
		
		else:
			raise Exception( f"Specified quantity: {quantity} does not exist..." )


	def get_ca_coordinates( self ):
		"""
		Get the coordinates for all Ca atoms of all residues.
		"""
		coords_dict = {}
		for residue, chain_id in self.get_residues():
			coords = self.extract_quantity( residue, "coords" )
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
			plddt = self.extract_quantity( residue, "plddt" )
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



