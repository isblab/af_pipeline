import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon, permutation_test
import argparse
from Bio import SeqIO
import json
import random
import os
import time

from typing import Dict, List

from main import Interaction, SaveConfidentPredictions


def make_prot_pair( self, p1, p2, sep ):
	return f"{p1}{sep}{p2}"


def get_prot_from_pair( self, prot_pair, sep ):
	return prot_pair.split( sep )


class ClsyAnalysis():
	def __init__( self ):
		self.base_dir = "../CLSY_R_SPlab/"

		# Defined in forward.
		self.contact_threshold = [8]
		self.plddt_cutoff = [70]
		self.pae_cutoff = [5]
		# Dict containing the sequences and lengths of all proteins of interest.
		self.prot_dict = {}

		# Single FASTA file containing all required protein sequences.
		self.fasta_path = f"{self.base_dir}Sequence_CLSYs_SHHs.fasta"

		# Dir containing the AF3 predictions.
		self.os_af3_preds_dir = f"{self.base_dir}Os_preds/"
		self.at_af3_preds_dir = f"{self.base_dir}At_preds/"

		# Dir to save predicted files.
		self.output_dir = f"{self.base_dir}output/"
		if not os.path.exists( self.output_dir ):
			os.makedirs( self.output_dir )

		# Set the seed for PRNG.
		np.random.seed( 1 )
		random.seed( 1 )


	def forward( self ):
		"""
		Determine whether there is a significant difference between the 
			inetractions formed by the wt_GATA and mut_GATA with DNA.
		"""
		self.initialize_prot_of_interest()
		self.read_input_fasta()
		self.get_combinations()
		self.create_af3_input()


	def initialize_prot_of_interest( self ):
		"""
		Get all the protein of interest.
		For some proteins an alternate identifiers is also specified.
		"""
		self.proteins_of_interest = {
					"Os": [
						["CLSY1", ""],
						["CLSY3", ""],
						["CLSY4", ""],
						["SHH1", ""],
						["SHH1-L1", "SHH1-like1"],
						["SHH2-L1", "SHH2like1"],
						["SHH2-L2", "SHH2like2"]
					],
					"At": [
						["CLSY1", ""],
						["CLSY2", ""],
						["CLSY3", ""],
						["CLSY4", ""],
						["SHH1", ""],
						["SHH2", ""]
					]
		}



	def read_input_fasta( self ):
		"""
		Parse the FASTA file containing CLSY and SHH sequences for 
			Orzya sativa (Os) and Arabidopsis thaliana (At).
		"""
		print( "Reading input FASTA file..." )
		self.prot_dict = {k:{} for k in ["Os", "At"]}

		for record in SeqIO.parse( self.fasta_path, "fasta" ):
			header = record.id
			# Neglesctng the '*' at he end, where present.
			seq = str( record.seq )
			seq = seq[:-1] if "*" in seq else seq
			length = len( seq )

			for org in self.proteins_of_interest.keys():
				for protein in self.proteins_of_interest[org]:
					prot, alt_prot = protein

					# If alternate identifier is present.
					if len( alt_prot ) != 0:
						if alt_prot in header:
							self.prot_dict[org][prot] = {
										"seq": seq,
										"len": length
							}

					else:
						if org in header and prot in header:
							self.prot_dict[org][prot] = {
											"seq": seq,
											"len": length
								}
		# for k in self.prot_dict:
		# 	print( k, " ------------" )
		# 	for k1 in self.prot_dict[k]:
		# 		print( k1 )
		# 		print( self.prot_dict[k][k1]["seq"] )
		# 		print( self.prot_dict[k][k1]["len"] )


	def combinations( self, prot1: List, prot2: List ):
		"""
		Create all-vs-all combinations of prot1 and prot2.
		"""
		combos = []
		for p1 in prot1:
			for p2 in prot2:
				combos.append( 
						make_prot_pair( p1, p2 )
						 )

		return combos



	def get_combinations( self ):
		"""
		Create all required protein pairs to be analyzed.
		"""
		self.prot_combos = {}
		
		for org in self.proteins_of_interest:
			
			clsy, shh = [], []
			for protein in self.proteins_of_interest[org]:
				prot, _ = protein
				if "CLSY" in prot:
					clsy.append( prot )
				elif "SHH" in prot:
					shh.append( prot )
				else:
					raise Exception( "Foreign protein encountered..." )

			self.prot_combos[org] = self.combinations( clsy, shh )
		# for k in self.prot_combos:
		# 	print( self.prot_combos[k] )
		# 	print( len( self.prot_combos[k] ) )



	def create_af3_input( self ):
		"""
		Given the input sequence for wt-GATA and e-GATA,
			create input JSON file for AF3.
		Create a JSON file with 20 entries with different seeds.
		"""
		print( "Creating AF3 input files..." )
		seed = 1.0
		
		for org in self.prot_combos:
			af3_batch = []

			for comb in self.prot_combos[org]:
				print( comb )
				entry_id = comb
				p1, p2 = get_prot_from_pair( comb )

				seq1 = self.prot_dict[org][p1]["seq"]
				seq2 = self.prot_dict[org][p2]["seq"]
					
				af3_entry = {}
				af3_entry["name"] = entry_id
				af3_entry["modelSeeds"] = [seed]
				af3_entry["sequences"] = [
									{
									"proteinChain": {
										"sequence": seq1,
										"count": 1
									} },
									{
									"proteinChain": {
										"sequence": seq2,
										"count": 1
									} }				]

				af3_batch.append( af3_entry )

			file = f"{org}_af3_batch_{seed}.json"
			if not os.path.exists( f"{self.base_dir}{file}" ):
				with open( f"{self.base_dir}{file}", "w" ) as w:
					json.dump( af3_batch, w )
			else:
				print( f"AF3 input already created for {org}..." )



	def run_analysis( self ):
		"""
		Obtain contact maps for the best model for all predictions.
		For both organisms Os and At:
			Get the paths to the structure and data files.
			Get the contact maps and confidence metrics: pLDDT and PAE.
			Get confident interactions for different values of:
				Contact threshold, pLDDT and PAE cutoff.
			Save the predictions for all prot pairs of an organism as a .npy file.
		"""
		for org in self.prot_combos:
			for pair in self.prot_combos[org]:



	def yield_struct_data_files( self, org: str, prot_pair: str ):
		"""
		Obtain structure and adat file paths for the 
			best model for the given entry.
		"""
		if org == "Os":
			path = self.os_af3_preds_dir
		else:
			path = self.at_af3_preds_dir

		struct_file = f"{path}fold_{prot_pair}/fold_{prot_pair}_model_{i}.cif"
		data_file = f"{path}fold_{prot_pair}/fold_{prot_pair}_full_data_{i}.json"

		return struct_file, data_file



	def get_interacting_region( self, org: str, prot_pair: str ):
		"""
		Create an dict containing regions for which interaction map is required.
		By construction, prot1 and prot2 will always correspond to chain A and B respectively.
		"""
		p, p2 = get_prot_from_pair( prot_pair )
		p1_len = self.prot_dict[org][p1]["len"]
		p2_len = self.prot_dict[org][p2]["len"]

		interacting_region = {"A": [1, p1_len], "B": [1, p2_len]}

		return interacting_region



	def set_thresholds( self, obj: Interaction, cthresh: int, pldd: int, pae: int ):
		"""
		Set the thresholds for the Interaction object.
		"""
		obj.contact_threshold = cthresh
		obj.plddt_cutoff = plddt
		obj.pae_cutoff = pae



	def get_confident_interactions( self, obj: Interaction, 
										contact_map: np.array, 
										plddt1: np.array, pldd2: np.array, 
										pae: np.array ):
		"""
		Get confident interactions across different metric cutoffs:
			Contact threshold
			pLDDT
			PAE
		"""
		contact_threshold = self.contact_threshold
		plddt_cutoff = self.plddt_cutoff
		pae_cutoff = self.pae_cutoff

		confident_interactions = {}

		for cthresh in contact_threshold:
			confident_interactions[cthresh] = {}
			for plddt in plddt_cutoff:
				confident_interactions[cthresh][plddt] = {}
				for pae_ in pae_cutoff:
					confident_interactions[cthresh][plddt][pae] = {}
					self.set_thresholds( obj, cthresh, pldd, pae_ )

					plddt_matrix, pae_matrix = obj.apply_confidence_cutoffs( plddt1, pldd2, pae )
					confident_contact_map = interaction_map * plddt_matrix * pae_matrix

					confident_interactions[cthresh][plddt][pae] = confident_contact_map

		return confident_interactions



	def get_contact_maps_from_struct( self, org: str, prot_pair: str ):
		"""
		Obtain contact maps for the best model for the given entry.
		"""
		struct_file, data_file = self.yield_struct_data_files( org, prot_pair )

		obj = Interaction( struct_file, data_file )

		interacting_region = self.get_interacting_region( org, prot_pair )

		contact_map, plddt1, pldd2, pae = obj.get_interaction_data( interacting_region )

		confident_interactions = self.get_confident_interactions( obj, contact_map, plddt1, pldd2, pae )

		return confident_interactions






if __name__ == "__main__":
	tic = time.time()
	ClsyAnalysis().forward()
	toc = time.time()
	t = toc - tic
	print( "\n------------------------------------------------\n" )
	print( f"Time taken: {t/60} minutes or {t} seconds" )
	print( "May the Force be with you.." )
