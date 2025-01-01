import numpy as np
import matplotlib.pyplot as plt
import argparse
from Bio import SeqIO
import json
import random
import os

from typing import Dict, List

from main import Interaction, SaveConfidentPredictions



class GataRaheMeraDil():
# class GataAnalysis():
	def __init__( self ):
		self.base_dir = "../GATA_NV_SRlab/"

		# FASTA files for the input proteins.
		self.wt_gata_fasta = f"{self.base_dir}gata3_wt.fasta"
		self.mut_gata_fasta = f"{self.base_dir}gata3_mut.fasta"
		self.dna_gata_fasta = f"{self.base_dir}gata3_dna.fasta"

		# Dir containing the AF3 predictions.
		self.wt_gata_dna_preds_dir = f"{self.base_dir}wt_gata_af3/"
		self.mut_gata_dna_preds_dir = f"{self.base_dir}mut_gata_af3/"

		# Set the seed for PRNG.
		np.random.seed( 1 )
		random.seed( 1 )

		# if not self.args.cj and not self.args.a:
		# 	raise Exception( "Either of the flags '-cj' or '-a' must be provided..." )


	def tu_hi_meri_manzil( self ):
	# def forward( self ):
		"""
		Determine whether there is a significant difference between the 
			inetractions formed by the wt_GATA and mut_GATA with DNA.
		"""
		self.kahin_beete_na_ye_raatein()
		self.kahin_beete_na_ye_din()
		self.arre_pyar_krne_waale()
		# self.get_confident_predictions()



	def kahin_beete_na_ye_raatein( self ):
		"""
		Parse the FASTA files for wt-GATA and mut-GATA, and GATA-DNA
			to obtain their sequences.
		"""
		self.lengths_dict = {}
		for record in SeqIO.parse( self.wt_gata_fasta, "fasta" ):
			self.wt_gata = str( record.seq )
			self.lengths_dict["wt_gata"] = len( self.wt_gata )

		for record in SeqIO.parse( self.mut_gata_fasta, "fasta" ):
			self.mut_gata = str( record.seq )
			self.lengths_dict["mut_gata"] = len( self.mut_gata )

		for record in SeqIO.parse( self.dna_gata_fasta, "fasta" ):
			self.dna_gata = str( record.seq )
			self.lengths_dict["dna_gata"] = len( self.dna_gata )



	def kahin_beete_na_ye_din( self ):
		"""
		Given the input sequence for wt-GATA and e-GATA,
			create input JSON file for AF3.
		Create a JSON file with 20 entries with different seeds.
		"""
		self.prot_ids = {}
		for prot, seq in zip( ["wt_gata", "mut_gata"], [self.wt_gata, self.mut_gata] ):
			self.prot_ids[f"{prot}"] = []
			
			af3_batch = []
			seeds = np.arange( 1, 4001, 200 )
			for seed in seeds:
				entry_id = f"{prot}_{seed}"
				af3_entry = {}
				af3_entry["name"] = entry_id
				af3_entry["modelSeeds"] = [float( seed )]
				af3_entry["sequences"] = [
									{
									"proteinChain": {
										"sequence": seq,
										"count": 2
									} },
									{
									"dnaSequence": {
										"sequence": self.dna_gata,
										"count": 2
									} }
				]

				af3_batch.append( af3_entry )

				self.prot_ids[f"{prot}"].append( entry_id )

			file = f"{prot}_af3_batch_{seeds[0]}-{seeds[-1]}"
			if not os.path.exists( f"{self.base_dir}{file}.json" ):
				with open( f"{self.base_dir}{id_}.json", "w" ) as w:
					json.dump( af3_batch, w )
			else:
				print( "AF3 input already created..." )



	def arre_pyar_krne_waale( self ):
	# def get_confident_predictions( self ):
		"""
		For both wt-GATA and mut-GATA, get confident predictions or load from disk if already existing.
		"""
		self.interactions = {}
		for prot in ["wt_gata"]:
			out_file = f"{self.base_dir}{prot}-dna.npy"
			if not os.path.exists( out_file ):
				gata_dna_interaction = self.jalne_waale_chahe( self.prot_ids[prot] )
				# gata_dna_interaction = self.extract_wt_mut_predictions( self.prot_ids[prot] )

				self.interactions[f"{prot}_dna"] = gata_dna_interaction
				np.save( out_file, gata_dna_interaction )
			else:
				self.interactions[f"{prot}_dna"] = np.load( out_file )


		p = self.interactions[f"{prot}_dna"]
		p = np.sum( p, axis = 0 )
		plt.imshow( np.where( p >0, 1, 0 ) )
		plt.show()



	def jalne_waale_chahe( self, entry_ids ):
	# def extract_wt_mut_predictions( self, entry_ids ):
		"""
		Get the confident interactions as a binary contact map.
		"""
		gata_dna = []

		# for files in self.yield_struct_data_files( entry_ids ):
		for files in self.jal_jal_marenge( entry_ids ):
			confident_interactions = self.milke_jo_dhadke_hain_do_dil( *files )
			# confident_interactions = self.get_confident_interactions( *files )
			gata_dna.append(
					self.concat_interactions( region = 0, confident_interactions = confident_interactions )
					)

		gata_dna = np.stack( gata_dna )
		return gata_dna



	def jal_jal_marenge( self, entry_ids: List ):
	# def yield_struct_data_files( self, entry_ids: List ):
		"""
		Yield the AF3 pred files for wt-GATA-DNA and mut-GATA-DNA.
		"""
		for entry in entry_ids:
			print( f"Entry: {entry}..." )
			for i in range( 5 ):
				struct_file = f"{self.wt_gata_dna_preds_dir}fold_{entry}/fold_{entry}_model_{i}.cif"
				data_file = f"{self.wt_gata_dna_preds_dir}fold_{entry}/fold_{entry}_full_data_{i}.json"

				yield struct_file, data_file



	def hardam_ye_kahenge( self, region: int ):
	# def get_interacting_region( self, region: int ):
		"""
		Get the regions of interest for both proteins for which 
			the confident interactions are required.
		Regions of interest:
			GATA( dna binding domains ) and DNA
		"""
		# GATA DNA binding domain (263-287 and 317-341).
		dbd = [263, 341]
		dna = [1, self.lengths_dict["dna_gata"]]

		for chain1 in ["A", "B"]:
			for chain2 in ["C", "D"]:
				interacting_region = {
							chain1: dbd,
							chain2: dna
								}
				yield interacting_region



	def milke_jo_dhadke_hain_do_dil( self, struct_file: str, data_file: str ):
	# def get_confident_interactions( self, struct_file: str, data_file: str ):
		"""
		Obtain confident interactions, given the structure and data file.
		"""
		# Initialize the Interactions class.
		obj = Interaction( struct_file, data_file )

		confident_interactions = []
		
		# for interacting_regions in self.get_interacting_region( 0 ):
		for interacting_regions in self.hardam_ye_kahenge( 0 ):

			pair_interactions = obj.get_confident_interactions( interacting_regions )
			confident_interactions.append( pair_interactions )

		return confident_interactions



	def concat_interactions( self, region, confident_interactions ):
		"""
		Based on the interacting regions, concat the pairwise confident predictions.
		"""
		p1_d1, p1_d2, p2_d1, p2_d2 = confident_interactions

		p1_d = np.concatenate( [p1_d1, p1_d2], axis = 1 )
		p2_d = np.concatenate( [p2_d1, p2_d2], axis = 1 )

		combined_map = np.concatenate( [p1_d, p2_d], axis = 0 )

		return combined_map



# struct_file = "./example/fold_e_gata3_dna_binding_domain_motif_c_terminus_model_0.cif"
# data_file = "./example/fold_e_gata3_dna_binding_domain_motif_c_terminus_full_data_0.json"

# interacting_regions = {"A": [1, 10], "B": [1, 10]}

# obj = Interaction( struct_file, data_file )

# print( obj.get_confident_interactions( interacting_regions ) )

# obj = SaveConfidentPredictions( struct_file, data_file, "trial.pdb" )
# obj.save_confident_regions()

"""
We need a statistical test for assessing the significance 
	in the difference between the interactions formed in 
	the WT and mutant GATA binding to DNA.

Assumptions for the test:
	1. Paired data.
	2. Should work for interval scaling.
	3. No assumtion of normal data distribution.
Wilcoxn signed-ranked test meets these criterion.
"""

if __name__ == "__main__":
	GataRaheMeraDil().tu_hi_meri_manzil()
	print( "May the Force be with you.." )
	# GataAnalysis().forward()
# 	parser = argparse.ArgumentParser(
# 						prog = "GATA analysis",
# 						description = "Perform analysis for wt-GATA and e-GATA interaction with DNA."
# 		)

# 	parser.add_argument( 
# 						"--wt_gata", "-wt", dest = "wt", 
# 						help = "wt-GATA FASTA file path.", 
# 						type = str, required = True )
# 	parser.add_argument( 
# 						"--e_gata", "-e", dest = "e", 
# 						help = "e-GATA FASTA file path.", 
# 						type = str, required = True )
# 	parser.add_argument( 
# 						"--dna", "-d", dest = "d", 
# 						help = "DNA sequence FASTA file path.", 
# 						type = str, required = True )
# 	parser.add_argument( 
# 						"--create_json", "-cj", dest = "cj", 
# 						help = "Create input JSON file for running AF3.", 
# 						action = "store_true", type = bool, required = False, default = False )
# 	parser.add_argument( 
# 						"--analysis", "-a", dest = "a", 
# 						help = "Perform analysis for GATA interaction with DNA.", 
# 						action = "store_true", type = bool, required = False, default = False )



