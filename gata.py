import numpy as np
import argparse
from Bio import SeqIO
import json
import random

from typing import Dict, List

from main import Interaction, SaveConfidentPredictions



class GataRaheMeraDil():
# class GataAnalysis():
	def __init__( self ):
		# self.args = args
		self.create_json = True
		# Identifiers for all input entries for AF3 for wt-GATA-DNA.
		self.wt_gata_ids = []
		# Identifiers for all input entries for AF3 for mut-GATA-DNA.
		self.mut_gata_ids = []
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
		if self.create_json:
			self.kahin_beete_na_ye_raatein()
			self.kahin_beete_na_ye_din()
		else:
			# self.arre_pyaar_krne_waale()
			pass


	def kahin_beete_na_ye_raatein( self ):
	# def parse_fasta( self ):
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
	# def create_json_input( self ):
		"""
		Given the input sequence for wt-GATA and e-GATA,
			create input JSON file for AF3.
		Create a JSON file with 20 entries with different seeds.
		"""
		for prot, seq in zip( ["wt_gata", "mut_gata"], [self.wt_gata, self.mut_gata] ):
			i = 0
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

				id_ = {prot}_af3_batch_{seeds[0]}-{seeds[-1]}
				if prot == "wt_gata":
					self.wt_gata_ids.append( id_ )
				else:
					self.mut_gata_ids.append( id_ )

			with open( f"{self.base_dir}{id_}.json", "w" ) as w:
				json.dump( af3_batch, w )


	def jalne_waale_chahe( self, entry_ids: List ):
	# def yield_struct_data_files( self, entry_ids ):
		"""
		Yield the AF3 pred files for wt-GATA-DNA and mut-GATA-DNA.
		"""
		for entry in entry_ids:
			for i in range( 5 ):
				struct_file = f"{self.wt_gata_dna_preds_dir}fold_{ntry}_model_{i}.cif"
				data_file = f"{self.wt_gata_dna_preds_dir}fold_{ntry}_full_data_{i}.json"

				yield struct_file, data_file


	def jalne_waale_chahe( self, region: int ):
	# def get_interacting_region( self, region: int ):
		"""
		Get the regions of interest for both proteins for which 
			the confident interactions are required.
		Regions of interest:
			0 --> GATA( dna binding domains ) and DNA
			1 --> GATA( dna binding domains ) and GATA( C-term )
			2 --> GATA( C-term ) and DNA
		"""
		# GATA DNA binding domain (263-287 and 317-341).
		dbd = [263, 341]
		dna = [1, self.lengths_dict["dna_gata"]]
		c_term_wt = [342, self.lengths_dict["wt_gata"]]
		c_term_mut = [342, self.lengths_dict["mut_gata"]]

		if region == 0:
			interactting_regions = {
						"A": dbd,
						"B": dbd,
						"C": dna,
						"D": dna
							}


	def arre_pyaar_hi_krenge( self, interactting_regions: Dict, struct_file: str, data_file: str ):
	# def get_confident_interactions( self ):
		"""
		Obtain confident interactions, given the structure and data file.
		"""

		obj = Interaction( struct_file, data_file )
		confident_interactions = obj.get_confident_interactions( interactting_regions )

		return confident_interactions


	def arre_pyaar_krne_waale( self ):
	# def get_wt_mut_predictions( self ):
		"""
		For both the wt-GATA and mut-GATA, get the confident interactions as a binary contact map.
		"""
		wt_gata_dna, mut_gata_dna = [], []

		wt_gata_dna.append( 
						self.pyaar_hi_krenge( 
											self.jalne_waale_chahe( 0 ),
											self.jal_jal_marenge( self.wt_gata_ids )
										 )
		 )

		wt_gata_dna = np.stack( wt_gata_dna )
		# self.get_confident_interactions( 
		# 								self.yield_struct_data_files( self.wt_gata_ids )
		# 								 )




# struct_file = "./example/fold_e_gata3_dna_binding_domain_motif_c_terminus_model_0.cif"
# data_file = "./example/fold_e_gata3_dna_binding_domain_motif_c_terminus_full_data_0.json"

# interactting_regions = {"A": [1, 10], "B": [1, 10]}

# obj = Interaction( struct_file, data_file )

# print( obj.get_confident_interactions( interactting_regions ) )

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



