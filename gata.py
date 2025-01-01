import numpy as np
import argparse
from Bio import SeqIO
import json
import random

from main import Interaction, SaveConfidentPredictions



class GataRaheMeraDil():
# class GataAnalysis():
	def __init__( self ):
		# self.args = args
		self.create_json = True
		self.base_dir = "../GATA_NV_SRlab/"
		self.wt_gata_fasta = f"{self.base_dir}gata3_wt.fasta"
		self.mut_gata_fasta = f"{self.base_dir}gata3_mut.fasta"
		self.dna_gata_fasta = f"{self.base_dir}gata3_dna.fasta"

		# Set the seed for PRNG.
		np.random.seed( 1 )
		random.seed( 1 )

		# if not self.args.cj and not self.args.a:
		# 	raise Exception( "Either of the flags '-cj' or '-a' must be provided..." )


	def tu_hi_meri_manzil( self ):
	# def forward( self ):
		"""
		"""
		if self.create_json:
			self.kahin_beete_na_ye_raatein()
			self.kahin_beete_na_ye_din()
		else:
			pass


	def kahin_beete_na_ye_raatein( self ):
	# def parse_fasta( self ):
		"""
		Parse the FASTA files for wt-GATA and mut-GATA, and GATA-DNA
			to obtain their sequences.
		"""
		for record in SeqIO.parse( self.wt_gata_fasta, "fasta" ):
			self.wt_gata = str( record.seq )

		for record in SeqIO.parse( self.mut_gata_fasta, "fasta" ):
			self.mut_gata = str( record.seq )

		for record in SeqIO.parse( self.dna_gata_fasta, "fasta" ):
			self.dna_gata = str( record.seq )


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

			with open( f"{self.base_dir}{prot}_af3_batch_{seeds[0]}-{seeds[-1]}.json", "w" ) as w:
				json.dump( af3_batch, w )



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



