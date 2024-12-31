import numpy as np
import argparse
from Bio import SeqIO
import json
import random

from main import Interaction, SaveConfidentPredictions



class GataRaheMeraDil():
	def __init__( self, args ):
		self.args = args

		# Set the seed for PRNG.
		np.random.seed( 1 )
		random.seed( 1 )

		if not self.args.cj and not self.args.a:
			raise Exception( "Either of the flags '-cj' or '-a' must be provided..." )


	def tu_hi_meri_manzil( self ):
		"""
		Given the input sequence for wt-GATA and e-GATA,
			create input JSON file for AF3.
		Create a JSON file with 20 entries with different seeds.
		"""
		for prot, seq in zip( ["wt_gata", "e_gata"], [self.wt_gata, self.e_gata] ):
			af3_batch = []
			seeds = np.arange( 1, 401, 20 )
			for seed in seeds:
				name = f"{prot}_{seed}"
				af3_entry = {}
				af3_entry["name"] = entry_id
				af3_entry["modelSeeds"] = []
				af3_entry["sequences"] = [
									{
									"proteinChain": {
										"sequence": self.wt_gata,
										"count": 2
									},
									"dnaSequence": {
										"sequence": self.e_gata,
										"count": 2
									}
									}
				]

				af3_batch.append( af3_entry )

			with open( f"{prot}_af3_batch.json", "w" ) as w:
				json.dump( af3_batch )


	def kahin_beete_na_ye_raatein( self ):
		"""
		Parse the FASTA files for wt-GATA and e-GATA
			to obtain the sequences for both.
		"""
		for record in SeqIO.parse( self.args.wt, "fasta" ):
			self.wt_gata = record.seq

		for record in SeqIO.parse( self.args.e, "fasta" ):
			self.e_gata = record.seq


struct_file = "./example/fold_e_gata3_dna_binding_domain_motif_c_terminus_model_0.cif"
data_file = "./example/fold_e_gata3_dna_binding_domain_motif_c_terminus_full_data_0.json"

interactting_regions = {"A": [1, 10], "B": [1, 10]}

obj = Interaction( struct_file, data_file, interactting_regions )

print( obj.get_confident_interactions() )

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

# if __name__ == "__main__":
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



