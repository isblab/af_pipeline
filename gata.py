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


class GataRaheMeraDil():
# class GataAnalysis():
	def __init__( self ):
		self.base_dir = "../GATA_NV_SRlab/"

		# Defined in forward.
		self.contact_threshold = [8]
		self.plddt_cutoff = [70]
		self.pae_cutoff = [5]
		self.significance_test = "permut_test"

		# FASTA files for the input proteins.
		self.wt_gata_fasta = f"{self.base_dir}gata3_wt.fasta"
		self.mut_gata_fasta = f"{self.base_dir}gata3_mut.fasta"
		self.dna_gata_fasta = f"{self.base_dir}gata3_dna.fasta"

		# Dir containing the AF3 predictions.
		self.wt_gata_dna_preds_dir = f"{self.base_dir}wt_gata_af3/"
		self.mut_gata_dna_preds_dir = f"{self.base_dir}mut_gata_af3/"

		# Dir to save predicted files.
		self.output_dir = f"{self.base_dir}output/"
		if not os.path.exists( self.output_dir ):
			os.makedirs( self.output_dir )
		# Output CSV file.
		self.output_file = f"{self.output_dir}result.csv"

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

		self.output = {k:[] for k in ["contact_dist", "pLDDT_cutoff", "PAE_cutoff", 
										"contacts_wt_GATA-DNA", "contacts_mut_GATA-DNA", 
										"statistic", "p-value"]}

		contact_threshold = self.contact_threshold
		plddt_cutoff = self.plddt_cutoff
		pae_cutoff = self.pae_cutoff
		
		for cthresh in contact_threshold:
			self.contact_threshold = cthresh
			for plddt in plddt_cutoff:
				self.plddt_cutoff = plddt
				for pae in pae_cutoff:
					print( f"Running: Contaact threshold: {cthresh} \t pLDDT: {plddt} \t PAE: {pae}" )
					self.pae_cutoff = pae

					self.output["contact_dist"].append( cthresh )
					self.output["pLDDT_cutoff"].append( plddt )
					self.output["PAE_cutoff"].append( pae )
					
					stat, pval = self.arre_pyaar_krne_waale()
					# stat, pval = self.run_analysis()
					self.output["contacts_wt_GATA-DNA"].append( np.count_nonzero( self.interactions["wt_gata_dna"] ) )
					self.output["contacts_mut_GATA-DNA"].append( np.count_nonzero( self.interactions["mut_gata_dna"] ) )
					self.output["statistic"].append( stat )
					self.output["p-value"].append( pval )

		self.write_output()


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


	def arre_pyaar_krne_waale( self ):
	# def run_analysis( self ):
		"""
		Parse all the predictions.
			Obtain all interacting pairs at the specified confience cutoffs.
		Aggregate the contacts.
		Assess statistical significance.
		Write output in a CSV file.
		"""
		suffix = f"{self.contact_threshold}_{self.plddt_cutoff}_{self.pae_cutoff}"
		# Files storing predicted interactions for wt-GATA and mut-GATA with DNA.
		self.wt_gata_dna_interactions = f"{self.output_dir}wt_gata-dna_{suffix}.npy"
		self.mut_gata_dna_interactions = f"{self.output_dir}mut_gata-dna_{suffix}.npy"

		# Plot the aggregated contact maps for wt-GATA-DNA and mut-GATA-DNA.
		self.interaction_map = f"{self.output_dir}interaction_map_{suffix}.png"

		self.pyaar_hi_krenge()
		# self.get_confident_predictions()
		stat, pval = self.test_statistical_significance()
		# if self.significance_test != "wilcox":
		# 	self.aggregate_predictions()

		self.plot_contact_map()

		return stat, pval



	def pyaar_hi_krenge( self ):
	# def get_confident_predictions( self ):
		"""
		For both wt-GATA and mut-GATA, get confident predictions or load from disk if already existing.
		"""
		self.interactions = {}
		for prot in ["wt_gata", "mut_gata"]:
			out_file = self.wt_gata_dna_interactions if prot == "wt_gata" else self.mut_gata_dna_interactions
			
			if not os.path.exists( out_file ):
				gata_dna_interaction = self.jalne_waale_chahe( prot, self.prot_ids[prot] )
				# gata_dna_interaction = self.extract_wt_mut_predictions( self.prot_ids[prot] )

				self.interactions[f"{prot}_dna"] = gata_dna_interaction
				np.save( out_file, gata_dna_interaction )
			else:
				self.interactions[f"{prot}_dna"] = np.load( out_file )



	def jalne_waale_chahe( self, prot: str, entry_ids: List ):
	# def extract_wt_mut_predictions( self, entry_ids ):
		"""
		Get the confident interactions as a binary contact map.
		"""
		gata_dna = []

		# for files in self.yield_struct_data_files( entry_ids ):
		for files in self.jal_jal_marenge( prot, entry_ids ):
			confident_interactions = self.milke_jo_dhadke_hain_do_dil( *files )
			gata_dna.append(
					self.kahin_beete_na( region = 0, confident_interactions = confident_interactions )
					)
					# kahin beete na ye raatein. Kahin beete na ye din.
			# confident_interactions = self.get_confident_interactions( *files )
			# gata_dna.append(
			# 		self.concat_interactions( region = 0, confident_interactions = confident_interactions )
			# 		)

		gata_dna = np.stack( gata_dna )
		return gata_dna



	def jal_jal_marenge( self, prot: str, entry_ids: List ):
	# def yield_struct_data_files( self, entry_ids: List ):
		"""
		Yield the AF3 pred files for wt-GATA-DNA and mut-GATA-DNA.
		"""
		for entry in entry_ids:
			print( f"Entry: {entry}..." )
			if prot == "wt_gata":
				path = self.wt_gata_dna_preds_dir
			else:
				path = self.mut_gata_dna_preds_dir
			
			for i in range( 5 ):
				struct_file = f"{path}fold_{entry}/fold_{entry}_model_{i}.cif"
				data_file = f"{path}fold_{entry}/fold_{entry}_full_data_{i}.json"

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
		# Adjust the cutofs as required.
		obj.contact_threshold = self.contact_threshold
		obj.plddt_cutoff = self.plddt_cutoff
		obj.pae_cutoff = self.pae_cutoff

		confident_interactions = []
		
		# for interacting_regions in self.get_interacting_region( 0 ):
		for interacting_regions in self.hardam_ye_kahenge( 0 ):

			pair_interactions = obj.get_confident_interactions( interacting_regions )
			confident_interactions.append( pair_interactions )

		return confident_interactions



	def kahin_beete_na( self, region, confident_interactions ):
	# def concat_interactions( self, region, confident_interactions ):
		"""
		Based on the interacting regions, concat the pairwise confident predictions.
		"""
		p1_d1, p1_d2, p2_d1, p2_d2 = confident_interactions

		p1_d = np.concatenate( [p1_d1, p1_d2], axis = 1 )
		p2_d = np.concatenate( [p2_d1, p2_d2], axis = 1 )

		combined_map = np.concatenate( [p1_d, p2_d], axis = 0 )

		return combined_map


	def aggregate_predictions( self ):
		"""
		Given contact maps across all predictions, create a summed contact map.
		"""
		for prot in self.interactions:
			self.interactions[prot] = np.sum( self.interactions[prot], axis = 0 )


	def plot_contact_map( self ):
		"""
		Plot binary contact maps for the wt-GATA and mut-GATA complex with DNA.
		"""
		fig, ax = plt.subplots( 2, 1, figsize = ( 10, 10 ) )

		p = self.interactions[f"wt_gata_dna"]
		ax[0].imshow( np.where( p > 0, 1, 0 ).T )
		ax[0].set_title( "wt_GATA-dbd with DNA" )
		ax[0].set_ylabel( "DNA" )
		ax[0].set_xlabel( "wt_GATA-dbd" )
		ax[0].set_xticklabels( np.arange( 263, 341 ) )

		p = self.interactions[f"mut_gata_dna"]
		ax[1].imshow( np.where( p > 0, 1, 0 ).T )
		ax[1].set_title( "mut_GATA-dbd with DNA" )
		ax[0].set_ylabel( "DNA" )
		ax[1].set_xlabel( "mut_GATA-dbd" )
		ax[1].set_xticklabels( np.arange( 263, 341 ) )
		plt.savefig( self.interaction_map, dpi = 300 )
		plt.close()



	def test_statistical_significance( self ):
		"""
		We need a statistical test for assessing the significance 
			in the difference between the interactions formed in 
			the WT and mutant GATA binding to DNA.

		Assumptions for the test:
			1. Paired data.
			2. Should work for interval scaling.
			3. No assumtion of normal data distribution.

		Null: there is no significant difference between the 
			confidently predicted interactions for the wt-GATA and mut-GATA with DNA.
		
		1. Wilcoxn signed-ranked test
			wilcox --> Discarding all zero-differences.
			two-sided --> the distribution underlying the difference is not symmetric about zero.
			However, this test is sensitive to sparsity.
		2. Permutation test
		"""
		self.aggregate_predictions()
		x = self.interactions["wt_gata_dna"]
		x = x.reshape( -1 )
		y = self.interactions["mut_gata_dna"]
		y = y.reshape( -1 )
		print( x.shape )
		def statistic( x, y ):
			return np.mean( x - y )

		result = permutation_test( ( x, y ), statistic, permutation_type = "samples", n_resamples = 5000 )

		stat = result.statistic
		pval = result.pvalue

		if self.significance_test == "wilcox":

			x = self.interactions["wt_gata_dna"]
			x = x.reshape( -1 )
			print( np.count_nonzero( x ) )
			y = self.interactions["mut_gata_dna"]
			y = y.reshape( -1 )
			print( np.count_nonzero( y ) )

			result = wilcoxon( x, y, zero_method = "wilcox", alternative = "two-sided" )
			stat = result.statistic
			pval = result.pvalue

			print( f"Statistic: {stat} \t p-value: {pval}" )
		return stat, pval


	def write_output( self ):
		"""
		Write the output to a CSV file.
		"""
		s = ""
		for k in self.output:
			s += f"{k}: {self.output[k]}"
		print( s )
		df = pd.DataFrame( self.output )
		df.to_csv( self.output_file )



if __name__ == "__main__":
	tic = time.time()
	GataRaheMeraDil().tu_hi_meri_manzil()
	toc = time.time()
	t = toc - tic
	print( "\n------------------------------------------------\n" )
	print( f"Time taken: {t/60} minutes or {t} seconds" )
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



