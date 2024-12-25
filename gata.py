from main import Interaction

struct_file = "./example/fold_e_gata3_dna_binding_domain_motif_c_terminus_model_0.cif"
data_file = "./example/fold_e_gata3_dna_binding_domain_motif_c_terminus_full_data_0.json"

interactting_regions = {"A": [1, 10], "B": [1, 10]}

obj = Interaction( struct_file, data_file, interactting_regions )

print( obj.get_confident_interactions() )


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

