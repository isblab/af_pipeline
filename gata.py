from main import Interaction

struct_file = "./example/fold_e_gata3_dna_binding_domain_motif_c_terminus_model_0.cif"
data_file = "./example/fold_e_gata3_dna_binding_domain_motif_c_terminus_full_data_0.json"

interactting_regions = {"A": [1, 10]}

Interaction( struct_file, data_file, interactting_regions )
