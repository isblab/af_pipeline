import os
from typing import Dict
import warnings
import numpy as np
import pandas as pd
from af_pipeline._initialize import _Initialize

class Interaction(_Initialize):
    """ Class to handle interaction data for the predicted structure. \n
    One can obtain:
    1. Interaction map: A binary contact map or distance map.
    2. Restraints: Contacts as pairwise residues in a DataFrame format.
    3. Interacting patches: contiguous regions in the interaction map obtained in (1).
    """

    def __init__(
        self,
        struct_file_path: str,
        data_file_path: str,
        af_offset: dict | None = None,
        output_dir: str = "./output/af_output",
        idr_chains: list = [],
        **kwargs,
    ):

        super().__init__(
            struct_file_path=struct_file_path,
            data_file_path=data_file_path,
            af_offset=af_offset,
            **kwargs,
        )

        dir_name = os.path.basename(struct_file_path).split(".")[0]
        output_dir = os.path.join(output_dir, f"{dir_name}_patches")

        self.interaction_map_type = "contact"  # Either contact/distance.
        self.contact_threshold = 8  # Distance threshold in (Angstorm) to define a contact between residue pairs.
        self.plddt_cutoff = 70  # pLDDT cutoff to consider a confident prediction.
        self.idr_plddt_cutoff = 50  # pLDDT cutoff for IDR chains.
        self.pae_cutoff = 5 # PAE cutoff to consider a confident prediction.
        self.output_dir = output_dir
        self.idr_chains = idr_chains # List of chains that are disordered

        self.save_plot = False
        self.save_table = False

    def create_regions_of_interest(self):
        """
        Create regions of interest for all possible chain pairs.

        Returns:
            regions_of_interest (list): list of regions of interest
        """

        regions_of_interest = []
        token_chain_ids = self.token_chain_ids
        chain_pairs = set()

        for chain1 in set(token_chain_ids):
            for chain2 in set(token_chain_ids):
                if chain1 != chain2:
                    pair = tuple(sorted([chain1, chain2]))
                    chain_pairs.add(pair)

        chain_pairs = list(chain_pairs)

        for chain1, chain2 in chain_pairs:

            ch1_start = self.renumber.renumber_chain_res_num(
                chain_res_num=1,
                chain_id=chain1
            )
            ch1_end = self.renumber.renumber_chain_res_num(
                chain_res_num=self.lengths_dict[chain1],
                chain_id=chain1
            )
            ch2_start = self.renumber.renumber_chain_res_num(
                chain_res_num=1,
                chain_id=chain2
            )
            ch2_end = self.renumber.renumber_chain_res_num(
                chain_res_num=self.lengths_dict[chain2],
                chain_id=chain2
            )

            region_of_interest = {
                chain1: [ch1_start, ch1_end],
                chain2: [ch2_start, ch2_end],
            }

            regions_of_interest.append(region_of_interest)

        return regions_of_interest
