import numpy as np
import pandas as pd
from itertools import combinations, product
from collections import defaultdict
from tqdm import tqdm
import warnings

class RigidBodyAssessment:

    def __init__(
        self,
        rb_dict: dict,
        num_to_idx: dict,
        idx_to_num: dict,
        contact_map: np.ndarray,
        plddt_list: np.ndarray,
        pae: np.ndarray,
        lengths_dict: dict,
        save_path: str,
        **kwargs,
    ):
        """ Initialize the RigidBodyAssessment class.

        Args:
            rb_dict (dict): Dictionary of rigid bodies, where each rigid body is a dictionary with chain IDs as keys and residue numbers as values.
            num_to_idx (dict): Residue number to index mapping.
            idx_to_num (dict): Index to residue number mapping.
            contact_map (np.ndarray): Binary contact map of the structure.
            plddt_list (np.ndarray): pLDDT scores for each residue in the structure.
            pae (np.ndarray): Predicted Alignment Error (PAE) matrix.
            lengths_dict (dict): Dictionary containing the lengths of each chain in the structure.
            save_path (str): Path to save the assessment results.
        """

        self.rb_dict = rb_dict
        self.num_to_idx = num_to_idx
        self.idx_to_num = idx_to_num
        self.symmetric_pae = kwargs.get("symmetric_pae", True)
        self.as_average = kwargs.get("as_average", False)
        self.idr_chains = kwargs.get("idr_chains", [])
        self.protein_chain_map = kwargs.get("protein_chain_map", {})

        self.unique_chains = self.get_unique_chains()
        self.chain_pairs = self.get_chain_pairs()
        self.rb_res_binary_map = self.get_rb_res_binary_map(lengths_dict=lengths_dict)
        self.rb_res_pairs = self.get_rb_res_pairs()

        self.per_chain_plddt = self.get_per_chain_plddt(plddt_list=plddt_list)
        self.per_chain_avg_plddt = self.get_per_chain_avg_plddt()
        self.pairwise_pae = self.get_pairwise_pae(pae=pae)

        self.interface_res_pairs = self.get_interface_res_pairs(contact_map=contact_map)
        self.per_chain_interface_residues = self.get_per_chain_interface_residues()
        self.num_contacts = self.get_num_contacts()
        self.num_interface_residues = self.get_num_interface_residues()

        self.pairwise_ipae = self.get_pairwise_ipae(pae=pae)

        self.per_chain_iplddt = self.get_per_chain_iplddt(plddt_list=plddt_list)
        self.per_chain_avg_iplddt = self.get_per_chain_average_iplddt()

        self.pairwise_avg_iplddt = self.get_pairwise_avg_iplddt()

        self.pairwise_min_pae = self.get_pairwise_min_pae()
        self.pairwise_avg_pae = self.get_pairwise_avg_pae(symmetric_pae=self.symmetric_pae)
        self.pairwise_avg_ipae = self.get_pairwise_avg_ipae(symmetric_pae=self.symmetric_pae)

        self.overall_assessment = self.get_overall_assessment()

        self.save_path = save_path

    # @time_it
    def save_rb_assessment(self):
        """ Save the assessment of the rigid bodies to an Excel file.

        The assessment includes:
        - Per chain assessment: Average pLDDT, Average iLDDT, Number of interface residues, Chain type (IDR or R)
        - Per chain pair assessment: Number of interface residues, Number of contacts, Average PAE, Average iPAE, Minimum PAE, Average iLDDT for each chain, Chain type (IDR or R) for each chain
        - Overall assessment: Average pLDDT, Average iLDDT, Number of interface residues, Chain type (IDR or R)

        The assessment is saved in an Excel file with three sheets:
        - "Chain Wise Assessment": Contains per chain assessment data.
        - "Chain Pairwise Assessment": Contains per chain pair assessment data.
        - "Overall Assessment": Contains overall assessment data.
        """

        chain_wise_assessment_rows = []
        chain_pairwise_assessment_rows = []
        overall_assessment_rows = []

        for chain_id in self.unique_chains:
            chain_wise_assessment_rows.append({
                "Chain ID": chain_id,
                "Protein Name": self.protein_chain_map.get(chain_id, None),
                "Average pLDDT": self.per_chain_avg_plddt[chain_id],
                "Average ipLDDT": self.per_chain_avg_iplddt.get(chain_id, np.nan),
                "Number of Interface Residues": len(self.per_chain_interface_residues[chain_id]),
                "Chain Type": "IDR" if chain_id in self.idr_chains else "R",
            })

        for chain_pair in self.chain_pairs:
            chain1, chain2 = chain_pair

            if self.as_average:

                if self.symmetric_pae:
                    chain_pairwise_assessment_rows.append({
                        "Chain Pair": f"{chain1}-{chain2}",
                        "Protein Name 1": self.protein_chain_map.get(chain1, "Unknown"),
                        "Protein Name 2": self.protein_chain_map.get(chain2, "Unknown"),
                        "Number of Interface Residues": self.num_interface_residues[chain_pair],
                        "Number of Contacts": self.num_contacts[chain_pair],
                        "Average PAE": self.pairwise_avg_pae[chain_pair],
                        "Average iPAE": self.pairwise_avg_ipae[chain_pair],
                        "Minimum PAE": self.pairwise_min_pae[chain_pair],
                        "Average ipLDDT chain1": self.pairwise_avg_iplddt[chain_pair].get(chain1, np.nan),
                        "Average ipLDDT chain2": self.pairwise_avg_iplddt[chain_pair].get(chain2, np.nan),
                        "Chain Type 1": "IDR" if chain1 in self.idr_chains else "R",
                        "Chain Type 2": "IDR" if chain2 in self.idr_chains else "R",
                    })
                else:
                    chain_pairwise_assessment_rows.append({
                        "Chain Pair": f"{chain1}-{chain2}",
                        "Protein Name 1": self.protein_chain_map.get(chain1, "Unknown"),
                        "Protein Name 2": self.protein_chain_map.get(chain2, "Unknown"),
                        "Number of Interface Residues": self.num_interface_residues[chain_pair],
                        "Number of Contacts": self.num_contacts[chain_pair],
                        "Average PAE ij": self.pairwise_avg_pae[chain_pair]["ij"],
                        "Average PAE ji": self.pairwise_avg_pae[chain_pair]["ji"],
                        "Average iPAE ij": self.pairwise_avg_ipae[chain_pair]["ij"] if chain_pair in self.pairwise_avg_ipae else np.nan,
                        "Average iPAE ji": self.pairwise_avg_ipae[chain_pair]["ji"] if chain_pair in self.pairwise_avg_ipae else np.nan,
                        "Minimum PAE ij": self.pairwise_min_pae[chain_pair]["ij"] if chain_pair in self.pairwise_min_pae else np.nan,
                        "Minimum PAE ji": self.pairwise_min_pae[chain_pair]["ji"] if chain_pair in self.pairwise_min_pae else np.nan,
                        "Average ipLDDT chain1": self.pairwise_avg_iplddt[chain_pair].get(chain1, np.nan),
                        "Average ipLDDT chain2": self.pairwise_avg_iplddt[chain_pair].get(chain2, np.nan),
                        "Chain Type 1": "IDR" if chain1 in self.idr_chains else "R",
                        "Chain Type 2": "IDR" if chain2 in self.idr_chains else "R",
                    })

            else:

                if self.symmetric_pae:
                    for res1_idx, res2_idx in self.interface_res_pairs[chain_pair]:
                        ipae_val = (
                            self.pairwise_ipae[chain_pair]["ij"].get((res1_idx, res2_idx), np.nan) +
                            self.pairwise_ipae[chain_pair]["ji"].get((res2_idx, res1_idx), np.nan)
                        ) / 2
                        chain_pairwise_assessment_rows.append({
                            "Chain Pair": f"{chain1}-{chain2}",
                            "Protein Name 1": self.protein_chain_map.get(chain1, "Unknown"),
                            "Protein Name 2": self.protein_chain_map.get(chain2, "Unknown"),
                            "Residue Pair": f"{res1_idx}-{res2_idx}",
                            "iPAE": ipae_val,
                            "ipLDDT res1": self.per_chain_iplddt[chain1].get(res1_idx, np.nan),
                            "ipLDDT res2": self.per_chain_iplddt[chain2].get(res2_idx, np.nan),
                            "Chain Type 1": "IDR" if chain1 in self.idr_chains else "R",
                            "Chain Type 2": "IDR" if chain2 in self.idr_chains else "R",
                        })
                else:
                    for res1_idx, res2_idx in self.interface_res_pairs[chain_pair]:
                        ipae_ij = self.pairwise_ipae[chain_pair]["ij"].get((res1_idx, res2_idx), np.nan)
                        ipae_ji = self.pairwise_ipae[chain_pair]["ji"].get((res2_idx, res1_idx), np.nan)
                        chain_pairwise_assessment_rows.append({
                            "Chain Pair": f"{chain1}-{chain2}",
                            "Protein Name 1": self.protein_chain_map.get(chain1, "Unknown"),
                            "Protein Name 2": self.protein_chain_map.get(chain2, "Unknown"),
                            "Residue Pair": f"{res1_idx}-{res2_idx}",
                            "iPAE ij": ipae_ij,
                            "iPAE ji": ipae_ji,
                            "ipLDDT res1": self.per_chain_iplddt[chain1].get(res1_idx, np.nan),
                            "ipLDDT res2": self.per_chain_iplddt[chain2].get(res2_idx, np.nan),
                            "Chain Type 1": "IDR" if chain1 in self.idr_chains else "R",
                            "Chain Type 2": "IDR" if chain2 in self.idr_chains else "R",
                        })

        overall_assessment_keys = {
            "Number of Chains": "num_chains",
            "Number of Interacting Chain Pairs": "num_interacting_chain_pairs",
            "Number of Interface Residues": "num_interface_residues",
            "Number of Contacts": "num_contacts",
            "Average ipLDDT": "avg_iplddt",
            "Average IDR ipLDDT": "avg_idr_iplddt",
            "Average iPAE ij": "avg_ipae_ij",
            "Average iPAE ji": "avg_ipae_ji",
        }

        for col_head, key in overall_assessment_keys.items():
            if self.overall_assessment.get(key, np.nan) is not np.nan:
                overall_assessment_rows.append({
                    "Key": col_head,
                    "Value": self.overall_assessment.get(key)
                })

        chain_pairwise_assessment_df = pd.DataFrame(chain_pairwise_assessment_rows)
        chainwise_assessment_df = pd.DataFrame(chain_wise_assessment_rows)
        overall_assessment_df = pd.DataFrame(overall_assessment_rows)

        df_dict = {
            "chain_pairwise_assessment": chain_pairwise_assessment_df,
            "chainwise_assessment": chainwise_assessment_df,
            "overall_assessment": overall_assessment_df,
        }

        for k, df_ in df_dict.items():
            df_dict[k] = df_.fillna(np.nan)
            df_dict[k] = df_.map(lambda x: round(x, 2) if isinstance(x, (int, float)) else x)

        with pd.ExcelWriter(self.save_path, engine='openpyxl', mode='w') as writer:
            for sheet_name, df in df_dict.items():

                if df.empty:
                    warnings.warn(f"Skipping empty DataFrame for sheet: {sheet_name}")
                    continue

                df.to_excel(
                    writer,
                    sheet_name=sheet_name,
                    index=False,
                )

    # @time_it
    def get_unique_chains(self):
        """Get unique chains in the rigid body.

        Returns:
            unique_chains (list): List of unique chain IDs in the rigid body.
        """

        unique_chains = [
            chain_id
            for chain_id in self.rb_dict.keys()
            if len(self.rb_dict[chain_id]) > 0
        ]

        return unique_chains

    # @time_it
    def get_chain_pairs(self):
        """Get all unique chain pairs in the rigid body.

        Returns:
            chain_pairs (list): List of tuples containing unique chain pairs.
            Each tuple contains two chain IDs.
        """

        chain_pairs = list(combinations(self.unique_chains, 2))

        return [tuple(pair) for pair in chain_pairs]

    # @time_it
    def get_rb_res_binary_map(self, lengths_dict):
        """Get a binary map of residues in the rigid body.

        Returns:
            rb_res_binary_map (np.ndarray): A binary map of residues in the rigid body.
            The shape is (total_length, total_length) where total_length is the sum of lengths of all chains.
            The value is 1 if the residue is part of the rigid body, 0 otherwise.
        """

        total_len = lengths_dict.get("total", 0)
        rb_res_binary_map = np.zeros((total_len, total_len), dtype=int)
        all_rb_interface_res_idxs = []

        for chain_id, atom_name_token_list in self.rb_dict.items():

            res_idxs = [
                self.num_to_idx[chain_id][token_num][atom_name]
                for atom_name, token_num in atom_name_token_list
            ]
            all_rb_interface_res_idxs.extend(res_idxs)

        all_rb_interface_res_idxs = np.unique(all_rb_interface_res_idxs)

        rb_res_binary_map[
            np.ix_(all_rb_interface_res_idxs, all_rb_interface_res_idxs)
        ] = 1

        return rb_res_binary_map

    # @time_it
    def get_rb_res_pairs(self):
        """Get all unique residue pairs in the rigid body.

        Returns:
            rb_res_pairs (defaultdict): A dictionary where keys are chain pairs (tuples) and values are lists of residue index pairs.
            Each residue index pair is a tuple of indices from the two chains in the rigid body.
        """

        rb_res_pairs = defaultdict(list)

        for chain_pair in self.chain_pairs:

            chain1, chain2 = chain_pair

            atom_name_token_list1 = self.rb_dict[chain1]
            atom_name_token_list2 = self.rb_dict[chain2]

            res1_idxs = [
                self.num_to_idx[chain1][token_num][atom_name]
                for atom_name, token_num in atom_name_token_list1
            ]
            res2_idxs = [
                self.num_to_idx[chain2][token_num][atom_name]
                for atom_name, token_num in atom_name_token_list2
            ]

            # Create pairs of residues from the two chains
            pairs = list(product(res1_idxs, res2_idxs))
            rb_res_pairs[chain_pair].extend(pairs)

        return rb_res_pairs

    # @time_it
    def get_interface_res_pairs(
        self,
        contact_map: np.ndarray,
    ):
        """ Get interface residue pairs from the contact map.

        Args:
            contact_map (np.ndarray): A binary contact map where 1 indicates a contact between residues and 0 indicates no contact.

        Returns:
            interface_res_pairs (defaultdict): A dictionary where keys are chain pairs (tuples) and values are lists of residue index pairs.
        """

        interface_res_pairs = defaultdict(list)
        contacting_res_indices = np.argwhere(contact_map == 1)

        for chain1, chain2 in tqdm(self.chain_pairs):

            atom_name_token_list1 = self.rb_dict[chain1]
            atom_name_token_list2 = self.rb_dict[chain2]

            res1_idxs = [
                self.num_to_idx[chain1][token_num][atom_name]
                for atom_name, token_num in atom_name_token_list1
            ]
            res2_idxs = [
                self.num_to_idx[chain2][token_num][atom_name]
                for atom_name, token_num in atom_name_token_list2
            ]

            mask1 = np.isin(contacting_res_indices[:, 0], res1_idxs)
            mask2 = np.isin(contacting_res_indices[:, 1], res2_idxs)

            mask = mask1 & mask2
            contacting_res_pairs = set(
                map(tuple, contacting_res_indices[mask])
            )

            if len(contacting_res_pairs) > 0:
                interface_res_pairs[(chain1, chain2)].extend(
                    list(contacting_res_pairs)
                )

        return interface_res_pairs

    # @time_it
    def get_per_chain_interface_residues(self):
        """Get interface residues for each chain.

        Returns:
            per_chain_interface_residues (defaultdict): A dictionary where keys are chain IDs and values are lists of residue indices.
            Each list contains the indices of residues that are part of any of the interfaces that the chain is involved in.
        """

        per_chain_interface_residues = defaultdict(list)

        for chain_pair, interacting_res_pairs in self.interface_res_pairs.items():

            chain1, chain2 = chain_pair

            for res1_idx, res2_idx in interacting_res_pairs:

                per_chain_interface_residues[chain1].append(
                    res1_idx
                ) if res1_idx not in per_chain_interface_residues[chain1] else None

                per_chain_interface_residues[chain2].append(
                    res2_idx
                ) if res2_idx not in per_chain_interface_residues[chain2] else None

        return per_chain_interface_residues

    # @time_it
    def get_num_interface_residues(self):
        """Get the number of interface residues for each chain pair.

        Returns:
            num_interface_residues (defaultdict): A dictionary where keys are chain pairs (tuples) and values are the number of unique interface residues.
            Each key is a tuple of two chain IDs, and the value is the count of unique residues that interact between those chains.
        """

        num_interface_residues = defaultdict(int)

        for chain_pair, interacting_res_pairs in self.interface_res_pairs.items():

            unique_interface_residues = np.unique(np.array(interacting_res_pairs).flatten())
            num_interface_residues[chain_pair] = len(unique_interface_residues)

        return num_interface_residues

    # @time_it
    def get_num_contacts(self):
        """Get the number of contacts for each chain pair.

        Returns:
            num_contacts (defaultdict): A dictionary where keys are chain pairs (tuples) and values are the number of contacts.
            Each key is a tuple of two chain IDs, and the value is the count of contacts between those chains.
        """

        num_contacts = defaultdict(int)

        for chain_pair, interacting_res_pairs in self.interface_res_pairs.items():

            num_contacts[chain_pair] = len(interacting_res_pairs)

        return num_contacts

    # @time_it
    def get_per_chain_plddt(self, plddt_list):
        """ Get per-chain pLDDT scores from a list of pLDDT scores.

        Args:
            plddt_list (list): A list of pLDDT scores for all residues in the structure.

        Returns:
            per_chain_plddt (defaultdict): A dictionary where keys are chain IDs and values are numpy arrays of pLDDT scores for residues in that chain.
        """

        per_chain_plddt = defaultdict(np.ndarray)

        for chain_id, atom_name_token_list in self.rb_dict.items():

            res_idxs = [
                self.num_to_idx[chain_id][token_num][atom_name]
                for atom_name, token_num in atom_name_token_list
            ]

            plddt_scores = np.array(plddt_list)[res_idxs]
            per_chain_plddt[chain_id] = plddt_scores

        return per_chain_plddt

    # @time_it
    def get_per_chain_avg_plddt(self):
        """ Get the average pLDDT score for each chain.

        Returns:
            per_chain_avg_plddt (dict): A dictionary where keys are chain IDs and values are the average pLDDT scores for that chain.
        """

        return {
            chain_id: np.mean(plddt_scores)
            for chain_id, plddt_scores in self.per_chain_plddt.items()
        }

    # @time_it
    def get_per_chain_iplddt(self, plddt_list):
        """ Get per-chain ipLDDT scores from a list of pLDDT scores.

        Args:
            plddt_list (list): A list of pLDDT scores for all residues in the structure.

        Returns:
            per_chain_iplddt (defaultdict): A dictionary where keys are chain IDs and values are dictionaries mapping residue indices to their pLDDT scores.
        """

        per_chain_iplddt = defaultdict(dict)

        for chain_id, interface_res_idxs in self.per_chain_interface_residues.items():

            for res_idx in interface_res_idxs:

                per_chain_iplddt[chain_id][res_idx] = plddt_list[res_idx]

        return per_chain_iplddt

    # @time_it
    def get_per_chain_average_iplddt(self):
        """ Get the average ipLDDT score for each chain.

        Returns:
            per_chain_avg_iplddt (dict): A dictionary where keys are chain IDs and values are the average ipLDDT scores for that chain.
        """

        return {
            chain_id: np.mean(list(iplddt_scores.values()))
            for chain_id, iplddt_scores in self.per_chain_iplddt.items()
        }

    # @time_it
    def get_pairwise_pae(self, pae):
        """ Get pairwise PAE values for each chain pair.

        Args:
            pae (np.ndarray): A 2D numpy array representing the predicted aligned error (PAE) matrix.

        Returns:
            pairwise_pae (defaultdict): A dictionary where keys are chain pairs (tuples) and values are dictionaries containing PAE values for residue pairs.
        """

        pairwise_pae = defaultdict(np.ndarray)

        for chain_pair in self.chain_pairs:

            rb_chain_pair_res = self.rb_res_pairs[chain_pair]

            rb_pae_vals_ij = [
                pae[res1_idx, res2_idx] for res1_idx, res2_idx in rb_chain_pair_res
            ]

            rb_pae_vals_ji = [
                pae[res2_idx, res1_idx] for res1_idx, res2_idx in rb_chain_pair_res
            ]

            if len(rb_pae_vals_ij) > 0:
                pairwise_pae[chain_pair] = {
                    "ij": rb_pae_vals_ij,
                    "ji": rb_pae_vals_ji,
                }

        return pairwise_pae

    # @time_it
    def get_pairwise_avg_pae(self, symmetric_pae: bool = True):
        """ Get the average PAE for each chain pair.

        Args:
            symmetric_pae (bool, optional): If True, calculates the average PAE symmetrically for both directions (ij and ji).

        Returns:
            pairwise_avg_pae (defaultdict): A dictionary where keys are chain pairs (tuples) and values are the average PAE values.
        """

        if symmetric_pae:
            pairwise_avg_pae = defaultdict(float)
        else:
            pairwise_avg_pae = defaultdict(dict)

        for chain_pair in self.chain_pairs:
            if symmetric_pae:
                pairwise_avg_pae[chain_pair] = (
                    np.mean(
                        self.pairwise_pae[chain_pair]["ij"] +
                        self.pairwise_pae[chain_pair]["ji"]
                    ) / 2
                )
            else:
                pairwise_avg_pae[chain_pair]["ij"] = np.mean(self.pairwise_pae[chain_pair]["ij"])
                pairwise_avg_pae[chain_pair]["ji"] = np.mean(self.pairwise_pae[chain_pair]["ji"])

        return pairwise_avg_pae

    # @time_it
    def get_pairwise_min_pae(self, symmetric_pae: bool = True):
        """Get the minimum PAE for each chain pair.

        Args:
            symmetric_pae (bool, optional): If True, calculates the minimum PAE symmetrically for both directions (ij and ji).

        Returns:
            pairwise_min_pae (defaultdict): A dictionary where keys are chain pairs (tuples) and values are the minimum PAE values.
            If symmetric_pae is True, the minimum PAE is calculated as the minimum of both directions (ij and ji).
            If symmetric_pae is False, the minimum PAE is calculated separately for each direction.
        """

        if symmetric_pae:
            pairwise_min_pae = defaultdict(float)
        else:
            pairwise_min_pae = defaultdict(dict)

        for chain_pair, pae_dict in self.pairwise_pae.items():
            if symmetric_pae:
                pairwise_min_pae[chain_pair] = np.min(
                    [np.min(pae_dict["ij"]), np.min(pae_dict["ji"])]
                )
            else:
                pairwise_min_pae[chain_pair]["ij"] = np.min(pae_dict["ij"])
                pairwise_min_pae[chain_pair]["ji"] = np.min(pae_dict["ji"])

        return pairwise_min_pae

    # @time_it
    def get_pairwise_ipae(self, pae):
        """ Get pairwise iPAE values for each chain pair.

        Args:
            pae (np.ndarray): A 2D numpy array representing the predicted aligned error (PAE) matrix.

        Returns:
            pairwise_ipae (defaultdict): A dictionary where keys are chain pairs (tuples) and values are dictionaries containing iPAE values for residue pairs.
        """

        pairwise_ipae = defaultdict(dict)

        for chain_pair, interacting_res_pairs in self.interface_res_pairs.items():

            pairwise_ipae[chain_pair] = {
                "ij" : {
                    (res1_idx, res2_idx): pae[res1_idx, res2_idx]
                    for res1_idx, res2_idx in interacting_res_pairs
                },
                "ji" : {
                    (res2_idx, res1_idx): pae[res2_idx, res1_idx]
                    for res1_idx, res2_idx in interacting_res_pairs
                }
            }

        return pairwise_ipae

    # @time_it
    def get_pairwise_avg_ipae(self, symmetric_pae: bool = True):
        """ Get the average iPAE for each chain pair.

        Args:
            symmetric_pae (bool, optional): If True, calculates the average iPAE symmetrically for both directions (ij and ji).

        Returns:
            pairwise_avg_ipae (defaultdict): A dictionary where keys are chain pairs (tuples) and values are the average iPAE values.
        """

        if symmetric_pae:
            pairwise_avg_ipae = defaultdict(float)
        else:
            pairwise_avg_ipae = defaultdict(dict)

        for chain_pair, ipae_dict in self.pairwise_ipae.items():

            if symmetric_pae:
                pairwise_avg_ipae[chain_pair] = (
                    np.mean(list(ipae_dict["ij"].values()) + list(ipae_dict["ji"].values())) / 2
                )
            else:
                pairwise_avg_ipae[chain_pair]["ij"] = np.mean(list(ipae_dict["ij"].values()))
                pairwise_avg_ipae[chain_pair]["ji"] = np.mean(list(ipae_dict["ji"].values()))

        return pairwise_avg_ipae

    # @time_it
    def get_pairwise_avg_iplddt(self):
        """ Get the average ipLDDT for each chain pair.

        Returns:
            pairwise_avg_iplddt (defaultdict): A dictionary where keys are chain pairs (tuples) and values are dictionaries containing average ipLDDT values for each chain in the pair.
        """

        pairwise_avg_iplddt = defaultdict(dict)

        for chain_pair, interacting_res_pairs in self.interface_res_pairs.items():

            chain1, chain2 = chain_pair

            iplddt1_values = [
                self.per_chain_iplddt[chain1].get(res1_idx, np.nan)
                for res1_idx, res2_idx in interacting_res_pairs
            ]

            iplddt2_values = [
                self.per_chain_iplddt[chain2].get(res2_idx, np.nan)
                for res1_idx, res2_idx in interacting_res_pairs
            ]

            pairwise_avg_iplddt[chain_pair][chain1] = np.mean(iplddt1_values)
            pairwise_avg_iplddt[chain_pair][chain2] = np.mean(iplddt2_values)

        return pairwise_avg_iplddt

    # @time_it
    def get_overall_assessment(self):
        """ Get overall assessment of the rigid body.

        Returns:
            overall_assessment (dict): A dictionary containing overall statistics about the rigid body.
            It includes the number of chains, number of interacting chain pairs, number of interface residues,
            number of contacts, average ipLDDT, average IDR ipLDDT, average iPAE ij, and average iPAE ji.
        """

        overall_assessment = {}

        overall_assessment["num_chains"] = len(self.unique_chains)

        overall_assessment["num_interacting_chain_pairs"] = len(self.interface_res_pairs)

        overall_assessment["num_interface_residues"] = sum(
            len(res_list)
            for res_list in self.per_chain_interface_residues.values()
        )

        overall_assessment["num_contacts"] = sum(
            len(contact_pairs)
            for contact_pairs in self.interface_res_pairs.values()
        )

        global_iplddt_scores = [
            iplddt
            for iplddt_scores in self.per_chain_iplddt.values()
            for iplddt in iplddt_scores.values()
        ]

        overall_assessment["avg_iplddt"] = (
            np.mean(global_iplddt_scores) if global_iplddt_scores else np.nan
        )

        global_idr_iplddt_scores = [
            iplddt
            for chain_id, iplddt_scores in self.per_chain_iplddt.items()
            for iplddt in iplddt_scores.values()
            if chain_id in self.idr_chains
        ]

        overall_assessment["avg_idr_iplddt"] = (
            np.mean(global_idr_iplddt_scores) if global_idr_iplddt_scores else np.nan
        )

        global_ipae_ij_scores = [
            ipae
            for ipae_dict in self.pairwise_ipae.values()
            for ipae in ipae_dict["ij"].values()
        ]

        global_ipae_ji_scores = [
            ipae
            for ipae_dict in self.pairwise_ipae.values()
            for ipae in ipae_dict["ji"].values()
        ]

        overall_assessment["avg_ipae_ij"] = (
            np.mean(global_ipae_ij_scores) if global_ipae_ij_scores else np.nan
        )

        overall_assessment["avg_ipae_ji"] = (
            np.mean(global_ipae_ji_scores) if global_ipae_ji_scores else np.nan
        )

        return overall_assessment
