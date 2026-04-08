"""
[rigid_bodies](https://github.com/isblab/af_pipeline/tree/main/af_pipeline/rigid_bodies)
=============================================================

- Module to extract confident region from AlphaFold predictions

- This module contains classes and methods for extracting **confidently predicted
  regions (rigid bodies)** from AlphaFold predictions based on the assessment
  metrics provided by AlphaFold: PAE and pLDDT.

- It uses graph-based community detection algorithms[^tristan] to identify
  communities of residues that are confidently predicted; a set of residues
  within each community satisfying a pLDDT cutoff constitute a rigid body.

[^tristan]: Tristan Croll, "Graph-based community clustering approach to extract protein domains from a predicted aligned error matrix": https://github.com/tristanic/pae_to_domains

## Organization

The module is organized into the following submodules:

- **rigid_bodies**: Contains the [`RigidBodies`](https://github.com/isblab/af_pipeline/blob/89b33286ff81b1f4075c65fd6d5dc88296d83c15/af_pipeline/rigid_bodies/rigid_bodies.py#L71) class, which can be used for extracting confidently predicted regions from AlphaFold predictions based on the assessment metrics provided by AlphaFold: PAE and pLDDT.

- **rigid_body_assessment**: Contains the [`RigidBodyAssessment`](https://github.com/isblab/af_pipeline/blob/89b33286ff81b1f4075c65fd6d5dc88296d83c15/af_pipeline/rigid_bodies/rigid_body_assessment.py#L1279) class, which can be used for assessing the quality of the extracted rigid bodies using various metrics such as average pLDDT, average PAE, and contact probabilities.

```mermaid
---
config:
  class:
    hideEmptyMembersBox: true
---
classDiagram

    class RigidBodies {
        + str library
        + float pae_cutoff
        + float plddt_cutoff
        + Optional[float] plddt_cutoff_idr
        + Optional[List] idr_chains
        + Optional[int] pae_power
        + Optional[float] resolution
        + Optional[int] random_seed
        - \_\_init__(self, library, plddt_cutoff, pae_cutoff, **kwargs) None
        + check_is_set_up(self)
        + set_attributes_from(self, instance)
        + extract_rigid_bodies(self, pae_matrix, min_res, min_proteins, plddt_filter) list[dict[str, list[tuple[str, int]]]]
        - \_convert_domain_to_dict(self, pseudo_domain) dict[str, list[tuple[str, int]]]
        - \_filter_by_plddt(self, domain_dict, token_plddts) dict[str, list[tuple[str, int]]]
        - @staticmethod \_filter_by_domain_size(rb_dict, min_res, min_proteins) dict$
        - @staticmethod \_keep_residue_numbers_only(rigid_bodies) list[dict[str, list[int]]]$
        + save_rigid_bodies(self, domains, output_dir, rb_out_fmt, save_structure, rb_struct_fmt, filter_struct_by_plddt, protein_chain_map)
        + show_rigid_bodies_on_pae_matrix(self, domains, output_dir)
        + assess_rigid_bodies(self, domains, output_dir, protein_chain_map, symmetric_pae, as_average)
    }

    link RigidBodies "rigid_bodies/rigid_bodies.html#RigidBodies" "link to RigidBodies class documentation"

    RigidBodies ..|> PAEPatches

    link PAEPatches "rigid_bodies/rigid_bodies.html#PAEPatches" "link to PAEPatches class documentation"

    RigidBodies ..> pae_to_domains

    link pae_to_domains "pae_to_domains/pae_to_domains.html" "link to pae_to_domains module documentation"

```

```mermaid
---
config:
  class:
    hideEmptyMembersBox: true
---
classDiagram
    class RigidBodyAssessment {
        + Dict[str, List[Tuple[str, int]]] rb_dict
        + Dict[str, Dict[int, Dict[str, int]]] num_to_idx
        + Dict[int, Dict[str, str | int]] idx_to_num
        + np.ndarray contact_map
        + np.ndarray plddt_list
        + np.ndarray pae
        + Dict[str, int] lengths_dict
        + bool symmetric_pae
        + bool as_average
        + List[str] idr_chains
        + Dict[str, str] protein_chain_map
        + List[str] unique_chains
        + List[Tuple[str, str]] chain_pairs
        + np.ndarray rb_mask
        + dict overall_assessment
        + str save_path
        - \_\_init__(self, rb_dict, as_average, symmetric_pae, **kwargs) None
        + check_is_set_up(self)
        + perform_assessment(self)
        + set_attributes_from(self, instance)
        + save_rb_assessment(self, rb_c_assess, rb_cp_assess, overall_assessment, save_path)
        + get_overall_assessment(self, rb_c_assess, rb_cp_assess)
    }

    link RigidBodyAssessment "rigid_bodies/rigid_body_assessment.html#RigidBodyAssessment" "link to RigidBodyAssessment class documentation"

    RigidBodyAssessment .. RigidBodies

    link RigidBodies "rigid_bodies/rigid_bodies.html#RigidBodies" "link to RigidBodies class documentation"

    RigidBodyAssessment ..|> RigidBodyChainAssessment

    link RigidBodyChainAssessment "rigid_bodies/rigid_body_assessment.html#RigidBodyChainAssessment" "link to RigidBodyChainAssessment class documentation"

    RigidBodyAssessment ..|> RigidBodyChainPairAssessment

    link RigidBodyChainPairAssessment "rigid_bodies/rigid_body_assessment.html#RigidBodyChainPairAssessment" "link to RigidBodyChainPairAssessment class documentation"

    RigidBodyAssessment ..|> _Mask

```

```mermaid
---
config:
  class:
    hideEmptyMembersBox: true
---
classDiagram

    class RigidBodyChainAssessment {
        + bool as_average
        + List[str] unique_chains
        + Dict[int, Dict[str, int]] idx_to_num
        + np.ndarray chain_mask_stack_1d
        + np.ndarray rb_mask_1d
        + np.ndarray contact_map_mask_1d
        + np.ndarray plddt_list
        + List[str] idr_chains
        + Dict[str, str] protein_chain_map
        + Dict[str, list | float] per_chain_plddt
        + Dict[str, list | float] per_chain_iplddt
        + Dict[str, List[int] | int] per_chain_interface_res
        - \_\_init__(self, _mask, as_average) None
        + get_per_chain_plddt(self, only_avg, only_interface) Dict[str, float | List[float]]
        + get_per_chain_interface_residues(self, only_count) Dict[str, int | List[int]]
        + get_per_chain_residues(self, only_count) Dict[str, int | List[int]]
        + get_chain_attr(self, chain_id, attr_name) float | str | int
        + get_res_attr(self, chain_id, res_idx, attr_name) float | str | int
        + get_chain_assessment(self) pd.DataFrame
    }

    link RigidBodyChainAssessment "rigid_bodies/rigid_body_assessment.html" "link to RigidBodyChainAssessment class documentation"
```

```mermaid
---
config:
  class:
    hideEmptyMembersBox: true
---
classDiagram

    class RigidBodyChainPairAssessment {
        + bool as_average
        + List[str] unique_chains
        + List[Tuple[str, str]] chain_pairs
        + Dict[int, Dict[str, int]] idx_to_num
        + bool symmetric_pae
        + np.ndarray chain_mask_stack_1d
        + np.ndarray chain_mask_stack_2d
        + np.ndarray chain_pair_mask_stack_1d
        + np.ndarray chain_pair_mask_stack_2d
        + np.ndarray rb_mask_1d
        + np.ndarray rb_mask_2d
        + np.ndarray contact_map_mask_1d
        + np.ndarray contact_map_mask_2d
        + np.ndarray plddt_list
        + np.ndarray pae
        + np.ndarray avg_pae
        + List[str] idr_chains
        + Dict[str, str] protein_chain_map
        + Dict[Tuple[str, str], List[int] | int] chain_pair_interface_res
        + Dict[Tuple[str, str], List[int] | int] chain_pair_contacts
        + Dict[Tuple[str, str], List[float] | float] chain_pair_iplddt
        + Dict[Tuple[str, str], List[float] | float] chain_pair_pae
        + Dict[Tuple[str, str], List[float] | float] chain_pair_ipae
        + Dict[Tuple[str, str], List[float] | float] chain_pair_pae_ij
        + Dict[Tuple[str, str], List[float] | float] chain_pair_pae_ji
        + Dict[Tuple[str, str], List[float] | float] chain_pair_ipae_ij
        + Dict[Tuple[str, str], List[float] | float] chain_pair_ipae_ji
        - \_\_init__(self, _mask, as_average) None
        + get_chain_pair_attr(self, chain_pair, attr_name)
        + get_res_pair_attr(self, chain_pair, res_pair, attr_name)
        + get_chain_pair_interface(self, per_chain, only_count) Dict[Tuple[str, str], List[int] | int | Tuple[int, int]]
        + get_chain_pair_plddt(self, only_avg, only_interface)
        + get_chain_pair_pae(self, only_avg, only_interface, symmetric)
        + get_chain_pair_assessment(self)
    }

    link RigidBodyChainPairAssessment "rigid_bodies/rigid_body_assessment.html" "link to RigidBodyChainPairAssessment class documentation"
```

```mermaid
---
config:
  class:
    hideEmptyMembersBox: true
---
classDiagram

    class PAEPatches {
        + dict num_to_idx
        + np.ndarray pae
        + dict lengths_dict
        + int rb_idx
        + dict | None af_offset
        - \_\_init__(self, num_to_idx, pae, lengths_dict, rb_idx, af_offset) None
        + extract_pae_patches(self, rb_dict) list[list]
        + plot_pae_patches(self, patches, output_dir)
    }

    class _Mask {
        + Dict[str, List[int]] rb_dict
        + Dict[str, Dict[int, Dict[str, int]]] num_to_idx
        + Dict[int, Dict[str, int]] idx_to_num
        + Dict[str, int] lengths_dict
        + np.ndarray contact_map
        + np.ndarray plddt_list
        + bool symmetric_pae
        + List[str] idr_chains
        + Dict[str, str] protein_chain_map
        + np.ndarray pae
        + np.ndarray avg_pae
        + np.ndarray contact_map_mask_2d
        + np.ndarray contact_map_mask_1d
        + np.ndarray rb_mask_2d
        + np.ndarray rb_mask_1d
        + np.ndarray chain_mask_stack_1d
        + np.ndarray chain_mask_stack_2d
        + np.ndarray chain_pair_mask_stack_1d
        + np.ndarray chain_pair_mask_stack_2d
        - \_\_init__(self, rb_dict, pae, avg_pae, plddt_list, contact_map, lengths_dict, num_to_idx, idx_to_num, symmetric_pae, idr_chains, protein_chain_map) None
        + get_unique_chains(self) List[str]
        + get_chain_pairs(self) List[Tuple[str, str]]
        + transform_rb_dict_to_idxs(self, rb_dict) Dict[str, List[int]]
        + get_rb_mask(self, lengths_dict, dimensions) np.ndarray
        + get_chain_mask(self, chain_id, lengths_dict, dimensions) np.ndarray
        + get_chain_pair_mask(self, chain_id_1, chain_id_2, lengths_dict, dimensions) np.ndarray
        + get_chain_mask_stack(self, dimensions) np.ndarray
        + get_chain_pair_mask_stack(self, dimensions) np.ndarray
        + @staticmethod sanity_check_mask_dimensions(dimensions)$
    }

    link PAEPatches "rigid_bodies/rigid_body_assessment.html" "link to PAEPatches class documentation"

    link _Mask "rigid_bodies/rigid_body_assessment.html" "link to _Mask class documentation"
```

## Prerequisites

- **best_af_predictions.json**: obtained by running [`rank_af_predictions.py`](https://github.com/isblab/af_pipeline/tree/main/examples/rank_af_predictions.py).
  This file contains information about the best prediction for each job.

- **Structure predictions**: obtained by running the prediction jobs on AlphaFold server
  or AlphaFold2 or ColabFold.

## Usage

- Please refer to the [examples directory](https://github.com/isblab/af_pipeline/tree/main/examples) for sample scripts and
  config file.

- Use the following command to run the example script for extracting rigid bodies:
```
python extract_rigid_bodies.py \\
    -i ./input/best_af_predictions.json \\
    -o ./output/rigid_bodies \\
    --plddt_cutoff 70 \\
    --pae_cutoff 5 \\
    --pae_power 1 \\
    --resolution 0.5 \\
    --library networkx \\
    --min_res 10 \\
    --min_proteins 1 \\
    --apply_plddt_filter True
```

## Workflows

- Workflow for extracting rigid bodies from AlphaFold predictions:

```mermaid
sequenceDiagram
    autonumber
    participant User
    participant Initialize
    participant RigidBodies

    rect rgb(210, 250, 200)
    note over User, RigidBodies: Set up RigidBodies instance.
    User->>Initialize:
    User->>RigidBodies: parameters & instance of Initialize
    end

    rect rgb(240, 255, 255)
    note over User, RigidBodies: Extract and save rigid bodies.
    User->>RigidBodies: extract_rigid_bodies()
    User->>RigidBodies: save_rigid_bodies()
    end

    rect rgb(245, 255, 200)
    note over User, RigidBodies: Assess and visualize rigid bodies.
    User->>RigidBodies: assess_rigid_bodies()
    User->>RigidBodies: show_rigid_bodies_on_pae_matrix()
    end
```

"""