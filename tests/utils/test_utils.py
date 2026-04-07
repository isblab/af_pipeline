from af_pipeline.utils.file_utils import read_fasta

def test_read_fasta():
    """Test the read_fasta function."""

    fasta_path = "tests/test_data/test_fasta.fasta"
    sequences = read_fasta(fasta_path)

    assert isinstance(sequences, dict), "Output should be a dictionary."
    assert len(sequences) == 2, "There should be 2 sequences in the output."
    assert sequences["protein_1"] == "ACYTQGGG", "Sequence for protein_1 does not match expected output."
    assert sequences["protein_2"] == "LLKPPPFCC", "Sequence for protein_2 does not match expected output."