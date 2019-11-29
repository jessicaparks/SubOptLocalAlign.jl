# Read Inputs
# ===========
#
# Functions to read fasta file and substitution matrix inputs for alignments.



"""
ReadInput(fp)
read input id and sequence from fasta file

Arguments:
fp (str): fasta filepath

Returns:
sequence_identifier
sequence
"""
function ReadInput(fp)
    (seq_id, seq) = readfasta(fp)[1]
    seq_id = split(seq_id, " ")[1]
    seq = uppercase(seq)
    return (seq_id, seq)
end


"""
ReadSubMatrix(fp)
read substitution matrix header and data from file

Arguments:
fp (str): filepath for sub matrix, default = BLOSUM62

Returns:
headers
matrix of substitution scores
"""
function ReadSubMatrix(fp="data/BLOSUM62.txt")
    (sub_matrix, sub_header) = readdlm(fp, header=true)
    sub_header = Dict([(uppercase(x[1]), i) for (i, x) in enumerate(sub_header)])
    return (sub_header, sub_matrix)
end