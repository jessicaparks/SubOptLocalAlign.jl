# Read Inputs
# ===========
#
# Functions to read fasta file and substitution matrix inputs for alignments.


using DelimitedFiles
using FastaIO


"""
    read_input(fp::AbstractString)

Read input id and sequence from fasta file at path `fp`. Return a tuple of the
`sequence identifier` and the `sequence`.
"""
function read_input(fp::AbstractString)
    (seq_id, seq) = readfasta(fp)[1]
    seq_id = split(seq_id, " ")[1]
    seq = uppercase(seq)
    return (seq_id, seq)
end


"""
    read_sub_matrix(fp::AbstractString)

Read substitution matrix header and data from file at path `fp`. Return an array of the
headers and a matrix of the substitution scores.  
  
The input file should be formatted with the headers on the first line but not in the first
column and use space-delimitation.
"""
function read_sub_matrix(fp::AbstractString)
    (sub_matrix, sub_header) = readdlm(fp, header=true)
    sub_header = Dict([(uppercase(x[1]), i) for (i, x) in enumerate(sub_header)])
    return (sub_header, sub_matrix)
end

const BLOSUM62 = read_sub_matrix("data/BLOSUM62")
