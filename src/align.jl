# Align
# =====
#
# Interface function for alignments.


include("read_inputs.jl")
include("local_align.jl")
include("global_align.jl")
include("align_viz.jl")


"""
    align(fp1::AbstractString, fp2::AbstractString)
    align(; thresh::Union{Int,Bool}=false, submatfp::Union{AbstractString,Bool}=false,
          gap_open::Int=-12, gap_extend::Int=-4, dedup::Bool=true,
          dedup_method::AbstractString="score", global_align::Bool=true,
          global_gap_open::Int=-12, global_gap_extend::Int=-4, global_end_gap_open::Int=0,
          global_end_gap_extend::Int=0, print::Bool=true, figure::Bool=true,
          figure_type::AbstractString="interactive", figurewidth::Int=1000)

Identify local alignments between sequence pair at filepaths `fp1` and `fp2`.  
  
Optionally, visualize these alignments through a printout or graph; and optionally include
a global alignment in this graph for comparison.

...
## Arguments:
- `fp1::AbstractString`: fasta file containing sequence 1.
- `fp2::AbstractString`: fasta file containing sequence 2.
- `thresh::Union{Int,Bool}=false`: score threshold for the suboptimal local alignments;
    calculated if set to `false`.
- `submatfp::Union{AbstractString,Bool}=false`: filepath for substitution matrix.
- `gap_open::Int=-12`: penalty for gap opening.
- `gap_extend::Int=-4`: penalty for gap extension.
- `dedup::Bool=true`: whether to deduplicate the local alignments.
- `dedup_method::AbstractString="score"`: method of local alignment deduplication; for
    `score`, keeps all highest-scoring alignments along each unique path; for `length`,
    keeps only the longest alignment along each unique path.
- `global_align::Bool=true`: whether to identify a global alignment.
- `global_gap_open::Int=-12`: penalty for gap opening in the global alignment.
- `global_gap_extend::Int=-4`: penalty for gap extension in the global alignment.
- `global_end_gap_open::Int=0`: penalty for end gap opening in the global alignment.
- `global_end_gap_extend::Int=0`: penalty for end gap extension in the global alignment.
- `print::Bool=true`: whether to print out a color-coded view of the local alignments.
- `figure::Bool=true`: whether to return a graph of the alignments.
- `figure_type::AbstractString="interactive"`: the type of the graph; `interactive` or
    `static`.
- `figurewidth::Int=1000`: the figure width of the alignment graph.
"""
function align(fp1::AbstractString, fp2::AbstractString;
               thresh::Union{Int,Bool}=false, submatfp::Union{AbstractString,Bool}=false,
               gap_open::Int=-12, gap_extend::Int=-4,
               dedup::Bool=true, dedup_method::AbstractString="score",
               global_align::Bool=true, global_gap_open::Int=-12, global_gap_extend::Int=-4,
               global_end_gap_open::Int=0, global_end_gap_extend::Int=0,
               print::Bool=true, figure::Bool=true, figure_type::AbstractString="interactive",
               figurewidth::Int=1000)

    seqA_id, seqA = ReadInput(fp1)
    seqB_id, seqB = ReadInput(fp2)
    if submatfp!=false
        sub_header, sub_matrix = ReadSubMatrix(fp=submatfp)
    else
        sub_header, sub_matrix = BLOSUM62
    end

    sm, tm, am = LocalAlign(seqA, seqB, sub_header, sub_matrix;
        thresh=thresh, gap_open=gap_open, gap_extend=gap_extend,
        dedup=dedup, dedup_method=dedup_method)
    
    if global_align
        gsm, gtm, gam = GlobalAlign(seqA, seqB, sub_header, sub_matrix;
            gap_open=global_gap_open, gap_extend=global_gap_extend,
            end_gap_open=global_end_gap_open, end_gap_extend=global_end_gap_extend)
    end
    
    if print
        LocalAlignPrint(seqA_id, seqB_id, am2, sub_header, sub_matrix)
    end
    
    if figure
        if figure_type=="interactive"
            plt = iLocalAlignViz(seqA, seqB, seqA_id, seqB_id, am, sub_header, sub_matrix;
                figurewidth=figurewidth, global_alignment=gam)
        end
        if figure_type=="static"
            plt = sLocalAlignViz(seqA, seqB, seqA_id, seqB_id, am, sub_header, sub_matrix;
                figurewidth=figurewidth, global_alignment=gam)
        end
        return plt
    end
    
end
