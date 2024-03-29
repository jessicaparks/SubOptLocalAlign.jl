# Align
# =====
#
# Interface function for alignments.


"""
    align(fp1::String, fp2::String; <keyword arguments>)
  
Identify local alignments between the sequence pair at filepaths `fp1` and `fp2`.  
  
Optionally, visualize these alignments through a printout or graph; and optionally include
a global alignment and/or a custom user-entered alignment in this graph for comparison.
The printout will print to the screen, and the graph will be returned as a plot.

---
## Arguments:
- `fp1::String`: fasta file containing sequence 1.
- `fp2::String`: fasta file containing sequence 2.
- `thresh::Union{Int,Bool}=false`: score threshold for the suboptimal local alignments;
    calculated if set to `false`.
- `submatfp::Union{String,Bool}=false`: filepath for substitution matrix.
- `gap_open::Int=-12`: penalty for gap opening.
- `gap_extend::Int=-4`: penalty for gap extension.
- `dedup::Bool=true`: whether to deduplicate the local alignments.
- `dedup_method::String="score"`: method of local alignment deduplication; for
    `score`, keeps all highest-scoring alignments along each unique path; for `length`,
    keeps only the longest alignment along each unique path.
- `global_alignment::Bool=true`: whether to identify a global alignment.
- `global_gap_open::Int=-12`: penalty for gap opening in the global alignment.
- `global_gap_extend::Int=-4`: penalty for gap extension in the global alignment.
- `global_end_gap_open::Int=0`: penalty for end gap opening in the global alignment.
- `global_end_gap_extend::Int=0`: penalty for end gap extension in the global alignment.
- `custom_alignment::Union{Array{Any,1},Bool}=false`: custom alignment to include in the
    output graph, entered by user.
- `print::Bool=true`: whether to print out a color-coded view of the local alignments.
- `figure::Bool=true`: whether to return a graph of the alignments.
- `figure_type::String="interactive"`: the type of the graph; `interactive` or
    `static`.
- `figurewidth::Int=1000`: the figure width of the alignment graph.
"""
function align(fp1::String, fp2::String;
               thresh::Union{Int,Bool}=false, submatfp::Union{String,Bool}=false,
               gap_open::Int=-12, gap_extend::Int=-4,
               dedup::Bool=true, dedup_method::String="score",
               global_alignment::Bool=true, global_gap_open::Int=-12, global_gap_extend::Int=-4,
               global_end_gap_open::Int=0, global_end_gap_extend::Int=0,
               custom_alignment::Union{Array{Any,1},Bool}=false,
               print::Bool=true, figure::Bool=true, figure_type::String="interactive",
               figurewidth::Int=1000)

    seqA_id, seqA = read_input(fp1)
    seqB_id, seqB = read_input(fp2)
    if submatfp!=false
        sub_header, sub_matrix = read_sub_matrix(submatfp)
    else
        sub_header, sub_matrix = BLOSUM62
    end

    println("aligning $seqA_id and $seqB_id ...")
    sm, tm, am = local_align(seqA, seqB, sub_header, sub_matrix;
        thresh=thresh, gap_open=gap_open, gap_extend=gap_extend,
        dedup=dedup, dedup_method=dedup_method)
    println(length(am), " local alignments found.\n")
    
    if global_alignment
        println("calculating global alignment ...")
        gsm, gtm, gam = global_align(seqA, seqB, sub_header, sub_matrix;
            gap_open=global_gap_open, gap_extend=global_gap_extend,
            end_gap_open=global_end_gap_open, end_gap_extend=global_end_gap_extend)
    end
    
    if print
        println("\n\n", "local alignments:", "\n")
        local_align_print(seqA_id, seqB_id, am, sub_header, sub_matrix)
        if global_alignment
            println("\n\n", "global alignment:", "\n")
            local_align_print(seqA_id, seqB_id, [gam], sub_header, sub_matrix)
        end
    end
    
    if figure
        if figure_type=="interactive"
            plt = i_local_align_viz(seqA, seqB, seqA_id, seqB_id, am, sub_header, sub_matrix;
                figurewidth=figurewidth, global_alignment=gam, custom_alignment=custom_alignment)
        end
        if figure_type=="static"
            plt = s_local_align_viz(seqA, seqB, seqA_id, seqB_id, am, sub_header, sub_matrix;
                figurewidth=figurewidth, global_alignment=gam, custom_alignment=custom_alignment)
        end
        return plt
    end
    
end
