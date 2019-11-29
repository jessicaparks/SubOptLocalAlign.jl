# Align
# =====
#
# Interface function for alignments.


include("read_inputs.jl")
include("local_align.jl")
include("global_align.jl")
include("align_viz.jl")


"""
"""
function align(fp1, fp2;
               thresh=false, submatfp="data/BLOSUM62",
               gap_open=-12, gap_extend=-4,
               dedup=true, dedup_method="score",
               global_align=true, global_gap_open=-12, global_gap_extend=-4,
               global_end_gap_open=0, global_end_gap_extend=0,
               print=true, figure=true, figure_type="interactive",
               figurewidth=1000)

    seqA_id, seqA = ReadInput(fp1)
    seqB_id, seqB = ReadInput(fp2)
    sub_header, sub_matrix = ReadSubMatrix(fp=submatfp)

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