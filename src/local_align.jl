# Local Align
# ===========
#
# Identifies sub-optimal local alignments for a pair of sequences.



"""
    traceback(align_pos::CartesianIndex{2},
        trace_matrix::Array{Array{Int,1},2}, score_matrix::Array{Int,2},
        seqA::String, seqB::String,
        alignments::Array{Tuple{String,String},1},
        scores::Array{Array{Int,1},1},
        positions::Array{Array{Tuple{Int,Int},1},1})
    traceback(; thresh::Int=0, dedup::Bool=true,
        current_alignment::Tuple{String,String}=("",""),
        current_scores::Array{Int,1}=Array{Int,1}(),
        current_positions::Array{Tuple{Int,Int},1}=Array{Tuple{Int,Int},1}(),
        maxpathscore::Union{Int,Bool}=false)
  
Trace path back through matrix to identify alignments. Push output to `alignments`,
`scores`, and `positions`.  
  
---
# Arguments:
- `align_pos::CartesianIndex{2}`: coordinates of the alignment end position, and then
    each current position as the traceback continues.
- `trace_matrix::Array{Array{Int,1},2}`: alignment traceback matrix.
- `score_matrix::Array{Int,2}`: score matrix for the alignment.
- `seqA::String`: first input sequence.
- `seqB::String`: second input sequence.
- `alignments::Array{Tuple{String,String},1}`: populate with output alignments.
- `scores::Array{Array{Int,1},1}`: populate with output scores.
- `positions::Array{Array{Tuple{Int,Int},1},1}`: populate with output positions.
- `thresh::Int=0`: score threshold for the alignment.
- `dedup::Bool=true`: whether alignments should be deduplicated.
- `current_alignment::Tuple{String,String}=("","")`: sequences of the current alignment.
- `current_scores::Array{Int,1}=Array{Int,1}()`: scores at each position along the current alignment.
- `current_positions::Array{Tuple{Int,Int},1}=Array{Tuple{Int,Int},1}()`: coordinates of
    each position in the current alignment.
- `maxpathscore::Union{Int,Bool}=false`: the maximum score along the traceback path,
    starts at `false`.
"""
function traceback(align_pos::CartesianIndex{2}, trace_matrix::Array{Array{Int,1},2},
                   score_matrix::Array{Int,2}, seqA::String, seqB::String,
                   alignments::Array{Tuple{String,String},1}, scores::Array{Array{Int,1},1},
                   positions::Array{Array{Tuple{Int,Int},1},1};
                   thresh::Int=0, dedup::Bool=true,
                   current_alignment::Tuple{String,String}=("",""),
                   current_scores::Array{Int,1}=Array{Int,1}(),
                   current_positions::Array{Tuple{Int,Int},1}=Array{Tuple{Int,Int},1}(),
                   maxpathscore::Union{Int,Bool}=false)

    (a, b) = Tuple(align_pos)
    trace = trace_matrix[a,b]
    score = score_matrix[a,b]
    if maxpathscore==false
        maxpathscore=score
    end
    if trace==[1] || trace==[]
        push!(alignments, current_alignment)
        push!(scores, current_scores)
        push!(positions, current_positions)
    else
        if !dedup && maxpathscore-score>=thresh
            push!(alignments, current_alignment)
            push!(scores, [x-score for x in current_scores])
            push!(positions, current_positions)
        end
        if 2 in trace
            traceback(
                (a-1,b-1), trace_matrix, score_matrix, seqA, seqB, alignments, scores, positions,
                thresh=thresh, dedup=dedup, 
                current_alignment=(string(seqA[a-1], current_alignment[1]), string(seqB[b-1], current_alignment[2])),
                current_scores=append!([score], deepcopy(current_scores)),
                current_positions=append!([(a, b)], deepcopy(current_positions)),
                maxpathscore=maxpathscore)
        end
        if 3 in trace
            traceback(
                (a-1,b), trace_matrix, score_matrix, seqA, seqB, alignments, scores, positions,
                thresh=thresh, dedup=dedup, 
                current_alignment=(string(seqA[a-1], current_alignment[1]), string("-", current_alignment[2])),
                current_scores=append!([score], deepcopy(current_scores)),
                current_positions=append!([(a, b)], deepcopy(current_positions)),
                maxpathscore=maxpathscore)
        end
        if 4 in trace
            traceback(
                (a,b-1), trace_matrix, score_matrix, seqA, seqB, alignments, scores, positions,
                thresh=thresh, dedup=dedup, 
                current_alignment=(string("-", current_alignment[1]), string(seqB[b-1], current_alignment[2])),
                current_scores=append!([score], deepcopy(current_scores)),
                current_positions=append!([(a, b)], deepcopy(current_positions)),
                maxpathscore=maxpathscore)
        end
    end
end


"""
    deduplicate(align_pos::CartesianIndex{2}, trace_matrix::Array{Array{Int,1},2},
        score_matrix::Array{Int,2}, thresh::Int, lenA::Int, lenB::Int)
  
Identify whether the alignment end position has any parents in the trace matrix where a
parent is a next alignment position that would contain this position in its alignment.
Return `true` if the alignment position is contained within a parent alignment, and
`false` if not.  
  
---
# Arguments:
- `align_pos::CartesianIndex{2}`: coordinates of the alignment position.
- `trace_matrix::Array{Array{Int,1},2}`: alignment traceback matrix.
- `score_matrix::Array{Int,2}`: score matrix for the alignment.
- `thresh::Int`: alignment score threshold.
- `lenA::Int`: length of seqA + 1 (the max size of the alignment in the a direction).
- `lenB::Int`: length of seqB + 1 (the max size of the alignment in the b direction).
"""
function deduplicate(align_pos::CartesianIndex{2}, trace_matrix::Array{Array{Int,1},2},
                     score_matrix::Array{Int,2}, thresh::Int, lenA::Int, lenB::Int)

    (a, b) = Tuple(align_pos)
    if a==lenA && b==lenB
        return false
    elseif a==lenA
        return ((4 in trace_matrix[a,b+1]) && ((score_matrix[a,b+1]>=thresh) || 
                                              DeDup((a,b+1), trace_matrix, score_matrix, thresh, lenA, lenB)))
    elseif b==lenB
        return ((3 in trace_matrix[a+1,b]) && ((score_matrix[a+1,b]>=thresh) || 
                                              DeDup((a+1,b), trace_matrix, score_matrix, thresh, lenA, lenB)))
    else
        return (((2 in trace_matrix[a+1,b+1]) && ((score_matrix[a+1,b+1]>=thresh) || 
                                              DeDup((a+1,b+1), trace_matrix, score_matrix, thresh, lenA, lenB))) ||
                ((3 in trace_matrix[a+1,b]) && ((score_matrix[a+1,b]>=thresh) || 
                                              DeDup((a+1,b), trace_matrix, score_matrix, thresh, lenA, lenB))) ||
                ((4 in trace_matrix[a,b+1]) && ((score_matrix[a,b+1]>=thresh) || 
                                              DeDup((a,b+1), trace_matrix, score_matrix, thresh, lenA, lenB))))
    end
end


"""
    local_align(seqA::String, seqB::String,
        sub_header::Dict{Char,Int}, sub_matrix::Array{Float64,2})
    local_align(; thresh::Union{Int,Bool}=false, gap_open::Int=-12, gap_extend::Int=-4,
        dedup::Bool=true, dedup_method::String="score")
  
Identify suboptimal alignments of SeqA and SeqB scoring above threshhold. Return the
`score matrix`, `trace matrix`, and an array of alignments, containing `start coordinate`,
`end coordinate`, `score`, and `aligned sequences` for each.  
  
---
# Arguments:
- `seqA::String`: first input sequence.
- `seqB::String`: second input sequence.
- `sub_header::Dict{Char,Int}`: substitution matrix headers.
- `sub_matrix::Array{Float64,2}`: substitution matrix.
- `thresh::Union{Int,Bool}=false`: threshold for the alignment score, `false` triggers an
    automatically-set threshold.
- `gap_open::Int=-12`: penalty for gap opening.
- `gap_extend::Int=-4`: penalty for gap extension.
- `dedup::Bool=true`: whether to filter duplicate alignments (on the same path).
- `dedup_method::String="score"`: method for deduplicating, either longest or highest-scoring
    alignment ("length" or "score").
"""
function local_align(seqA::String, seqB::String,
                     sub_header::Dict{Char,Int}, sub_matrix::Array{Float64,2};
                     thresh::Union{Int,Bool}=false, gap_open::Int=-12, gap_extend::Int=-4,
                     dedup::Bool=true, dedup_method::String="score")
    
    # instantiate scoring matrix and trace matrix
    score_matrix = fill(0,(length(seqA)+1,length(seqB)+1))
    trace_matrix = fill(Array{Int,1}(),(length(seqA)+1,length(seqB)+1))
    
    # fill scoring matrix
    for i in 2:length(seqA)+1
        for j in 2:length(seqB)+1
            # calculate score from each direction
            matchscore = score_matrix[i-1,j-1] + sub_matrix[sub_header[seqA[i-1]],sub_header[seqB[j-1]]]
            hgapscore = score_matrix[i-1,j] + gap_extend + (3 in trace_matrix[i-1,j] ? 0 : gap_open)
            vgapscore = score_matrix[i,j-1] + gap_extend + (4 in trace_matrix[i,j-1] ? 0 : gap_open)
            # find max score and save direction(s)
            score = max(0,matchscore,hgapscore,vgapscore)
            score_matrix[i,j] = score
            trace_matrix[i,j] = findall(x -> x == score, [0,matchscore,hgapscore,vgapscore])
        end
    end
    
    # traceback to identify alignments, optionally deduplicating
    if thresh==false
        # auto-set threshold if not set by the user
        thresh = sort(collect(Iterators.flatten(score_matrix)), rev=true)[10]
        println("Threshold assigned to ", thresh)
    end
    align_ends = findall(x -> x >= thresh, score_matrix)
    if dedup
        align_ends = [
            pos for pos in align_ends 
                if !deduplicate(pos, trace_matrix, score_matrix, thresh, length(seqA)+1, length(seqB)+1)]
    end
    all_alignments = []
    for align_pos in align_ends
        alignments = Array{Tuple{String,String},1}()
        pos_scores = Array{Array{Int,1},1}()
        positions = Array{Array{Tuple{Int,Int},1},1}()
        traceback(align_pos, trace_matrix, score_matrix, seqA, seqB, alignments, pos_scores, positions,
            thresh=thresh, dedup=dedup)
        # select highest-scoring alignment(s) for each alignment path, if de-duplicating by score
        if dedup && dedup_method == "score"
            new_alignments = []
            for (i, scores) in enumerate(pos_scores)
                append!(new_alignments,
                        [[[k-1 for k in positions[i][1]], [k-1 for k in positions[i][x]], scores[x],
                        (alignments[i][1][1:x], alignments[i][2][1:x])]
                        for x in argmax(scores)])
            end
        # otherwise, use the longest alignment for each alignment path
        else
            new_alignments = [
                [[k-1 for k in positions[i][1]], [k-1 for k in positions[i][end]], pos_scores[i][end], x]
                for (i, x) in enumerate(alignments)]
        end
        append!(all_alignments, new_alignments)
    end
    
    # remove any remaining child alignments if deduplicating by score
    if dedup && dedup_method == "score"
        aligned_seqs = [x[4] for x in all_alignments]
        all_alignments = [
            all_alignments[i] for (i,x) in enumerate(aligned_seqs)
            if !(sum([x==(k[1][1:min(length(x[1]),length(k[1]))], k[2][1:min(length(x[1]),length(k[1]))])
                      for k in aligned_seqs]) > 1)]
    end

    return (score_matrix, trace_matrix, all_alignments)
end
