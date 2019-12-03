# Local Align
# ===========
#
# Identifies sub-optimal local alignments for a pair of sequences.



"""
traceback(align_pos, trace_matrix, score_matrix, seqA, seqB, alignments, scores, positions,
          thresh, dedup, current_alignment, current_scores, current_positions, maxpathscore)
traceback to identify alignments

Arguments:
align_pos: coordinates of the alignment end position
trace_matrix: alignment traceback matrix
score_matrix: score matrix for the alignment
seqA: sequence A
seqB: sequence B
alignments: empty list to fill with resulting alignments
scores: empty list to fill with resulting scores
positions: empty list to fill with resulting positions
current_alignment: sequences in the current alignment at the position
current_scores: scores at each position along the current alignment
current_positions: coordinates of each position in the current alignment

Returns:
pushes updates to the alignments, scores, and positions arrays
"""
function traceback(align_pos, trace_matrix, score_matrix, seqA, seqB, alignments, scores, positions;
                   thresh=0, dedup=true,
                   current_alignment=("",""), current_scores=[], current_positions=[], maxpathscore=false)
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
            TraceBack(
                (a-1,b-1), trace_matrix, score_matrix, seqA, seqB, alignments, scores, positions,
                thresh=thresh, dedup=dedup, 
                current_alignment=(string(seqA[a-1], current_alignment[1]), string(seqB[b-1], current_alignment[2])),
                current_scores=append!([score], deepcopy(current_scores)),
                current_positions=append!([(a, b)], deepcopy(current_positions)),
                maxpathscore=maxpathscore)
        end
        if 3 in trace
            TraceBack(
                (a-1,b), trace_matrix, score_matrix, seqA, seqB, alignments, scores, positions,
                thresh=thresh, dedup=dedup, 
                current_alignment=(string(seqA[a-1], current_alignment[1]), string("-", current_alignment[2])),
                current_scores=append!([score], deepcopy(current_scores)),
                current_positions=append!([(a, b)], deepcopy(current_positions)),
                maxpathscore=maxpathscore)
        end
        if 4 in trace
            TraceBack(
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
deduplicate(align_pos, trace_matrix, score_matrix, thresh, lenA, lenB)
identify whether the alignment end position has any parents in the trace matrix
where a parent is a next alignment position that would contain this position in its alignment

Arguments:
align_pos: coordinates of the alignment position
trace_matrix: alignment traceback matrix
score_matrix: score matrix for the alignment
thresh: alignment threshold
lenA: length of seqA + 1 (the max size of the alignment in the a direction)
lenB: length of seqB + 1 (the max size of the alignment in the b direction)

Returns:
(bool): true if the alignment position is contained within a parent alignment, false if not;
"""
function deduplicate(align_pos, trace_matrix, score_matrix, thresh, lenA, lenB)
    (a, b) = Tuple(align_pos)
    if a==lenA && b==lenB
        return false
    elseif a==lenA
        # could remove the recursion because only adding gaps
        return ((4 in trace_matrix[a,b+1]) && ((score_matrix[a,b+1]>=thresh) || 
                                              DeDup((a,b+1), trace_matrix, score_matrix, thresh, lenA, lenB)))
    elseif b==lenB
        # could remove the recursion because only adding gaps
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
local_align(seqA, seqB, sub_header, sub_matrix, thresh, gap_open, gap_extend, dedup, dedup_method)
identify all suboptimal alignments of SeqA and SeqB scoring above threshhold

Arguments:
seqA: first input sequence
seqB: second input sequence
thresh (int): threshold for the alignment score, defaults to false (which triggers an automatic threshold)
sub_header (list): substitution matrix headers
sub_matrix (matrix): substitution matrix
gap_open (int): contribution to alignment score for each gap opening, default = -12
gap_extend (int): contribution to alignment score for each gap extension, default = -4
dedup (bool): filter "duplicate" alignments (on the same path) if true, do not filter if false
dedup_method (str): method for filtering, either longest or highest-scoring alignment ("length" or "score")

Returns:
score matrix
traceback matrix
(start position, end position, score, alignment) list for all suboptimal alignments
"""
function local_align(seqA, seqB, sub_header, sub_matrix;
                    thresh=false, gap_open=-12, gap_extend=-4, dedup=true, dedup_method="score")
    
    # instantiate scoring matrix and trace matrix
    score_matrix = fill(0,(length(seqA)+1,length(seqB)+1))
    trace_matrix = fill([],(length(seqA)+1,length(seqB)+1))
    
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
        alignments = []
        pos_scores = []
        positions = []
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
