# Global Align
# ============
#
# Identifies a best-scoring global alignment for a pair of sequences.



"""
    global_align(seqA::String, seqB::String,
        sub_header::Dict{Char,Int}, sub_matrix::Array{Int,2},
        gap_open::Int=-12, gap_extend::Int=-4,
        end_gap_open::Int=0, end_gap_extend::Int=0)
  
Global pairwise alignment of the two input sequences. Return `score matrix`, 
`traceback matrix`, and array of `start position`, `end position`, `score`, 
and `alignment sequences`.  
  
Defaults to semi-global alignment, with no penalties for end gaps. 
Note, this does not penalize for gaps on any sequence end, which results in 
very short overlap alignments if the sequences are very dissimilar.  
  
Only returns one best-scoring alignment, regardless of ties for best score.
"""
function global_align(seqA::String, seqB::String,
                      sub_header::Dict{Char,Int}, sub_matrix::Array{Int,2};
                      gap_open::Int=-12, gap_extend::Int=-4,
                      end_gap_open::Int=0, end_gap_extend::Int=0)
    
    m = length(seqA) + 1
    n = length(seqB) + 1
    score_matrix = fill(0,(m,n))
    trace_matrix = fill(0,(m,n))
    
    for i in 2:m
        score_matrix[i,1] = end_gap_open + ((i-1)*end_gap_extend)
        trace_matrix[i,1] = 4
    end
    for j in 2:n
        score_matrix[1,j] = end_gap_open + ((j-1)*end_gap_extend)
        trace_matrix[1,j] = 5
    end

    # fill scoring matrix
    for i in 2:m
        for j in 2:n
            # calculate score from each direction
            matchscore = score_matrix[i-1,j-1] + sub_matrix[sub_header[seqA[i-1]],sub_header[seqB[j-1]]]
            hgapscore = score_matrix[i-1,j] + gap_extend + (2 in trace_matrix[i-1,j] ? 0 : gap_open)
            vgapscore = score_matrix[i,j-1] + gap_extend + (3 in trace_matrix[i,j-1] ? 0 : gap_open)
            # find max score and save direction(s)
            score = max(matchscore,hgapscore,vgapscore)
            score_matrix[i,j] = score
            trace_matrix[i,j] = findall(x -> x == score,
                                        [matchscore,hgapscore,vgapscore])[1]
        end
    end    
    
    max_i = argmax(score_matrix[m,:])
    max_j = argmax(score_matrix[:,n])
    max_ij = argmax([score_matrix[m,max_i], score_matrix[max_j,n]])
    (i, j) = max_ij==1 ? [m,max_i] : [max_j,n]
    align_end = [i, j]
    best_score = score_matrix[i,j]
    
    alignA = ""
    alignB = ""
    while (i > 1) && (j > 1)
        p = trace_matrix[i,j]
        if p == 1
            i = i-1
            j = j-1
            alignA = string(seqA[i], alignA)
            alignB = string(seqB[j], alignB)
        elseif p == 2
            i = i-1
            alignA = string(seqA[i], alignA)
            alignB = string("-", alignB)
        elseif p == 3
            j = j-1
            alignA = string("-", alignA)
            alignB = string(seqB[j], alignB)
        else
            break
        end
    end
    align_start = [i, j]

    return (score_matrix, trace_matrix, [align_start, [x-1 for x in align_end], best_score, (alignA, alignB)])
    
end
