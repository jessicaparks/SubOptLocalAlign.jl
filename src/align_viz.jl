# Align Viz
# =========
#
# Visualizations for pairwise local alignments.



using PlotUtils
using Plots.PlotMeasures
using Plots
using PlotlyJS


"""
    print_in_color(text::String, rgbcolor::RGB{Float64})
    print_in_color(; colorbkgd::Bool=true, colortext::Bool=false, bold::Bool=false)
  
Print RGB-colored background or text.
"""
function print_in_color(text::String, rgbcolor::RGB{Float64};
                        colorbkgd::Bool=true, colortext::Bool=false, bold::Bool=false)
    r = Int(round(red(rgbcolor)*255))
    g = Int(round(green(rgbcolor)*255))
    b = Int(round(blue(rgbcolor)*255))
    background = colorbkgd ? "\e[48;2;$r;$g;$b;249m" : ""
    bold = bold ? "\e[1m" : ""
    textcolor = colortext ? "\e[38;2;$r;$g;$b;249m" : ""
    reset = "\e[0m"
    print(background, bold, textcolor, text, reset)
end


"""
    local_align_print(seqA_id::String, seqB_id::String,
        alignments::Array{Array{Any,1},1}, sub_header::Dict{Char,Int},
        sub_matrix::Array{Float64,2})
  
Print the local alignments between two sequences in colored text.

---
# Arguments:
- `seqA_id::String`: identifier for first input sequence.
- `seqB_id::String`: identifier for second input sequence.
- `alignments::Array{Array{Any,1},1}`: array of alignments; containing `start coordinate`, `end coordinate`,
    `score`, and `alignment sequences` for each
- `sub_header::Dict{Char,Int}`: substitution matrix headers.
- `sub_matrix::Array{Float64,2}`: substitution matrix.
"""
function local_align_print(seqA_id::String, seqB_id::String,
                           alignments::Array{Array{Any,1},1},
                           sub_header::Dict{Char,Int}, sub_matrix::Array{Float64,2})

    # sort alignments by start position and then by score
    alignments = sort(sort(alignments, by=x->x[3], rev=true), by=x->x[1], rev=false)

    # colors to use in print-out
    C(g::PlotUtils.ColorGradient) = RGB[g[i] for i in 3:28]
    colors = reverse(PlotUtils.cgrad(:curl) |> C)

    # lengths of id spacings
    idlength = max(length(seqA_id), length(seqB_id)) + 3
    seqA_spaces = join([" " for i in 1:idlength-length(seqA_id)])
    seqB_spaces = join([" " for i in 1:idlength-length(seqB_id)])

    for (i, (a_start, a_end, a_score, a)) in enumerate(alignments)
        # length of the alignment
        alength = length(a[1])

        # alignment position colors
        c = [if (a[1][j]=='-' || a[2][j]=='-') colors[1]
             else colors[Integer(1 * (14 + sub_matrix[sub_header[a[1][j]],sub_header[a[2][j]]]))] end
             for j in 1:alength]
        
        # print the alignment
        printstyled("Alignment $(i):\n"; bold=true, color=0)
        printstyled("score: $(a_score)\n"; color=0)
        printstyled("coord: $(a_start) to $(a_end)\n\n"; color=0)
        for i in 1:80-idlength:alength
            printstyled(seqA_id, seqA_spaces)
            for j in i:min(i+80-idlength, alength)
                print_in_color(string(a[1][j]), c[j])
            end
            printstyled("\n", seqB_id, seqB_spaces)
            for j in i:min(i+80-idlength, alength)
                print_in_color(string(a[2][j]), c[j])
            end
            print("\n\n")
        end
    end
end


"""
    i_local_align_viz(seqA::String, seqB::String,
        seqA_id::String, seqB_id::String,
        alignments::Array{Array{Any,1},1},
        sub_header::Dict{Char,Int}, sub_matrix::Array{Float64,2};
        global_alignment::Union{Array{Any,1},Bool}=false,
        figurewidth::Int=1000)
  
Interactive visualization of the local alignments between two sequences.

---
# Arguments:
- `seqA::String`: first input sequence.
- `seqB::String`: second input sequence.
- `seqA_id::String`: identifier for first input sequence.
- `seqB_id::String`: identifier for second input sequence.
- `alignments::Array{Array{Any,1},1}`: array of possible alignments; containing
    `start coordinate`, `end coordinate`, `score`, and `alignment sequences` for each.
- `sub_header::Dict{Char,Int}`: substitution matrix headers.
- `sub_matrix::Array{Float64,2}`: substitution matrix.
- `global_alignment::Union{Array{Any,1},Bool}=false`: a global alignment to be included in
    the plot; if `false`, no global alignment included; otherwise, contains
    `start coordinate`, `end coordinate`, `score`, and `alignment sequences`.
- `figurewidth::Int=1000`: width (in pixels) of the plot.
"""
function i_local_align_viz(seqA::String, seqB::String,
                           seqA_id::String, seqB_id::String,
                           alignments::Array{Array{Any,1},1},
                           sub_header::Dict{Char,Int}, sub_matrix::Array{Float64,2};
                           global_alignment::Union{Array{Any,1},Bool}=false,
                           figurewidth::Int=1000)

    # colors to use in figure
    C(g::PlotUtils.ColorGradient) = RGB[g[i] for i in 3:28]
    colors = reverse(PlotUtils.cgrad(:curl) |> C)
    
    # figure layout
    layout = PlotlyJS.Layout(xaxis_range=[0.5, length(seqA)+0.5],
                    xaxis_side="top",
                    yaxis_range=[length(seqB)+0.5, 0.5],
                    yaxis_scaleanchor="x", yaxis_scaleratio=1,
                    hovermode="closest",
                    width=figurewidth,
                    height=figurewidth*length(seqB)/length(seqA),
                    showlegend=false,
                    xaxis_title=seqA_id,
                    yaxis_title=seqB_id,
                    xaxis_tickvals=[i for i in 1:length(seqA)], xaxis_ticktext=split(seqA,""),
                    yaxis_tickvals=[i for i in 1:length(seqB)], yaxis_ticktext=split(seqB,""),
                    xaxis_showgrid=false, yaxis_showgrid=false
    )
    
    # compile traces of lines and points
    traces = []
    for (i, (a_start, a_end, a_score, a)) in enumerate(alignments)
        # coordinates for trace
        x = [a_start[1]+x-1 for x in cumsum([c != '-' for c in a[1]])]
        y = [a_start[2]+x-1 for x in cumsum([c != '-' for c in a[2]])]
        # marker sizes
        s = [if (a[1][j]=='-' || a[2][j]=='-') 8
             else abs(sub_matrix[sub_header[a[1][j]],sub_header[a[2][j]]])+1 end
             for j in 1:length(a[1])]
        # marker colors
        c = [if (a[1][j]=='-' || a[2][j]=='-') colors[1]
             else colors[Integer(1 * (14 + sub_matrix[sub_header[a[1][j]],sub_header[a[2][j]]]))] end
             for j in 1:length(a[1])]
        # text to display when hovering over points
        t = [if (a[1][j]=='-' || a[2][j]=='-') "gap"
             else sub_matrix[sub_header[a[1][j]],sub_header[a[2][j]]] end
             for j in 1:length(a[1])]
        d = [string(a[1][j],a[2][j]) for j in 1:length(a[1])]
        hm = "<b>aligned:</b> %{customdata}<br><b>score:</b> %{hovertext} / $a_score<br><i>(alignment $i)</i>"
        # line trace
        l_trace = PlotlyJS.scatter(x=x, y=y, mode="lines", name="", line_color="gray", line_width=5,
                                   opacity=0.75, hoverinfo="none")
        # marker trace
        m_trace = PlotlyJS.scatter(x=x, y=y, mode="markers", name="", hovertext=t, customdata=d, hovertemplate=hm,
                                   marker_size=s, marker_color=c, marker_symbol="hexagon", marker_opacity=1,
                                   marker_sizemode="area", marker_sizeref=0.05)
        push!(traces, l_trace)
        push!(traces, m_trace)
    end
    
    if global_alignment!=false
        (a_start, a_end, a_score, a) = global_alignment
        # coordinates for trace
        x = [a_start[1]+x-1 for x in cumsum([c != '-' for c in a[1]])]
        y = [a_start[2]+x-1 for x in cumsum([c != '-' for c in a[2]])]
        # line trace
        l_trace = PlotlyJS.scatter(x=x, y=y, mode="lines", name="", line_color="black", line_width=2,
                                   opacity=0.75, hoverinfo="none")
        push!(traces, l_trace)
    end

    # create plot
    plt = PlotlyJS.plot([x for x in traces], layout)
    return plt
end


"""
    s_local_align_viz(seqA::String, seqB::String,
        seqA_id::String, seqB_id::String,
        alignments::Array{Array{Any,1},1},
        sub_header::Dict{Char,Int}, sub_matrix::Array{Float64,2};
        global_alignment::Union{Array{Any,1},Bool}=false,
        figurewidth::Int=1000)
  
Static visualization of the local alignments between two sequences.
  
---
# Arguments:
- `seqA::String`: first input sequence.
- `seqB::String`: second input sequence.
- `seqA_id::String`: identifier for first input sequence.
- `seqB_id::String`: identifier for second input sequence.
- `alignments::Array{Array{Any,1},1}`: array of possible alignments; containing
    `start coordinate`, `end coordinate`, `score`, and `alignment sequences` for each.
- `sub_header::Dict{Char,Int}`: substitution matrix headers.
- `sub_matrix::Array{Float64,2}`: substitution matrix.
- `global_alignment::Union{Array{Any,1},Bool}=false`: a global alignment to be included in
    the plot; if `false`, no global alignment included; otherwise, contains
    `start coordinate`, `end coordinate`, `score`, and `alignment sequences`.
- `figurewidth::Int=1000`: width (in pixels) of the plot.
"""
function s_local_align_viz(seqA::String, seqB::String,
                           seqA_id::String, seqB_id::String,
                           alignments::Array{Array{Any,1},1},
                           sub_header::Dict{Char,Int}, sub_matrix::Array{Float64,2};
                           global_alignment::Union{Array{Any,1},Bool}=false,
                           figurewidth::Int=1000)

    # colors to use in figure
    C(g::PlotUtils.ColorGradient) = RGB[g[i] for i in 3:28]
    colors = reverse(PlotUtils.cgrad(:curl) |> C)
    
    # figure layout
    if length(seqA)<50
        xtickdata = ([i for i in 1:length(seqA)], split(seqA,""))
        ytickdata = ([i for i in 1:length(seqB)], split(seqB,""))
    else
        xtickdata = cat([1], [i for i in 50:50:length(seqA)], [length(seqA)], dims=1)
        ytickdata = cat([1], [i for i in 50:50:length(seqB)], [length(seqB)], dims=1)
    end
    plt = Plots.plot(legend=:none, xlabel=seqA_id, ylabel=seqB_id,
        xlims=(0.5,length(seqA)+0.5), ylims=(0.5,length(seqB)+0.5), yflip=true,
        size=(figurewidth+20, (length(seqB)*figurewidth/length(seqA))+20), xmirror=true,
        xticks=xtickdata,
        yticks=ytickdata,
        foreground_color_axis=:transparent,
        yrotation=90,
        left_margin=20px,
        bottom_margin=20px)
    
    # add traces of lines and points
    for (a_start, a_end, a_score, a) in alignments
        # coordinates for trace
        x = [a_start[1]+x-1 for x in cumsum([c != '-' for c in a[1]])]
        y = [a_start[2]+x-1 for x in cumsum([c != '-' for c in a[2]])]
        # marker sizes
        s = [if (a[1][j]=='-' || a[2][j]=='-') 8
             else abs(sub_matrix[sub_header[a[1][j]],sub_header[a[2][j]]])+1 end
             for j in 1:length(a[1])]
        # marker colors
        c = [if (a[1][j]=='-' || a[2][j]=='-') colors[1]
             else colors[Integer(1 * (14 + sub_matrix[sub_header[a[1][j]],sub_header[a[2][j]]]))] end
             for j in 1:length(a[1])]
        # plot line
        Plots.plot!(plt, x, y, seriestype=:line,
                    alpha=0.5, linecolor="gray", linewidth=4)
        # plot points
        Plots.plot!(plt, x, y, seriestype=:scatter,
                    alpha=1, markercolor=c,
                    markerstrokewidth = 0, markershape=:hexagon, markersize=s)
    end
    if global_alignment!=false
        (a_start, a_end, a_score, a) = global_alignment
        # coordinates for trace
        x = [a_start[1]+x-1 for x in cumsum([c != '-' for c in a[1]])]
        y = [a_start[2]+x-1 for x in cumsum([c != '-' for c in a[2]])]
        # plot line
        Plots.plot!(plt, x, y, seriestype=:line,
                    alpha=0.75, linecolor="black", linewidth=2)
    end
    return plt
end
