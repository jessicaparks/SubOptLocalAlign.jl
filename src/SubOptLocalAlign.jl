__precompile__()

"""
# SubOptLocalAlign
Identify and visualize sub-optimal local alignments.
"""
module SubOptLocalAlign

export BLOSUM62

include("read_inputs.jl")
include("local_align.jl")
include("global_align.jl")
include("align_viz.jl")
include("align.jl")

end
