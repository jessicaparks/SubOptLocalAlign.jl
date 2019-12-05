# SubOptLocalAlign

## Description

SubOptLocalAlign identifies sub-optimal local alignments for a pair of sequences and provides
visualizations for these alignments, including a highlighted print-out of the alignments and
static and interactive graphs of the alignments.

## Installation

Install from the Julia REPL or from Jupyter:
```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/jessicaparks/SubOptLocalAlign.jl"))
```
This assumes you already have [Julia](https://julialang.org/downloads/) installed. Output will display best in an interactive
environment, such as [Jupyter Lab](https://jupyterlab.readthedocs.io/en/stable/getting_started/installation.html),
which also requires [IJulia](https://github.com/JuliaLang/IJulia.jl).

## Usage

The main method of SubOptLocalAlign is `align`.

```julia
using SubOptLocalAlign

plot = SubOptLocalAlign.align("seqA.fasta", "seqB.fasta");
```

Other methods can be used in cases where only subsets of the functionality are needed. See
function documentation for description and arguments.
