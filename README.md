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
  
To view the interactive graph in Jupyter, the following are also required:  
(run from the terminal and refresh any open JupyterLab pages afterward)  
```bash
jupyter labextension install @jupyterlab/plotly-extension
jupyter labextension install @webio/jupyter-lab-provider
```
Check currently installed JupyterLab extensions by:
```bash
jupyter labextension list
```

## Usage

The main method of SubOptLocalAlign is `align`. See the function documentation for a complete
list of arguments.

```julia
using SubOptLocalAlign

plot = SubOptLocalAlign.align("seqA.fasta", "seqB.fasta");
```

To use without interactive graphing, which will produce a static figure:

```julia
using SubOptLocalAlign

plot = SubOptLocalAlign.align("seqA.fasta", "seqB.fasta", figure_type="static");
```

Other methods can be used in cases where only subsets of the functionality are needed. See
function documentation for description and arguments.
