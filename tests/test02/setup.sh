#!/bin/bash
#
#  setup Enzyme, KA and Onan as local copies and on correct branches and with patches
#
\rm -fr .julia
\rm Project.toml
\rm Manifest.toml
export JULIA_DEPOT_PATH=`pwd`/.julia
export JULIA_PATH=/Applications/Julia-1.9.app/Contents/Resources/julia/bin/julia
${JULIA_PATH} --project=. -e 'using Pkg;Pkg.add("Plots");Pkg.add("LaTeXStrings")'
${JULIA_PATH} --project=. -e 'using Pkg;Pkg.add(name="EnzymeCore",rev="main")'
${JULIA_PATH} --project=. -e 'using Pkg;Pkg.add(name="Enzyme",rev="main")'
### ${JULIA_PATH} --project=. -e 'using Pkg;Pkg.add(url="https://github.com/christophernhill/KernelAbstractions.jl",rev="vc/rules_cnh_fix_Project")'
### ${JULIA_PATH} --project=. -e 'using Pkg;Pkg.add(url="https://github.com/christophernhill/OceanLES.jl",rev="cnh/onan_ka_enz")'
${JULIA_PATH} --project=. -e 'using Pkg;Pkg.add("Oceananigans")'
exit
#
# ${JULIA_PATH} --project=. chris02.jl
#
