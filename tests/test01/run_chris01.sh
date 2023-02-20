#!/bin/bash
#
#  setup Enzyme, KA and Onan as local copies and on correct branches and with patches
#
\rm -fr .julia
\rm Project.toml
\rm Manifest.toml
export JULIA_DEPOT_PATH=`pwd`/.julia
export JULIA_PATH=/Applications/Julia-1.9.app/Contents/Resources/julia/bin/julia
${JULIA_PATH} --project=. chris01.jl
