An Enzyme + Oceananigans setup based on 

  1. simple 1d diffusion equation
  2. no kernel abstractions
  3. enzyme 0.11.xxxx

This setup produces a adjoint sensitivity plot
that mostly matches what we expect analytically.
There is a small mystery at the periodic overlaps,
but otherwise the solution fits the mathematics.

This setup provides a useful reference to compare
with enzyme 0.11 and beyond modifications for supporting
kernel abstractions.
