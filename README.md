# NHSWM2DCurved

This MATLAB/GNU Octave code repository can be used to simulate weakly non-hydrostatic shallow water flow in closed bodies of water. Visualization of simulation results is also facilitated by these scripts. See below for details.

## Getting Started

    $ git clone https://github.com/dsteinmo/NHSWM2DCurved.git
    $ cd NHSWM2DCurved
    $ matlab -nodesktop -nosplash          # OR: octave-cli

## Theoretical description

See the publication referenced [below](#citing-this-work) for theoretical background. For a more detailed look, see [chapter 4 from Derek Steinmoeller's Ph.D. Thesis (University of Waterloo, 2014)](https://raw.githubusercontent.com/dsteinmo/NHSWM2DCurved/refs/heads/master/docs/Steinmoeller_Derek_ch_4.pdf).

## Visualization

FIXME: Add details.

## Citing this work

Results from using the code herein was published in 'Ocean Modelling' in 2016. BibTeX:
```
@article{steinmoeller2016discontinuous,
  title={Discontinuous Galerkin methods for dispersive shallow water models in closed basins: Spurious eddies and their removal using curved boundary methods},
  author={Steinmoeller, DT and Stastna, M and Lamb, KG},
  journal={Ocean Modelling},
  volume={107},
  pages={112--124},
  year={2016},
  publisher={Elsevier}
}
```