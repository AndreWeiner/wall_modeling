# Development of wall functions for turbulent flows

The repository is structured as follows:
- *test_cases*: OpenFOAM simulation setups
- *notebooks*: Jupyter notebooks for post-processing of results and other visualizations

## Getting started

Compiling the wall function and running the test cases requires an installation of [OpenFOAM-v2012](https://openfoam.com/). Other releases might work as well but have not been tested.

To run a test case, create a *run* folder (ignored by version control), copy the case from *test_cases* to *run*, and execute the *Allrun* script.

```
mkdir -p run
cp -r test_cases/turbulentFlatPlate run/
cd run/turbulentFlatPlate
./Allrun
```

## Singularity

TODO
- add singularity container to make results reproducible

## References

The basic ideas of this work were first developed in the context of mass transfer at rising bubbles. Instead of the turbulent viscosity, the molecular diffusivity of a chemical species is corrected to mitigate the so-called high-Schmidt number problem. The profile function's parameter is determined by the species mass in an interface cell.
```
@article{weiner2017,
title = {Advanced subgrid-scale modeling for convection-dominated species transport at fluid interfaces with application to mass transfer from rising bubbles},
journal = {Journal of Computational Physics},
volume = {347},
pages = {261-289},
year = {2017},
issn = {0021-9991},
doi = {https://doi.org/10.1016/j.jcp.2017.06.040},
author = {Andre Weiner and Dieter Bothe}
}
```

In the following reference, the ideas above were first implemented in OpenFOAM and applied to surfactant adsorption rising bubbles.
```
@article{pesci2018,
title={Computational analysis of single rising bubbles influenced by soluble surfactant},
volume={856},
DOI={10.1017/jfm.2018.723},
journal={Journal of Fluid Mechanics},
publisher={Cambridge University Press},
author={Pesci, Chiara and Weiner, Andre and Marschall, Holger and Bothe, Dieter},
year={2018},
pages={709–763}
}
```

Instead of an analytical profile function, one can also use a machine learning model to correct diffusive and convective fluxes.
```
@article{weiner2019,
author = {Weiner, Andre and Hillenbrand, Dennis and Marschall, Holger and Bothe, Dieter},
title = {Data-Driven Subgrid-Scale Modeling for Convection-Dominated Concentration Boundary Layers},
journal = {Chemical Engineering \& Technology},
volume = {42},
number = {7},
pages = {1349-1356},
keywords = {Data-driven modeling, High Schmidt numbers, Machine learning, Reactive mass transfer, Subgrid-scale modeling},
doi = {https://doi.org/10.1002/ceat.201900044},
year = {2019}
}
```

The following two PhD theses contain more detailed descriptions as well as more extensive parameter studies and applications than the references above. *pesci2019* focuses more on surfactants and mass transfer while *weiner2020* addresses the modeling principles and extensions to reactive mass transfer by means of machine learning.
```
@phdthesis{weiner2020,
          school = {Technical University of Darmstadt, Mathematical Modeling and Analysis},
         address = {Darmstadt},
           month = {February},
           title = {Modeling and simulation of convection-dominated species transfer at rising bubbles},
          author = {Andre Weiner},
            year = {2020},
             url = {http://tuprints.ulb.tu-darmstadt.de/11405/}
}
@phdthesis{pesci2019,
            year = {2019},
          author = {Chiara Pesci},
         address = {Darmstadt},
           title = {Computational Analysis of Fluid Interfaces Influenced by Soluble Surfactant},
           month = {July},
          school = {Technische Universit{\"a}t},
             url = {http://tuprints.ulb.tu-darmstadt.de/9303/}
}
```

## Get in touch

The easiest and quickest way to get in touch is to open an issue in this repository.