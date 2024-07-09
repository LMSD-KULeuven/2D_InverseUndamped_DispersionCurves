# Inverse approach for dispersion curve calculations

## Paper

This code is discussed in the educational paper entitled
``` 
A guide to numerical dispersion curve calculations: explanation, interpretation and basic Matlab code
```
which is available as a preprint on [arxiv](https://arxiv.org/abs/2311.09843) and [published](https://www.sciencedirect.com/science/article/pii/S0888327024002917) in the journal of Mechanical Systems and Signal Processing. 

```
@article{cool_guide_2024,
	title = {A guide to numerical dispersion curve calculations: {Explanation}, interpretation and basic {Matlab} code},
	volume = {215},
	issn = {0888-3270},
	shorttitle = {A guide to numerical dispersion curve calculations},
	url = {https://www.sciencedirect.com/science/article/pii/S0888327024002917},
	doi = {10.1016/j.ymssp.2024.111393},
	journal = {Mechanical Systems and Signal Processing},
	author = {Cool, Vanessa and Deckers, Elke and Van Belle, Lucas and Claeys, Claus},
	month = jun,
	year = {2024},
}
```

## Abstract

This Matlab script contains an implementation to compute the dispersion curves of 2D periodic media based on a 3D solid finite element unit cell model. The scripting is part of the corresponding educational paper which provides a guide for novel researcher to understand and reason upon dispersion curves.

With all included material, the results and plots of Section 5 in the corresponding paper can be obtained.

## Scripting

- ``src/DispersionCurveCalculation.m``
This is the main Matlab script which has to be run. It will compute the dispersion curves for 2D periodic media with the ‘inverse approach’. A detailed description of the Matlab code is provided in the corresponding paper.

- ``src/FEM_incompatible_modes.m``
This function computes the necessary adaptations to the element stiffness matrix (incompatible mode elements) to resolve the issue of shear locking in linear 3D solid finite elements.

## Use

The [main script (DispersionCurveCalculation.m)](./src/DispersionCurveCalculation.m) should be run in Matlab.

All user-inputs to define the problem are given inside the Matlab script (on [lines 4-25](./src/DispersionCurveCalculation.m#L4)).

<hr/>

All questions and feedback can be addressed to [Claus Claeys](mailto:claus.claeys@kuleuven.be).

We hope our tool can help researchers to better explore and understand dispersion curves during the analysis of periodic media.
