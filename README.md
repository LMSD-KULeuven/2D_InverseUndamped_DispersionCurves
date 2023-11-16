# Inverse approach for dispersion curve calculations

## Paper

```
Will be updated.
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
