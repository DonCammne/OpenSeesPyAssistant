# OpenSeesPyAssistant

The OpenSeesPyAssistant (OSPA) is a Python-based simulation tools library for 2D reinforced concrete, steel and composite structures. It is meant to assist the user modeling nonlinear structures with the interpreter OpenSeesPy. Various features are implemented (material models, members model, fibers section, etc) to help the user to create flexible, reliable, systematic and readable main code programs for nonlinear modeling

- [MIT License](https://choosealicense.com/licenses/mit/)
- [GitHub Pages](https://github.com/DonCammne/OpenSeesPyAssistant)
- [Online Documentaion](https://ospa.karma-riuk.com/namespaces.html)
- Author: Carmine Schipani

## Requirements

- [openseespy==3.3.0.1.1](https://pypi.org/project/openseespy/3.3.0.1.1/) (newer version removed functions used in the library)
- [numpy](https://pypi.org/project/numpy/)

## Features

- Member model:
    - Elastic element
    - Force-based element
    - Spring-based element with zero-length element
    - Gradient-inelastic flexibilita-based element
    - Panel zone member
- Material model:
    - Uniaxial bilinear
    - GMP 1970
    - UVC and VC
    - Modified IMK
    - Gupta 1999
    - Skiadopoulos 2021
    - Mander 1988 (confined and unconfined)
- Fiber section:
    - Rectangular reinforced concrete
    - Circluar reinforced concrete
    - I shape steel profile
- Analysis options (with automatic convergence analysis):
    - Gravity (vertical loading)
    - Lateral Force
    - Pushover
    - Loading Protocol
- Automatic units management
- Plotting functions (fibers, memebrs, material models)
- ID convention management
- Import/export of data from analysis to post-processing
- Discretizer for curves
- Quick geometry and frame model templates

For the post processing, three Matlab module are presented in the folder [MATLAB_postprocessing](https://github.com/DonCammne/OpenSeesPyAssistant/tree/main/MATLAB_postprocessing) (GitHub).


## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install OSPA using the following prompt command.

```bash
pip install OpenSeesPyAssistant
```

## Usage

The entire library is imported with this lines:

```python
from OpenSeesPyAssistant.Section import *
from OpenSeesPyAssistant.DataManagement import *
from OpenSeesPyAssistant.ErrorHandling import *
from OpenSeesPyAssistant.Units import *
from OpenSeesPyAssistant.Constants import *
from OpenSeesPyAssistant.Fibers import *
from OpenSeesPyAssistant.Connections import *
from OpenSeesPyAssistant.FunctionalFeatures import *
from OpenSeesPyAssistant.MemberModel import *
from OpenSeesPyAssistant.AnalysisAndPostProcessing import *
from OpenSeesPyAssistant.GeometryTemplate import *
```


## User Manual

An auto generated User Manual with Doxygen is available in the folder [user_manual](https://github.com/DonCammne/OpenSeesPyAssistant/tree/main/user_manual) (GitHub).


## Examples

An application of the OSPA library can be found in the folder [examples](https://github.com/DonCammne/OpenSeesPyAssistant/tree/main/examples) (GitHub).
The specimen studied is UT04 from [Shin 2017](https://repositories.lib.utexas.edu/handle/2152/47306) and an example of postprocessing with the Matlab module proposed is also availble [here](https://github.com/DonCammne/OpenSeesPyAssistant/blob/main/examples/UT04_LP/AnalysisResults.m).

## Library Status

The library is currently work in progress. Future implementaions:

- ZeroLength Sections
- Bond SP01
- Mass and dynamic analysis
- RBS
- Splacing
- Leaning columns
- Slab interaction (with the steel beam)