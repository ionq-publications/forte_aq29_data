# Supplemental Data for "Benchmarking a trapped-ion quantum computer with 30 qubits"

These data are used to reproduce the results in the [paper](https://arxiv.org/abs/2308.05071), "Benchmarking a trapped-ion quantum computer with 30 qubits". 

## License
The content in this repository is copyright IonQ Inc., all rights reserved. It is licensed under the [CC BY-NC 4.0 license](https://creativecommons.org/licenses/by-nc/4.0/legalcode.en).

## Circuits
This directory contains all of the application-oriented algorithms used in the paper: 
- `AE` - Amplitude Estimation
- `HSIM` - Hamiltonian Simulation
- `MC` - Monte Carlo Sampling
- `PE` - Phase Estimation
- `QFT` - Quantum Fourier Transform
- `VQE` - VQE Simulation

All circuits are provided in [OpenQasm 2.0](https://github.com/openqasm/openqasm/tree/OpenQASM2.x) as submitted to the system (input folder), as transpiled by qiskit (qiskit folder), and as transpiled by the IonQ compiler to the level of native gates (ionq folder). The naming convention is `ALG`N`QQ`V`II`_`V` where
- `ALG` - algorithm identifier as listed above
- `QQ` - number of qubits used in the circuit
- `II` - instance identifier as a number
- `V` - variant identifier (specified for the circuits compiled to native gates)

For example, `circuits/input/AEN04V02.qasm` is the second instance of the amplitude estimation algorithm on four qubits as submitted to the system, while `circuits/ionq/AEN04V02_2` is its second variant transpilation to native gates.

## Figures
This directory contains jupyter notebooks to reproduce figures shown in the paper from the raw data. To run these notebooks, you will need a basic scientific python setup, including at least:
- numpy
- matplotlib
- pandas
- seaborn

In addition, we include code in to perform efficient, approximate plurality voting, as explained in the paper, in `pluralityvoite.py`. 

## Results
This folder contains the raw experimental and simulation results as json files. These files are read in and processed using the analysis routines in `figures`. We also include cached error mitigated fidelities for convenience, since plurality voting takes 5-10 minutes to regenerate these.
