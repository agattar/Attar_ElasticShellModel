This repository includes simulation and analyses codes for the Attar et al, 2023.


Nucleus_Generator_F.py -> Generates the data file containing the shell model for the nuclear lamina and the polymer blocks representing the chromatin. To generate the polymer blocks, it reads 'mm9_domains_200k.npy' and uses the Hi-C compartment sequences.

BondingMin.py -> It uses the data file generated from the minimization simulation and reassigns bonds for the simulation efficiency.

Dumpread_FS.py -> Reads the dump file generated using LAMMPS MD simulations and exports the average radial distance, density, and tensors.

FFT_steady.py -> Reads the dump file and exports a CSV file that contains the amplitudes and the corresponding wavenumbers.

Input files folder contaıns the input files required for running the simulations.

References:

Attar, A. G., Paturej, J., Banigan, E. J., & Erbas, A. (2023). Polymer modeling suggests correlations between chromatin phase separation and nuclear shape fluctuations. bioRxiv. https://doi.org/10.1101/2023.12.16.571697
