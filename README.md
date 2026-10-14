# Monitored free fermions

Simulation code from my MS (Research) thesis at the Indian Institute of Science,
Bengaluru (2025), supervised by Prof. Sumilan Banerjee: response and correlations of a
free-fermion tight-binding chain whose local occupations are continuously and weakly
measured, evolved along quantum trajectories.

More about the project: https://anirudhverma33.github.io/research/monitored-fermions/

This repository is a working copy and is still being organised.

## Contents

| File | What it does |
|---|---|
| `Measurement induced phase transition.ipynb`, `Measurement induced phase transition (1).py` | First trajectory simulations of the monitored chain and its entanglement |
| `initial.ipynb` | Early setup: chain, time step and trajectory evolution |
| `tight.ipynb` | Tight-binding chain dynamics |
| `entropycal.ipynb` | Entanglement entropy from the correlation matrix |
| `entopyhistogram.ipynb` | Distributions of entanglement entropy over trajectories |
| `correlators.ipynb` | Two-time correlation functions |
| `dendencor.ipynb`, `Den_Den_Corr.py` | Density–density correlations |
| `calculation.ipynb`, `dsa.ipynb` | Small scratch calculations |
| `Images/`, `Entanglement_800.png`, `Plots.pptx`, `plots2.pptx` | Figures and plots |

## Requirements

Python 3 with NumPy, SciPy, Matplotlib and tqdm. Some notebooks use `np.complex_`,
which was removed in NumPy 2.0, so use NumPy < 2.0 or replace it with `np.complex128`.
