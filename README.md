# MMI Failure for Graph States Data Repository

## Introduction
Data files accompanying the paper *MMI Failure for Graph States* (arXiv:XXXX.XXXXX) by Jesús Fuentes, Cynthia Keeler, Jason Pollack, and William Munizzi. This repository consolidates every artifact needed to reproduce the monogamy-of-mutual-information (MMI) study on stabilizer and graph states.

## Contents
1. `Data/MMI-Failure-Stars-n=8`
   - Complete data on the MMI violation of graph states at \(n = 8\), including canonical representatives and evaluation logs.
2. `Notebooks/MMIFailureGraphsState.nb`
   - A comprehensive functionality notebook that instructs and demonstrates the functions most relevant to the project, alongside additional exploratory tools.
   - Capabilities showcased in the notebook span state generation, tableaux representation creation and manipulation, entropy vector computation in tableau, graph, and foliage representations, entropy-vector symmetrization, unique-orbit extraction, entropy inequality verification, and Mathematica Graph object to TikZ conversion (with table generation and styling utilities).
3. `Package/MMIFailureGraphStates.wl`
   - A ready-to-deploy Mathematica package bundling all function definitions, usage messages, and helper utilities highlighted in the notebook.

## Usage Notes
- Open `Notebooks/MMIFailureGraphsState.nb` in Mathematica to follow guided examples and reproduce the main simulations from the paper.
- Load the package in a Mathematica session with:
  ```wl
  Get[FileNameJoin[{NotebookDirectory[], "MMIFailureGraphStates.wl"}]];
  ```
  After loading, evaluate `?SymbolName` (for example `?GenerateStabilizerStates`) to view usage instructions for public functions.

## Citation
If you use this repository, please cite the associated paper (*MMI Failure for Graph States*, arXiv:XXXX.XXXXX) and acknowledge the authors Jesús Fuentes, Cynthia Keeler, Jason Pollack, and William Munizzi. See `NOTICE`.
