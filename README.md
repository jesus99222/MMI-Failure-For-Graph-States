# Monogamy of Mutual Information in Graph States Data Repository

## Introduction
Data files accompanying the paper *Monogamy of Mutual Information in Graph States* ([arXiv:2511.19585](http://arxiv.org/abs/2511.19585)) by Jesus Fuentes, Cynthia Keeler, William Munizzi, and Jason Pollack. This repository consolidates every function needed to reproduce the monogamy of mutual information (MMI) study on stabilizer and graph states.

## Contents
1. `Data/MMI-Failure-Stars-n=8`
   - Complete data on the MMI violation of graph states at \(n = 8\), including canonical representatives and evaluation logs.
2. `Notebooks/MMIFailureGraphsState.nb`
   - A comprehensive functionality notebook that instructs and demonstrates the functions most relevant to the project, alongside additional exploratory tools.
   - Capabilities showcased in the notebook span state generation, tableaux representation creation and manipulation, entropy vector computation in tableau, graph, and foliage representations, entropy-vector symmetrization, unique-orbit extraction, entropy inequality verification, and Mathematica Graph object to TikZ conversion (with table generation and styling utilities).
3. `Package/MMI-Failure-For-Graph-States.wl`
   - A ready-to-deploy Mathematica package bundling all function definitions, usage messages, and helper utilities highlighted in the notebook.

## Usage Notes
- Open `Notebooks/MMIFailureGraphsState.nb` in Mathematica to follow guided examples and reproduce the main simulations from the paper.
- Load the package in a Mathematica session with:
  ```wl
  Get[FileNameJoin[{NotebookDirectory[], "MMI-Failure-For-Graph-States.wl"}]];
  ```
  After loading, evaluate `?SymbolName` (for example `?GenerateStabilizerStates`) to view usage instructions for public functions.

## Entropy Vector HTML Library
- Pre-rendered entropy tables and graph previews live in `entropy/`. Open `entropy/index.html` to browse the global catalog or `entropy/<table-id>/index.html` for a specific family (for example, `entropy/minimalStarNonTFailSvecs/index.html`).
- Each card links to a per-graph page containing the reduced entropy vector, a copy-to-clipboard shortcut, and a thumbnail of the underlying graph (all thumbnails are stored centrally in `entropy/images/`).

## Citation
If you use this repository, please cite the associated paper (*Monogamy of Mutual Information in Graph States*, [arXiv:2511.19585](http://arxiv.org/abs/2511.19585)
) and acknowledge the authors, Jesus Fuentes, Cynthia Keeler, William Munizzi, and Jason Pollack. See [NOTICE](NOTICE).
