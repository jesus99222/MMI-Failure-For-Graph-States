# Data Loading and Reproduction Guide

This folder contains Mathematica `.mx` dumps of symbols used for the n = 8 analysis in our paper. Each file was saved so that `Get["Data/<name>.mx"]` defines the symbol whose name matches the filename (e.g., `graphRepsSymmetryOrbitsEight.mx` defines `graphRepsSymmetryOrbitsEight`).

Note: Mathematica `.mx` files are version- and architecture-dependent. If `Get` fails due to version mismatch, regenerate in your environment or re-dump the data.

## Quick Start: Load All Data

From the repository root:

```wl
(* Load all .mx datasets in this directory *)
Scan[Get, FileNames["Data/*.mx"]]

(* Alternatively, from a notebook inside the repo: *)
dataDir = FileNameJoin[{NotebookDirectory[], "Data"}];
Scan[Get, FileNames["*.mx", dataDir]]

(* Load a single data instance *)
Get["Data/graphRepsSymmetryOrbitsEight.mx"]
```

After loading, you can directly use the symbols listed below.

## Symbol Inventory (what to call and what it contains)

- `graphRepsSymmetryOrbitsEight`
  - Graph representatives for each unique entropy vector at `n = 8`.
  - Constructed by: `generateSymmetryOrbitGraphRepresentatives[8, 182]`.
  - Type: list of 8-vertex graph adjacency representations (SparseArray objects usable with `AdjacencyGraph`).

- `starRepsSymmetryOrbitsEight`
  - Minimal-edge star-graph representatives of `graphRepsSymmetryOrbitsEight`.
  - Constructed by: `findStarGraphReps[graphRepsSymmetryOrbitsEight, 8]`.
  - Type: list of graph adjacencies (SparseArray objects usable with `AdjacencyGraph`).

- `uniqueEntropyVectorsEight`
  - Entropy vectors for `graphRepsSymmetryOrbitsEight`.
  - Constructed by: `ComputeEntropyVectorsAdj[graphRepsSymmetryOrbitsEight]`.
  - Type: list of reduced entropy vectors.

- `starRepsFailMMI`
  - Subset of `starRepsSymmetryOrbitsEight` that fail MMI for some disjoint triples.
  - Constructed via boolean mask `evalsBool` from `MMIChecker` on `uniqueEntropyVectorsEight` and `triplesEight = DisjointTriples[8]`:
    `Pick[starRepsSymmetryOrbitsEight, evalsBool, False]`.

- `starRepsSatisfyMMI`
  - Subset of `starRepsSymmetryOrbitsEightthat` that never fail any MMI.
  - Constructed by: `Pick[starRepsSymmetryOrbitsEight, evalsBool, True]`.

- `minimalStarSatisfyMMINontrivialInt`
  - For each graph in `starRepsSatisfyMMI`, a minimal-edge representative with a nontrivial intersection for some `(C, I, J)`.
  - Constructed by: `searchNonTrivialIntLC[starRepsSatisfyMMI, triplesEight]`.
  - Type: list of `{graphAdj, instance}` where `instance` is a triple of the form `(C, I, J)` with nontrivial intersection of the column spaces of `CI, CJ,` and `CK`, if it exists.

- `minimalStarFailMMINontrivialInt`
  - Same nontrivial intersection search for `starRepsFailMMI`.
  - Constructed by: `searchNonTrivialIntLC[starRepsFailMMI, triplesEight]`.

- `minimalStarFailMMINontrivialIntAbsent`
  - Those entries of `minimalStarFailMMINontrivialInt` with `instance == {}` (no nontrivial intersection).
  - Constructed by: `Select[minimalStarFailMMINontrivialInt, (#[[2]] == {}) &]`.

- `minimalStarFailMMINontrivialIntPresent`
  - Those entries with non-empty `instance` (nontrivial intersection exists).
  - Constructed by: `Complement[minimalStarFailMMINontrivialInt, minimalStarFailMMINontrivialIntAbsent]`.

- `minimalStarNonTFailCIJ`
  - For each graph in `minimalStarFailMMINontrivialIntPresent`, a failing `(C, I, J)` with maximal support size.
  - Constructed using `InstancesCIJ` and `MMICheckerInputSelectTriples`, selecting the failing triple with maximal `Length[Flatten[...]]`.
  - Type: list of `{graphAdj, triple}`.

- `minimalStarNoCIJFail`
  - For each graph in `minimalStarFailMMINontrivialIntAbsent`, a general `(I,J,K)` failing triple chosen from `triplesEight` with maximal support size.
  - Type: list of `{graphAdj, triple}`.

- `starRepsSatisfyMMICount`
  - `minimalStarSatisfyMMINontrivialInt` paired with counts `{Satisfies, Saturates, Fails}` of MMI outcomes across `triplesEight`, sorted by satisfies count and then by edge count.

- `minimalStarNonTFailCIJCount`
  - `minimalStarNonTFailCIJ` with counts `{Satisfies, Saturates, Fails}`, sorted by fails and then by edge count.

- `minimalStarTFailnoCIJCount`
  - `minimalStarNoCIJFail` with counts `{Satisfies, Saturates, Fails}`, sorted by fails and then by edge count.

All of the above symbols are available as `.mx` files in this directory and can be loaded with `Get`.

## Basic Usage Examples

```wl
(* Visualize a few graphs *)
Map[AdjacencyGraph, starRepsSymmetryOrbitsEight[[;; 5]]]

(* Inspect entropy vector size *)
Length[uniqueEntropyVectorsEight[[1]]]

(* Access a chosen (C, I, J) triple from the nontrivial-fail set *)
minimalStarNonTFailCIJ[[1, 2]]
```

## Reproducing Tables 10–12 (Appendix A.3)

Prerequisites: ensure the following functions are available in your session (e.g., load your project package/notebook):

```
generateSymmetryOrbitGraphRepresentatives, findStarGraphReps,
ComputeEntropyVectorsAdj, DisjointTriples, MMIChecker,
searchNonTrivialIntLC, InstancesCIJ, MMICheckerInputSelectTriples
```

1) Graph representatives and star reps at n = 8

```wl
graphRepsSymmetryOrbitsEight = generateSymmetryOrbitGraphRepresentatives[8, 182];
starRepsSymmetryOrbitsEight = findStarGraphReps[graphRepsSymmetryOrbitsEight, 8];
```

2) Identify MMI failing and satisfying cases

```wl
uniqueEntropyVectorsEight = ComputeEntropyVectorsAdj[graphRepsSymmetryOrbitsEight];
triplesEight = DisjointTriples[8];

evalsBool = ParallelMap[MMIChecker[#, triplesEight] &, uniqueEntropyVectorsEight];

starRepsFailMMI = Pick[starRepsSymmetryOrbitsEight, evalsBool, False];
starRepsSatisfyMMI = Pick[starRepsSymmetryOrbitsEight, evalsBool, True];
```

3) Minimal-edge star reps with nontrivial intersections

```wl
minimalStarSatisfyMMINontrivialInt = searchNonTrivialIntLC[starRepsSatisfyMMI, triplesEight];

minimalStarFailMMINontrivialInt = searchNonTrivialIntLC[starRepsFailMMI, triplesEight];

minimalStarFailMMINontrivialIntAbsent = Select[minimalStarFailMMINontrivialInt, (#[[2]] == {}) &];
minimalStarFailMMINontrivialIntPresent = Complement[minimalStarFailMMINontrivialInt, minimalStarFailMMINontrivialIntAbsent];
```

4) Table graphs (visual)

```wl
Print["Table 10 graphs"]; Map[AdjacencyGraph, minimalStarSatisfyMMINontrivialInt[[All, 1]]]
Print["Table 11 graphs"]; Map[AdjacencyGraph, minimalStarFailMMINontrivialIntPresent[[All, 1]]]
Print["Table 12 graphs"]; Map[AdjacencyGraph, minimalStarFailMMINontrivialIntAbsent[[All, 1]]]
```

5) Select representative failing triples

```wl
minimalStarNonTFailCIJ = Module[{evals, failed, lenghtsFailed},
  ParallelTable[
    triplesCenter = InstancesCIJ[starFailCIJ, triplesEight, 8];
    evals = MMICheckerInputSelectTriples[(ComputeEntropyVectorsAdj[{starFailCIJ}])[[1]], triplesCenter];
    failed = Pick[triplesCenter, evals, "Fails"]; 
    lenghtsFailed = Map[Length[Flatten[#]] &, failed];
    {starFailCIJ, First[Pick[failed, lenghtsFailed, Max[lenghtsFailed]]]},
    {starFailCIJ, (minimalStarFailMMINontrivialIntPresent[[All, 1]])}
  ]
];

minimalStarNoCIJFail = Module[{evals, failed, lenghtsFailed},
  ParallelTable[
    evals = MMICheckerInputSelectTriples[ComputeEntropyVectorsAdj[{starFailnoCIJ}][[1]], triplesEight];
    failed = Pick[triplesEight, evals, "Fails"]; 
    lenghtsFailed = Map[Length[Flatten[#]] &, failed];
    {starFailnoCIJ, First[Pick[failed, lenghtsFailed, Max[lenghtsFailed]]]},
    {starFailnoCIJ, minimalStarFailMMINontrivialIntAbsent[[All, 1]]}
  ]
];
```

6) Count outcomes and sort to match table ordering

```wl
satisfySvecs = ComputeEntropyVectorsAdj[minimalStarSatisfyMMINontrivialInt[[All, 1]]];
evalsSatisfy = ParallelMap[MMICheckerInputSelectTriples[#, triplesEight] &, satisfySvecs];
evalsSatisfyCount = ParallelMap[{{Count[#, "Satisfies"], Count[#, "Saturates"], Count[#, "Fails"]}} &, evalsSatisfy];
minimalStarSatisfyMMINontrivialIntCount = Join[minimalStarSatisfyMMINontrivialInt, evalsSatisfyCount, 2 ];
starRepsSatisfyMMICount = SortBy[minimalStarSatisfyMMINontrivialIntCount, {#[[3, 1]] &, Total@Flatten@#[[1]] &}];

FailsCIJSvecs = ComputeEntropyVectorsAdj[minimalStarNonTFailCIJ[[All, 1]]];
evalsFailsCIJ = ParallelMap[MMICheckerInputSelectTriples[#, triplesEight] &, FailsCIJSvecs];
evalsFailsCIJcount = ParallelMap[{{Count[#, "Satisfies"], Count[#, "Saturates"], Count[#, "Fails"]}} &, evalsFailsCIJ];
minimalStarNonTFailCIJCount = Join[minimalStarNonTFailCIJ, evalsFailsCIJcount, 2 ];
minimalStarNonTFailCIJCount = SortBy[minimalStarNonTFailCIJCount, {#[[3, 3]] &, Total@Flatten@#[[1]] &}];

FailsnoCIJSvecs = ComputeEntropyVectorsAdj[minimalStarNoCIJFail[[All, 1]]];
evalsFailsnoCIJ = ParallelMap[MMICheckerInputSelectTriples[#, triplesEight] &, FailsnoCIJSvecs];
evalsFailsnoCIJcount = ParallelMap[{{Count[#, "Satisfies"], Count[#, "Saturates"], Count[#, "Fails"]}} &, evalsFailsnoCIJ];
minimalStarTFailnoCIJCount = Join[minimalStarNoCIJFail, evalsFailsnoCIJcount, 2 ];
minimalStarTFailnoCIJCount = SortBy[minimalStarTFailnoCIJCount, {#[[3, 3]] &, Total@Flatten@#[[1]] &}];
```

### Table mapping
- Table 10 → `minimalStarSatisfyMMINontrivialInt`
- Table 11 → `minimalStarFailMMINontrivialIntPresent`
- Table 12 → `minimalStarFailMMINontrivialIntAbsent`

## Notes
- `ParallelMap` can be replaced with `Map` if no parallel kernels are available.
- The counts lists have the order `{Satisfies, Saturates, Fails}`.
- Graph adjacencies can be visualized via `AdjacencyGraph[...]`.

