# Data for MMI Failure on Graph States (n = 8)

This directory contains Mathematica `.mx` data files related to monogamy of mutual information (MMI) for graph states on eight vertices. All graphs are considered **up to graph isomorphism**, and entropy vectors are considered **up to subsystem-label symmetry**.

## Directory `MMI-Failure-Stars-n=8`

All files in this subdirectory concern graph states on eight qubits that are (or are locally Clifford–equivalent to) generalized star graphs under some partition.

- `graphRepsSymmetryOrbitsEight.mx`  
  A representative graph (up to isomorphism) for each **unique entropy vector** on eight qubits, modulo subsystem-label symmetry.

- `StarLCorbitsEight.mx`  
  Local-Clifford (LC) orbits of graph states on eight vertices, **excluding** any orbit representatives that are not generalized stars under any partition.

- `uniqueEntropyVectorsEight.mx`  
  All entropy vectors for graph states on eight vertices, listed once per equivalence class **up to subsystem-label symmetry**.

- `triplesEight.mx`  
  The set of all tripartite subsystems on eight qubits that **uniquely define an MMI inequality** (i.e., all triples specifying an MMI on eight-qubit systems).

- `starMMIFailLCorbisEight.mx`  
  LC orbits restricted to **star-graph representatives** of graph states on eight vertices that **fail** at least one MMI.

- `starMMISatisfyLCorbisEight.mx`  
  LC orbits restricted to **star-graph representatives** of graph states on eight vertices that **satisfy all** MMI inequalities.

- `minimalStarFailCIJCount.mx`  
  For each entropy-vector representative corresponding to a generalized star that **fails an MMI of type `MMI_CIJ`**, this file stores entries of the form  
  `{ adjacencyMatrix, failedMMI_CIJInstance, {satisfies, saturates, fails} }`,  
  where `{satisfies, saturates, fails}` is the outcome count for all MMI instances.

- `minimalStarFailOtherMMI.mx`  
  For each entropy-vector representative that **fails some MMI but none of type `MMI_CIJ`**, this file stores entries of the form  
  `{ adjacencyMatrix, failedMMIInstance, starPartition, {satisfies, saturates, fails} }`,  
  where `starPartition` is a partition under which the graph is a generalized star, and `{satisfies, saturates, fails}` is the outcome count for all MMI instances.

- `starRepsSatisfyMMICount.mx`  
  For each entropy-vector representative that **satisfies every MMI**, this file stores entries of the form  
  `{ adjacencyMatrix, starPartition, {satisfies, saturates, fails} }`,  
  where `starPartition` is a partition under which the graph is a generalized star, and `{satisfies, saturates, fails}` is the outcome count for all MMI instances.
