(* ::Package:: *)

BeginPackage["MMIFailureGraphStatesPackage`"];

  ClearAll[
      GenerateStabilizerStates,
      tableauToStabilizerGenerators,
      ComputeEntropyVectors,
      ComputeEntropyVectorsAdj,
      EquivalentEntropyVectors,
      UniqueOrbits,
      PermuteEntropyVector,
      DisjointTriples,
      MMIChecker,
      MMICheckerInputSelectTriples,
      MMICheckerEntropyInput,
      generateAllSimpleUndirectedGraphs,
      areLocalComplements,
      generateSymmetryOrbitGraphRepresentatives,
      allGeneralizedGraphsnoInnerNonTrivialInterInfo,
      SampleAllGeneralizedGraphsnoInnerNonTrivialInterInfoTriples,
      findStarGraphReps,
      EntropyVectorsFoliageRep,
      fourPartitions,
      applyPermutations,
      canSubsystems,
      validCenters,
      findMinimalEdgeRep,
      findBiggestCenter,
      componentDiameters,
      pathTieBreaker,
      searchNonTrivialIntLC,
      nontrivialInterAdjs,
      InstancesCIJ,
      TikZTableFromAdjTriples
  ];

 GenerateStabilizerStates::usage =
    "GenerateStabilizerStates[n] returns a list of n-qubit stabilizer tableaux (global phases ignored), one representative per stabilizer state.";

 tableauToStabilizerGenerators::usage =
    "tableauToStabilizerGenerators[tab] converts a stabilizer tableau tab to its list of Pauli-string generators (e.g. {{\"X\",\"I\",...},...}).";

 ComputeEntropyVectors::usage =
    "ComputeEntropyVectors[tableaux] computes the reduced-entropy vector for each stabilizer tableau in tableaux.";

 ComputeEntropyVectorsAdj::usage =
    "ComputeEntropyVectorsAdj[adjs] computes the reduced-entropy vector(s) for graph state(s) defined by adjacency matrix/matrices adjs. ComputeEntropyVectorsAdj[{AdjacencyMatrix@graph}] handles a single graph; a list of matrices yields a list of vectors.";

 EquivalentEntropyVectors::usage =
    "EquivalentEntropyVectors[n] returns {baseVector, eqVectors} encoding the action of qubit permutations on n-qubit entropy vectors for use with PermuteEntropyVector and UniqueOrbits.";

 UniqueOrbits::usage =
    "UniqueOrbits[vectors, baseVector, eqVectors] returns one canonical representative per symmetry orbit of the given entropy vectors under the permutation action defined by baseVector and eqVectors.";

 PermuteEntropyVector::usage =
    "PermuteEntropyVector[vec, baseVector, eqVectors] returns the full symmetry orbit of vec under the permutation group encoded by eqVectors (with baseVector used for indexing).";

 DisjointTriples::usage =
    "DisjointTriples[n] returns the list of triples {I,J,K} of pairwise-disjoint, nonempty subsets of Range[n], each defining an MMI instance.";

 MMIChecker::usage =
    "MMIChecker[entropyVector, triples] returns True if the entropyVector satisfies/saturates every MMI instance in triples, and False if any instance fails.";

 MMICheckerInputSelectTriples::usage =
    "MMICheckerInputSelectTriples[entropyVector, triples] returns a list of 'Satisfies' | 'Saturates' | 'Fails' labels, one per triple in triples.";

 MMICheckerEntropyInput::usage =
    "MMICheckerEntropyInput[entropyVector] returns 'Satisfies' | 'Saturates' | 'Fails' for every MMI instance on entropyVector.";

 generateAllSimpleUndirectedGraphs::usage =
    "generateAllSimpleUndirectedGraphs[n] returns all non-isomorphic simple undirected graphs on n vertices as Graph objects.";

 areLocalComplements::usage =
    "areLocalComplements[adj1, adj2] returns True if the adjacency matrices adj1 and adj2 are LC-equivalent (reachable via a sequence of local complementations), and False otherwise.";

 generateSymmetryOrbitGraphRepresentatives::usage =
    "generateSymmetryOrbitGraphRepresentatives[n, orbitCount] uses random sampling to produce a list of Graphs, one representative per entropy-vector symmetry orbit, until orbitCount distinct orbits are found. Set SeedRandom for reproducibility.";

 allGeneralizedGraphsnoInnerNonTrivialInterInfo::usage =
    "allGeneralizedGraphsnoInnerNonTrivialInterInfo[n] returns a list of {adjacency, triples} pairs for the generalized graphs (Sec. 4.3) on n, ignoring intra-region adjacencies, that exhibit at least one non-trivial {C, I, J} intersection.";

 SampleAllGeneralizedGraphsnoInnerNonTrivialInterInfoTriples::usage =
    "SampleAllGeneralizedGraphsnoInnerNonTrivialInterInfoTriples[n, samplesPerCenter] randomly samples generalized graphs (Sec. 4.3) on n and returns {adjacency, triples} pairs for those with at least one non-trivial intersection.";

 findStarGraphReps::usage =
    "findStarGraphReps[adjs, n] filters or extracts generalized star family representatives (in the sense of Sec. 4.3) from the adjacency matrix of graphs on n vertices (e.g., generalized star/starlike graphs).";

 fourPartitions::usage =
    "fourPartitions[n] returns all set partitions of Range[n] into four nonempty, disjoint blocks.";

 applyPermutations::usage =
    "applyPermutations[lists, vertices, perms] applies each permutation of vertices to every subsystem specification in lists, canonicalizing and de-duplicating the results.";

 canSubsystems::usage =
    "canSubsystems[list] sorts vertices within each block and sorts the blocks themselves, removing duplicate subsystem specifications.";

 validCenters::usage =
    "validCenters[adj, fourParts] returns all {C,{I,J,K}} choices from the supplied 4-partitions where C is a valid center whose complement blocks are pairwise nonadjacent; returns {\"NF\"} if none exist.";

 findMinimalEdgeRep::usage =
    "findMinimalEdgeRep[lcOrbitsWithCenters] selects, within each LC orbit, the representative with the fewest edges together with its center data.";

 findBiggestCenter::usage =
    "findBiggestCenter[{adj, centers}] picks the center partition with the largest |C|, returning {adj, {C, I, J}}.";

 componentDiameters::usage =
    "componentDiameters[adj] returns the diameter of each connected component of the graph with adjacency adj.";

 pathTieBreaker::usage =
    "pathTieBreaker[adj] returns {maxDiameter, countMax, secondMax, countSecond} for sorting graphs by component diameters.";

 searchNonTrivialIntLC::usage =
    "searchNonTrivialIntLC[adjs, triples] scans each graph's LC orbit to find a minimal-edge LC-equivalent representative that exhibits at least one non-trivial intersection from triples. Returns a list of {rep, triplesFound} pairs (triplesFound may be {}).";

 InstancesCIJ::usage =
    "InstancesCIJ[adj, triples, n] returns every triple of the form {C,I,J}, as defined in Section 4.2 for the graph state represented by adj.";

 TikZTableFromAdjTriples::usage =
    "TikZTableFromAdjTriples[data, opts] generates a TikZ table summarizing graphs, their associated {C, I, J} instances and evaluation labels. Options include Columns, RowVspace, MultiPage, Version, CanvasCM, PanelNumbering, PanelStyle, EvalLabelSize, EvalLabelOffsetCM, and ColorScheme.";


  Begin["`Private`"];

  (* ========================================= *)
  (* Tableaux Utilities                        *)
  (* ========================================= *)

  (* generateStabilizerGroup computes the complete stabilizer group by iterative multiplication of
  generators. *)
  generateStabilizerGroup[generators_,parallelQ_]:=Block[{g = Join[{IdentityMatrix[First[Dimensions[First[generators]]]]},generators],count=0,newElements=1},
	count=Length[g];
	
	If[parallelQ==True,
	While[newElements>0,
	g = DeleteDuplicates[Flatten[ParallelTable[g[[i]] . g[[j]],{i,Length[g]},{j,Length[g]}],1]];
	newElements = Length[g]-count;
	count = Length[g];
	],
	While[newElements>0,
	g = DeleteDuplicates[Flatten[Table[g[[i]] . g[[j]],{i,Length[g]},{j,Length[g]}],1]];
	newElements = Length[g]-count;
	count = Length[g];
	]];
 g];

  (* pauliStringToOperator converts symbolic Pauli strings to tensor product matrices. *)
  pauliStringToOperator[pauliString_] :=
    Apply[KroneckerProduct]@ReplaceAll[
      pauliString,
      {
        "I" -> PauliMatrix[0],
        "X" -> PauliMatrix[1],
        "Y" -> PauliMatrix[2],
        "Z" -> PauliMatrix[3],
        "-I" -> -PauliMatrix[0],
        "-X" -> -PauliMatrix[1],
        "-Y" -> -PauliMatrix[2],
        "-Z" -> -PauliMatrix[3],
        "iI" -> I PauliMatrix[0],
        "iX" -> I PauliMatrix[1],
        "iY" -> I PauliMatrix[2],
        "iZ" -> I PauliMatrix[3],
        "-iI" -> -I PauliMatrix[0],
        "-iX" -> -I PauliMatrix[1],
        "-iY" -> -I PauliMatrix[2],
        "-iZ" -> -I PauliMatrix[3]
      }
    ];

  (* stabilizerGeneratorsToState builds the stabilizer-state projector from generators. *)
  stabilizerGeneratorsToState [generatingSet_]:=Block[{stabGroup = generateStabilizerGroup[pauliStringToOperator/@generatingSet,False]},(1/2^Length[generatingSet])*Sum[stabGroup[[i]],{i,Length[stabGroup]}]]

  (* tableauToStabilizerGenerators converts binary tableaux to signed Pauli strings. *)
  tableauToStabilizerGenerators[tableau_]:=Block[{dim = Length[tableau],generatorList={}, tableauLastPhase},
	tableauLastPhase=ConstantArray[0, Dimensions[tableau]];
	If[Mod[Dimensions[tableau][[2]], 2]==0, tableauLastPhase=tableau,tableauLastPhase[[All, 1]] = tableau[[All, -1]] ; tableauLastPhase[[All,2;; ]] = tableau[[All, ;;-2]] ];
	If[OddQ[Length[tableauLastPhase[[1]]]],Table[AppendTo[generatorList,Table[{tableauLastPhase[[i,j]],tableauLastPhase[[i,j+dim]]},{j,2,dim+1}]],{i,dim}];
	Table[PrependTo[generatorList[[k]],tableauLastPhase[[k,1]]],{k,Length[generatorList]}],
	Table[AppendTo[generatorList,Table[{tableauLastPhase[[i,j]],tableauLastPhase[[i,j+dim]]},{j,dim}]],{i,dim}]];
	ReplaceAll[ReplaceAll[generatorList,{{0,0}->"I",{1,0}->"X",{1,1}->"Y",{0,1}->"Z"}],{0->"+1",1->"-1"}]
 ];

  (* stabilizerGeneratorsToTableau maps Pauli-string generators back to tableaux. *)
  stabilizerGeneratorsToTableau[generatingSet_]:=Block[{dim = Length[First[generatingSet]],negatives = First/@Join[Position[generatingSet,"-I"],Position[generatingSet,"-X"],Position[generatingSet,"-Y"],Position[generatingSet,"-Z"]],preTab =Table[Join[{0},generatingSet[[i]],generatingSet[[i]]],{i,Length[generatingSet]}],trTab },
	Table[preTab[[negatives[[k]],1]]=BitXor[preTab[[negatives[[k]],1]],1],{k,Length[negatives]}];
	preTab = ReplaceAll[preTab,{"-I"->"I","-X"->"X","-Y"->"Y","-Z"->"Z"}];
	trTab = Transpose[preTab];
	Table[trTab[[i]] =ReplaceAll[trTab[[i]],{"X"->1,"Y"->1,"I"->0,"Z"->0}] ,{i,2,dim+1}];
	Table[trTab[[j]] =ReplaceAll[trTab[[j]],{"Z"->1,"Y"->1,"X"->0,"I"->0}] ,{j,dim+1,2*dim+1}];
	Transpose[trTab][[All,2;;]]
 ];

  (* ========================================= *)
  (* Clifford Operators on Tableaux            *)
  (* ========================================= *)

  phasedCNOT[tableau_, control_ , target_, n_]:= Module[{newState = tableau},
	newState[[All, target]] = BitXor[tableau[[All, control]], tableau[[All, target]]  ];
	newState[[All, n + control]] = BitXor[tableau[[All, n +control]], tableau[[All, n +target]]  ];
	newState[[All, 2n + 1]] =  BitXor[ tableau[[All, 2n +1]], BitAnd[BitAnd[tableau[[All, control]], tableau[[All, n+target]]],BitXor[tableau[[All, target]],tableau[[All, n+control]],ConstantArray[1,n]]]];
	newState
 ];

  phasedHadamard[tableau_, qubit_, n_]:= Module[{newState=tableau},
	newState[[All, qubit]] = tableau[[All, n +qubit]];
	newState[[All, n+qubit]] = tableau[[All, qubit]];
	newState[[All, 2n + 1]] = BitXor[tableau[[All, 2n+1]], BitAnd[tableau[[All, qubit]], tableau[[All, n +qubit]]]];
	newState
 ];

  phasedPhase[tableau_, qubit_, n_]:= Module[{newState=tableau},
	newState[[All, n +qubit]] =BitXor[tableau[[All, n +qubit]], tableau[[All, qubit]]];
	newState[[All, 2n + 1]] = BitXor[tableau[[All, 2n+1]], BitAnd[tableau[[All, qubit]], tableau[[All, n +qubit]]]];
	newState
 ];

  CNOT[tableau_, control_ , target_, n_]:= Module[{newState = tableau},
	newState[[All, target]] = BitXor[tableau[[All, control]], tableau[[All, target]]  ];
	newState[[All, n + control]] = BitXor[tableau[[All, n +control]], tableau[[All, n +target]]  ];
	
	newState
 ];

  Hadamard[tableau_,qubit_,n_]:=Module[{newState=tableau},If[qubit==0,Return[newState];];
	newState[[All,qubit]]=tableau[[All,n+qubit]];
	newState[[All,n+qubit]]=tableau[[All,qubit]];
	newState
	];

  Phase[tableau_,qubit_,n_]:=Module[{newState=tableau},If[qubit==0,Return[newState];];
newState[[All,n+qubit]]=BitXor[tableau[[All,n+qubit]],tableau[[All,qubit]]];
newState];

  CZGate[tableau_, control_ , target_, n_]:= Module[{newState = tableau},
	newState[[All, n+target]] = BitXor[tableau[[All, n+target]], tableau[[All, control]]  ];
	newState[[All, n + control]] = BitXor[tableau[[All, target]], tableau[[All, n +control]]  ];
	newState
 ];

  CliffordGenerators[Dimension_]:=Block[{ntuples=Tuples[Range[Dimension],{2}],cnots, pandh},
	cnots=DeleteCases[Table[If[ntuples[[j,1]]!=ntuples[[j,2]],ntuples[[j]],Null],{j,Length[ntuples]}],Null];
	pandh=Range[1,Dimension];
	{cnots, pandh}
 ];

  (* ========================================= *)
  (* Stabilizer Group Operations               *)
  (* ========================================= *)

 g[x1_,z1_,x2_,z2_]:=Module[{out},
	If[x1==0&&z1==0,out=0,
	If[x1==1&&z1==1,out=z2-x2,
	If[x1==1&&z1==0,out=z2 ( 2 x2-1),
	out=x2 (1-2 z2)]]];
 out];

  rowsum[tableau_,h_,i_,n_]:=Module[{newRow,x1,z1,x2,z2,sum,parity},
	x1=tableau[[h,1;;n]];
	z1=tableau[[h,n+1;;2 n]];
	x2=tableau[[i,1;;n]];
	z2=tableau[[i,n+1;;2 n]];
	
	newRow=BitXor[tableau[[h,All]],tableau[[i,All]]];
	
	sum=Total[MapThread[g,{x1,z1,x2,z2}]];
	parity=Mod[2 Last[tableau[[h]]]+2 Last[tableau[[i]]]+sum,4];
	If[parity==0,newRow[[-1]]=0,newRow[[-1]]=1];
	newRow
 ];

  SymplecticMatrix[n_]:=ArrayFlatten[{{ConstantArray[0,{n,n}],IdentityMatrix[n]},{IdentityMatrix[n],ConstantArray[0,{n,n}]}}];

  isSymplectic[tableau_]:=Module[{n=Length[tableau],sMat,zeroMat},
	zeroMat=ConstantArray[0,{n,n}];
	sMat=SymplecticMatrix[n];
	Mod[Dot[tableau,sMat,Transpose[tableau]], 2]==zeroMat
 ];

  (* ========================================= *)
  (* Stabilizer State Generation               *)
  (* ========================================= *)

  GenerateVacuumTableau[n_]:=Table[Join[ConstantArray[0,n],ReplacePart[ConstantArray[0,n],i->1]],{i,1,n}]

  GeneratePlusTableau[n_]:= ArrayFlatten[{{IdentityMatrix[n], ConstantArray[0,{n,n}]}}];

  GHZTableau[n_]:=Module[{tableau,row1,rows},
	row1=Join[Table[1,{n}],Table[0,{n}]];
	rows=Table[Join[Table[0,{n}],
	Table[If[i==j||i==j+1,1,0],{i,1,n}] ],
	{j,1,n-1}];
	tableau=Prepend[rows,row1];
	Return[tableau]
 ];

  GenerateStabilizerStates[n_]:=Module[{gates=CliffordGenerators[n],stabStates,newStates,targetSize,parallelNewStates},
	stabStates={GenerateVacuumTableau[n]};
	targetSize=Product[2^(n-k)+1,{k,0,n-1}];
	newStates=stabStates;
	While[Length[stabStates]<targetSize,
	parallelNewStates=Join[
	Flatten[Table[CNOT[stab,gate[[1]],gate[[2]],n],{stab,newStates},{gate,gates[[1]]}],1],Flatten[Table[Hadamard[stab,qubit,n],{stab,newStates},{qubit,gates[[2]]}],1],
	Flatten[Table[Phase[stab,qubit,n],{stab,newStates},{qubit,gates[[2]]}],1]];
	newStates=Complement[RowReduce[#,Modulus->2]&/@parallelNewStates,stabStates];
	stabStates=Union[stabStates,newStates];
	;];
 stabStates];

  StabEvolution[tableaux_, gates_]:=Module[{n, newStates},
	n=Length[tableaux[[1]]];
	newStates=Join[
	If[gates[[1]]!={}, Flatten[Table[CNOT[stab,gate[[1]],gate[[2]],n],{stab,tableaux},{gate,gates[[1]]}],1], {}],If[gates[[2]]!={},Flatten[Table[Hadamard[stab,qubit,n],{stab,tableaux},{qubit,gates[[2]]}],1], {}],
	If[gates[[3]]!={}, Flatten[Table[Phase[stab,qubit,n],{stab,tableaux},{qubit,gates[[3]]}],1], {}]];
 newStates];

  EvolUnentangleGHZ[tableaux_]:=Module[{n, cnotsArgs, newStates},
	n=Length[tableaux[[1]]];
	cnotsArgs=Table[{1,i,n}, {i, 2, n}];
	newStates=Table[
	FoldList[
	Function[{state,op},If[op==="Hadamard",Hadamard[state,1,n],(*Apply Hadamard after all CNOTs*)CNOT[state,op[[1]],op[[2]],op[[3]]]]],tableau,cnotsArgs],
	{tableau,tableaux}];
	newStates=Flatten[newStates, {1,2}];
	newStates=RowReduce[#,Modulus->2]&/@newStates;
	DeleteDuplicates[newStates]
 ];

  EvolEntangleVacuum[tableaux_]:=Module[{n, cnotsArgs, newStates},
	n=Length[tableaux[[1]]];
	cnotsArgs=Table[{1,i,n}, {i, 2, n}];
	newStates=Table[FoldList[Function[{state,op},CNOT[state,op[[1]],op[[2]],op[[3]]] ],Hadamard[tableau,1,n],cnotsArgs],{tableau,tableaux}];
	newStates=Flatten[newStates, {1,2}];
	newStates=RowReduce[#,Modulus->2]&/@newStates;
	DeleteDuplicates[newStates]
 ];

  ApplyHPCircuit[circuit_List, stab_] := Module[{newStab = stab, n = Length[stab]},
    Do[
      Switch[circuit[[qubit]],
        1, newStab = Hadamard[newStab, qubit, n],
        2, newStab = Phase[newStab, qubit, n],
        _, Null
      ],
      {qubit, n}
    ];
    newStab
  ];

  applyAllHPs[tableau_, allHPs_List] := ApplyHPCircuit[#, tableau] & /@ allHPs;

  applyCZcircuit[tableau_, circuit_List] := Module[{newTableau = tableau, n = Length[tableau]},
    Do[
      newTableau = CZGate[newTableau, gate[[1]], gate[[2]], n],
      {gate, circuit}
    ];
    newTableau
  ];

  LCorbit[tableau_]:=Module[{allHPs, n =Length[tableau], newStates, localSphere},
	allHPs=Tuples[{0,1,2},n];
	localSphere={tableau};
	newStates=localSphere;
	While[Length[newStates]!=0,
	newStates=Join@@Map[applyAllHPs[#,allHPs]&,newStates];
	newStates=Complement[RowReduce[#,Modulus->2]&/@newStates,localSphere];
	localSphere=Union[localSphere,newStates];
	;];
 localSphere];

  (* ========================================= *)
  (* Graph States                              *)
  (* ========================================= *)

  generateGraphStates[n_]:=Module[{pairs, circuits, plusState, graphStates},
	pairs = Tuples[Range[1,n], 2] ;
	pairs =DeleteDuplicates[Select[Map[Sort,pairs],#[[1]]!=#[[2]]&]];
	
	circuits = Subsets[pairs, {1, Length[pairs]}];
	
	plusState = GeneratePlusTableau[n];
	
	graphStates =ParallelMap[applyCZcircuit[plusState, #]&, circuits];
	graphStates= Join[{plusState}, graphStates];
	
	graphStates
 ];


  degreeNGStates[tableaux_, deg_Integer?NonNegative] := Module[{adjs, degrees, hasdegreeN,
  positions},
    adjs = graphTabtoAdj /@ tableaux;
    degrees = vertexDegrees /@ adjs;
    hasdegreeN = If[MemberQ[#, deg], 1, 0] & /@ degrees;
    positions = Flatten[Position[hasdegreeN, 1]];
    tableaux[[positions]]
  ];

  graphTabtoAdj[tableau_] := tableau[[All, Length[tableau] + 1 ;;]];

  graphAdjtoTab[adj_] := Join[IdentityMatrix[Length[adj]], adj, 2];

  (* ========================================= *)
  (* Graph Generation & Operations             *)
  (* ========================================= *)

 generateSymmetryOrbitGraphRepresentatives[n_,nOrbits_]:=Module[{vertices=Range[n],possibleEdges, baseVector,graphReps={},orbits={},eqVectors,randomEdges,uniqueSigs,randomGraphs,sigGraphPairs,progressCount=0,progressFrac=0., graphNames, randomGraphsSvecs, keys, randomE , seen=<||>, newKeys, pos, newSeen, newOrbits, randomAdjs, sVecHash, CZs, initializationSize=4, maxCZsLenght=2},

	(*
	sVecHash[S] returns a stable hash for an entropy vector or adjacency signature.
		-S:entropy vector or data structure to hash.
	-SHA256 hash value used for orbit de-duplication.
	*)
	sVecHash[S_]:=Hash[S, "SHA256"];
	
	randomEdges[number_]:=Table[RandomSample[vertices, 2], number];
	
	{baseVector,eqVectors}=EquivalentEntropyVectors[n];
	
	randomAdjs =Table[
	upper=RandomInteger[{0,1},{n,n}];
	adj=UpperTriangularize[upper,1];
	adj=adj+Transpose[adj], {initializationSize}
	];
	
	randomGraphs=Join[{ AdjacencyMatrix@Graph[Range[n], {}],  AdjacencyMatrix@CompleteGraph[n]}, randomAdjs];
	
	Monitor[
	While[
	Length[orbits]!=nOrbits,
	randomGraphsSvecs=ComputeEntropyVectorsAdj[randomGraphs];
	keys =Map[sVecHash,randomGraphsSvecs];
	newKeys=Map[!KeyExistsQ[seen, #]&, keys];
	
	pos =Pick[Range[Length[randomGraphsSvecs]], newKeys];
	
	If[pos !={},
	
	
	sigGraphPairs =orbitSignature[randomGraphs[[pos]], baseVector,eqVectors];
	uniqueSigs=DeleteDuplicates[sigGraphPairs, (#1[[1]]==#2[[1]])&];
	newOrbits=uniqueSigs[[All,1]];
	
	newSeen=Thread[Map[sVecHash[#]&, Flatten[newOrbits,1]]->True];
	AssociateTo[seen,newSeen];
		   orbits= Join[orbits,newOrbits];
	graphReps= Join[graphReps,uniqueSigs[[All,2]]]
	
	];


randomE=randomEdges[RandomInteger[maxCZsLenght]];
CZs=ConstantArray[0, {n,n}];
CZs=ReplacePart[CZs,Thread[randomE->1]];
CZs =Transpose[CZs]+CZs;

randomGraphs=Map[Mod[(CZs+#),2]&, graphReps];
progressCount=Length[graphReps];
progressFrac=progressCount/nOrbits;
],

Column[{Row[{"Orbits found: ",Dynamic[progressCount,TrackedSymbols:>{progressCount}],"/",nOrbits}],
ProgressIndicator[Dynamic[progressFrac,TrackedSymbols:>{progressFrac}],{0,1}]}]];
SparseArray/@graphReps
];

  orbitSignature[adjs_, baseVector_, eqVectors_] := Module[{svecs, neworbits, signature},
    svecs = ComputeEntropyVectorsAdj[adjs];
    neworbits = Sort /@ (PermuteEntropyVector[#, baseVector, eqVectors] & /@ svecs);
    signature = Transpose[{neworbits, adjs}];
    DeleteDuplicatesBy[signature, First]
  ];

  degreeNGraphs[adjs_, deg_Integer?NonNegative] := Module[{degrees, hasdegreeN, positions},
    degrees = vertexDegrees /@ adjs;
    hasdegreeN = If[MemberQ[#, deg], 1, 0] & /@ degrees;
    positions = Flatten[Position[hasdegreeN, 1]];
    adjs[[positions]]
  ];

  adjacencyMatrices[n_Integer?Positive] := Module[{allMatrices},
    allMatrices = Tuples[{0, 1}, {n, n}];
    Select[allMatrices, (Tr[#] == 0 && Transpose[#] == #) &]
  ];

  generateAllSimpleUndirectedGraphs[n_Integer?Positive] := Module[{vertices, possibleEdges, edgeSets,
  allGraphs},
    vertices = Range[n];
    possibleEdges = Subsets[vertices, {2}];
    edgeSets = Subsets[possibleEdges];
    allGraphs = Graph[vertices, #] & /@ edgeSets;
    DeleteDuplicates[allGraphs, IsomorphicGraphQ]
  ];

  vertexDegrees[adj_] := Total[adj];

  canonicalForm[m_] := Transpose[Sort[Transpose[m]]];

  allIsomorphicLabeledGraphs[g_Graph] := Module[{n, labels, perms},
    n = VertexCount[g];
    labels = VertexList[g];
    perms = Permutations[labels];
    DeleteDuplicates[Relabel[g, Thread[labels -> #]] & /@ perms, IsomorphicGraphQ]
  ];

  LocalComplementation[adj_, vertex_] := Module[{neighbors, neighborAdj, rowC, localComplement},
    rowC = Normal[adj[[vertex]]];
    neighbors = Flatten[Position[rowC, 1]];
    If[neighbors =!= {},
      neighborAdj = adj[[neighbors, neighbors]];
      localComplement = adj;
      localComplement[[neighbors, neighbors]] = 1 - IdentityMatrix[Length[neighbors]] - neighborAdj;
      localComplement,
      adj
    ]
  ];

  CanonicalAdj[adj_] := AdjacencyMatrix@CanonicalGraph@AdjacencyGraph@adj;

 LCorbit[adj_]:=Module[{newGraphs,allLCs,canAdj =CanonicalAdj@adj,n=Length[adj], vertices},
	vertices=Range[n];
	
	allLCs={canAdj};
	newGraphs= allLCs;
	
	While[Length[newGraphs]>0,
	
	newGraphs =LocalComplementation@@@Tuples[{newGraphs,vertices}];
	newGraphs=DeleteDuplicates[Map[CanonicalAdj, newGraphs]];
	
	newGraphs=Complement[newGraphs,allLCs];
	
	allLCs=Union[newGraphs,allLCs];
	
	];
	allLCs
 ];

  areLocalComplements[adj1_, adj2_] := Module[{LCorbit1, evals},

	LCorbit1=LCorbit[adj1];
	
	MemberQ[LCorbit1, CanonicalAdj@adj2]
 ];

  (* ========================================= *)
  (* Foliage Partition                         *)
  (* ========================================= *)

  areVerticesRelated[adj_, n_, vertices_, pair_] := Module[{v1, v2, complement, pairs, evals},
    v1 = pair[[1]];
    v2 = pair[[2]];
    complement = Complement[vertices, {v1, v2}];
    pairs = Subsets[complement, {2}];
    evals = (adj[[v1, #[[1]]]]*adj[[v2, #[[2]]]] == adj[[v1, #[[2]]]]*adj[[v2, #[[1]]]]) & /@ pairs;
    Not[MemberQ[evals, False]]
  ];

  FoliagePartition[adj_] := Module[{n, vertices, connectedComponents, foliagePartition,
  possibleRelations, i, j},
    n = Length[adj[[1]]];
    vertices = Range[n];
    connectedComponents = ConnectedComponents[AdjacencyGraph[adj]];
    foliagePartition = ({#} &) /@ vertices;
    Do[
      possibleRelations = Subsets[component, {2}];
      Do[
        If[areVerticesRelated[adj, n, vertices, pair],
          i = First@FirstPosition[foliagePartition, pair[[1]]];
          j = First@FirstPosition[foliagePartition, pair[[2]]];
          foliagePartition = ReplacePart[
            foliagePartition,
            {i -> Union[foliagePartition[[i]], foliagePartition[[j]]], j -> Nothing}
          ];
        ],
        {pair, possibleRelations}
      ],
      {component, connectedComponents}
    ];
    foliagePartition
  ];

  (* ========================================= *)
  (* Entropy Vector Utilities                  *)
  (* ========================================= *)

  Bipartitions[n_]:= Module[{qubits, subsystems, bipartitions},
	qubits=ToString/@Range[n];
	subsystems=Subsets[qubits,{1,Floor[n/2]}];
	bipartitions=DeleteDuplicates[Map[Sort@{Sort@ToExpression@#,Sort@ToExpression@Complement[qubits,#]}&,subsystems]] ;
	bipartitions =Map[{Map[ToString,#[[1]]], Map[ToString,#[[2]]]}&, bipartitions]
	];

  ProjectionOperator[tableaux_, complements_List] := Module[{n, projections, colIndices, ctableau},
    n = ToExpression[Length[tableaux[[1]]]];
    projections = Table[
      colIndices = Flatten[{complement, complement + n}];
      Table[
        ctableau = tableau;
        ctableau[[All, colIndices]] = 0;
        ctableau,
        {tableau, tableaux}
      ],
      {complement, complements}
    ];
    Transpose[projections]
  ];

  ComputeEntropyVectors[tableaux_] := Module[{n, bipartitions, A, nA, B, PAS},
    n = ToExpression[Length[tableaux[[1]]]];
    bipartitions = Bipartitions[n];
    A = bipartitions[[All, 1]];
    nA = Length /@ A;
    B = ToExpression /@ bipartitions[[All, 2]];
    PAS = ProjectionOperator[tableaux, B];
    Table[
      MatrixRank[#, Modulus -> 2] & /@ PAS[[i]] - nA,
      {i, 1, Length[PAS]}
    ]
  ];

  ComputeReducedRankVector[tableaux_] := Module[{n, bipartitions, B, PAS},
    n = ToExpression[Length[tableaux[[1]]]];
    bipartitions = Bipartitions[n];
    B = ToExpression /@ bipartitions[[All, 2]];
    PAS = ProjectionOperator[tableaux, B];
    Table[
      MatrixRank[#, Modulus -> 2] & /@ PAS[[i]],
      {i, 1, Length[PAS]}
    ]
  ];

  ComputeFullRankVector[tableaux_] := Module[{n, bipartitions, A, B, PAS, nA, nB, dAB, reduced},
    n = ToExpression[Length[tableaux[[1]]]];
    bipartitions = Bipartitions[n];
    A = bipartitions[[All, 1]];
    B = ToExpression /@ bipartitions[[All, 2]];
    PAS = ProjectionOperator[tableaux, B];
    nA = Length /@ A;
    nB = Length /@ B;
    dAB = Reverse[nB - nA];
    Table[
      reduced = MatrixRank[#, Modulus -> 2] & /@ PAS[[i]];
      Join[reduced, Reverse[reduced] + dAB, {n}],
      {i, 1, Length[PAS]}
    ]
  ];

  connectionsAdj[adjs_, bipartitions_, n_Integer?Positive] := Module[{connAdjs, columns, rows, cadj},
    connAdjs = Table[
      columns = bipartition[[1]];
      rows = bipartition[[2]];
      Table[
        cadj = adj;
        cadj[[All, columns]] = 0;
        cadj[[rows, All]] = 0;
        cadj,
        {adj, adjs}
      ],
      {bipartition, bipartitions}
    ];
    Transpose[connAdjs]
  ];

  ComputeEntropyVectorsAdj[adjs_] := Module[{n, bipartitions, cADJS, Svectors},
    n = ToExpression[Length[adjs[[1]]]];
    bipartitions = ToExpression[Bipartitions[n]];
    cADJS = connectionsAdj[adjs, bipartitions, n];
    Svectors = Map[
      Map[MatrixRank[#, Modulus -> 2] &, #] &,
      cADJS
    ];
    Svectors
  ];

  ComputeEntropyVectorsFromFile[filePath_String, n_Integer?Positive] := Module[
    {bipartitions, chunkSize = 10000, allSvectors = {}, stream, chunk, cADJS, Stemp},
    bipartitions = ToExpression[Bipartitions[n]];
    stream = OpenRead[filePath];
    If[stream === $Failed,
      Print["Error: Could not open the file."];
      Return[$Failed];
    ];
    While[True,
      chunk = ReadList[stream, "Expression", chunkSize];
      If[chunk === {}, Break[]];
      cADJS = connectionsAdj[chunk, bipartitions, n];
      Stemp = Map[MatrixRank[#, Modulus -> 2] &, cADJS, {2}];
      allSvectors = Join[allSvectors, Stemp];
    ];
    Close[stream];
    Return[Union[allSvectors]];
  ];

  nontrivialInterAdjs[foliagePartition_, Ematrices_, bipartitions_] :=
    Module[{Ehats, Ahat, APhat},
      Ehats = Table[
        Ahat =
          Select[foliagePartition, (Intersection[bipartition[[1]], #] != {}) &];
        Ahat = Flatten[Map[Position[foliagePartition, #] &, Ahat]];
        APhat =
          Select[foliagePartition, (Intersection[bipartition[[2]], #] != {}) &];
        APhat = Flatten[Map[Position[foliagePartition, #] &, APhat]];
        Table[matrix[[Ahat, APhat]], {matrix, Ematrices}],
        {bipartition, bipartitions}
      ];
      Transpose[Ehats]
    ];

  EntropyVectorsFoliageRep[foliagePartition_, Ematrices_] := Module[
    {n, bipartitions, nonTInterAdjs, Svectors},
    n = ToExpression[Length[Flatten[foliagePartition]]];
    bipartitions = ToExpression[Bipartitions[n]];
    nonTInterAdjs = nontrivialInterAdjs[foliagePartition, Ematrices, bipartitions];
    Svectors = Map[
      Map[MatrixRank[#, Modulus -> 2] &, #] &,
      nonTInterAdjs
    ];
    Svectors
  ];

  (* ========================================= *)
  (* Non-trivial Column Intersections          *)
  (* ========================================= *)

  isColIntsNonTrivial[A_, B_, C_] := Module[{m, nA, nB, nC, blockMatrix, ns, xvecs, evals},
    {m, nA} = Dimensions[A];
    nB = Dimensions[B][[2]];
    nC = Dimensions[C][[2]];
    blockMatrix = ArrayFlatten[{{A, -B, ConstantArray[0, {m, nC}]}, {A, ConstantArray[0, {m, nB}],
  -C}}];
    ns = NullSpace[blockMatrix, Modulus -> 2];
    xvecs = ns[[All, ;; nA]];
    evals = Mod[A . Transpose[#], 2] & /@ xvecs;
    MemberQ[Flatten[evals], 1]
  ];

  nonTrivialColIntersections[adj_, triples_] := Module[{evals, n, vertices, positions},
    n = Length[adj[[1]]];
    vertices = Range[n];
    evals = isColIntsNonTrivial[adj[[#[[1]], #[[2]]]], adj[[#[[1]], #[[3]]]], adj[[#[[1]],
  Complement[vertices, Flatten[#]]]]] & /@ triples;
    positions = Flatten[Position[evals, True]];
    {MemberQ[evals, True], triples[[positions]]}
  ];

  (* ========================================= *)
  (* Section 4.3 Graph Families                *)
  (* ========================================= *)

  allBinaryMatrices[m_Integer?NonNegative, n_Integer?NonNegative] := Tuples[{0, 1}, {m, n}];

  allGeneralizedGraphsnoInner[n_Integer?Positive] := Module[{baseMatrix, sizesCenter, allGraphs},
    baseMatrix = ConstantArray[0, {n, n}];
    sizesCenter = Range[n - 3];
    allGraphs = ParallelTable[
      Module[{Cset, adjacenciesList},
        Cset = Range[sizeC];
        adjacenciesList = Tuples[{0, 1}, {sizeC, n - sizeC}];
        Last@Reap[
          Do[
            With[{temp = baseMatrix + Transpose[ReplacePart[baseMatrix, Thread[Cset -> adj]]]},
              Sow[temp]
            ],
            {adj, adjacenciesList}
          ]
        ]
      ],
      {sizeC, sizesCenter}
    ];
    Flatten[allGraphs, 2]
  ];

  allGeneralizedGraphsnoInnerNonTrivialInter[n_Integer?Positive] := Module[{baseMatrix, sizesCenter,
  disjointTriples, allGraphs},
    baseMatrix = ConstantArray[0, {n, n}];
    sizesCenter = Range[n - 3];
    SetSharedFunction[DisjointTriples, nonTrivialColIntersections];
    disjointTriples = DisjointTriples[n];
    allGraphs = Table[
      Module[{Cset, adjacenciesList, relevantTriples},
        Cset = Range[sizeC];
        adjacenciesList = Tuples[{0, 1}, {sizeC, n - sizeC}];
        relevantTriples = Select[disjointTriples, MemberQ[Cset]];
        relevantTriples = (Prepend[DeleteCases[#, Cset], Cset] &) /@ relevantTriples;
        Last@Reap[
          Do[
            With[{temp = baseMatrix + Transpose[ReplacePart[baseMatrix, Thread[Cset -> adj]]]},
              If[First@nonTrivialColIntersections[temp, relevantTriples], Sow[temp]]
            ],
            {adj, adjacenciesList}
          ]
        ]
      ],
      {sizeC, sizesCenter}
    ];
    Flatten[allGraphs, 2]
  ];

  allGeneralizedGraphsnoInnerNonTrivialInterInfo[n_Integer?Positive] := Module[{baseMatrix,
  sizesCenter, disjointTriples, allGraphs},
    baseMatrix = ConstantArray[0, {n, n}];
    sizesCenter = Range[n - 3];
    SetSharedFunction[DisjointTriples, nonTrivialColIntersections];
    disjointTriples = DisjointTriples[n];
    allGraphs = Table[
      Module[{Cset, adjacenciesList, relevantTriples},
        Cset = Range[sizeC];
        adjacenciesList = Tuples[{0, 1}, {sizeC, n - sizeC}];
        relevantTriples = Select[disjointTriples, MemberQ[Cset]];
        relevantTriples = (Prepend[DeleteCases[#, Cset], Cset] &) /@ relevantTriples;
        Last@Reap[
          Do[
            With[{temp = baseMatrix + Transpose[ReplacePart[baseMatrix, Thread[Cset -> adj]]], eval},
              eval = nonTrivialColIntersections[temp, relevantTriples];
              If[First@eval, Sow[{temp, eval[[2]]}]]
            ],
            {adj, adjacenciesList}
          ]
        ]
      ],
      {sizeC, sizesCenter}
    ];
    Flatten[allGraphs, 2]
  ];

  allGeneralizedGraphsnoInnerNonTrivialInterInfoParallel[n_Integer?Positive] := Module[{baseMatrix,
  sizesCenter, disjointTriples, allGraphs},
    baseMatrix = ConstantArray[0, {n, n}];
    sizesCenter = Range[n - 3];
    SetSharedFunction[DisjointTriples, nonTrivialColIntersections];
    disjointTriples = DisjointTriples[n];
    allGraphs = ParallelTable[
      Module[{Cset, adjacenciesList, relevantTriples},
        Cset = Range[sizeC];
        adjacenciesList = Tuples[{0, 1}, {sizeC, n - sizeC}];
        relevantTriples = Select[disjointTriples, MemberQ[Cset]];
        relevantTriples = (Prepend[DeleteCases[#, Cset], Cset] &) /@ relevantTriples;
        Last@Reap[
          Do[
            With[{temp = baseMatrix + Transpose[ReplacePart[baseMatrix, Thread[Cset -> adj]]], eval},
              eval = nonTrivialColIntersections[temp, relevantTriples];
              If[First@eval, Sow[{temp, eval[[2]]}]]
            ],
            {adj, adjacenciesList}
          ]
        ]
      ],
      {sizeC, sizesCenter}
    ];
    Flatten[allGraphs, 2]
  ];

  SampleAllGeneralizedGraphsnoInnerNonTrivialInterInfoTriples[n_Integer?Positive, numSamples_Integer?
  Positive] := Module[
    {baseMatrix, sizesCenter, disjointTriples, allGraphs},
    baseMatrix = ConstantArray[0, {n, n}];
    sizesCenter = Range[n - 3];
    SetSharedFunction[DisjointTriples, nonTrivialColIntersections];
    disjointTriples = DisjointTriples[n];
    allGraphs = Table[
      Module[{Cset, adjacenciesList, relevantTriples},
        Cset = Range[sizeC];
        adjacenciesList = Table[RandomInteger[{0, 1}, {sizeC, n - sizeC}], {numSamples}];
        relevantTriples = Select[disjointTriples, MemberQ[Cset]];
        relevantTriples = (Prepend[DeleteCases[#, Cset], Cset] &) /@ relevantTriples;
        Last@Reap[
          Do[
            With[{temp = baseMatrix + Transpose[ReplacePart[baseMatrix, Thread[Cset -> adj]]], eval},
              eval = nonTrivialColIntersections[temp, relevantTriples];
              If[First@eval, Sow[{temp, eval[[2]]}]]
            ],
            {adj, adjacenciesList}
          ]
        ]
      ],
      {sizeC, sizesCenter}
    ];
    Flatten[allGraphs, 2]
  ];

  (* ========================================= *)
  (* Entropy Vector Symmetrization             *)
  (* ========================================= *)

  EquivalentEntropyVectors[n_Integer?Positive] := Module[{qubits, baseVector, Qexchange, permuted,
  eqVectors},
    qubits = ToString /@ Range[n];
    baseVector = (Bipartitions[n])[[All, 1]];
    Qexchange = Rest[Permutations[qubits]];
    eqVectors = Table[
      permuted = baseVector /. Thread[qubits -> Qexchange[[i]]];
      permuted = SortBy[#, ToExpression] & /@ permuted;
      permuted = Map[
        If[MemberQ[baseVector, #],
          #,
          SortBy[Complement[qubits, #], ToExpression]
        ] &,
        permuted
      ];
      permuted,
      {i, Length[Qexchange]}
    ];
    eqVectors = Union[eqVectors];
    eqVectors = Complement[eqVectors, {baseVector}];
    {baseVector, eqVectors}
  ];

  permRules[Sv_, baseVector_] := permRules[Sv, baseVector] = Dispatch@Thread[baseVector -> Sv];

  PermuteEntropyVector[Sv_, baseVector_, eqVectors_] := DeleteDuplicates[eqVectors /. permRules[Sv,
  baseVector]];

  UniqueOrbits[SVectors_, baseVector_, eqVectors_] := Module[{i = 1, Sv, equivOrbits, uniqueOrbits},
    uniqueOrbits = SVectors;
    While[i <= Length[uniqueOrbits],
      Sv = uniqueOrbits[[i]];
      equivOrbits = PermuteEntropyVector[Sv, baseVector, eqVectors];
      equivOrbits = Complement[equivOrbits, {Sv}];
      uniqueOrbits = Complement[uniqueOrbits, equivOrbits];
      i++;
    ];
    uniqueOrbits
  ];

  ParallelUniqueOrbits[SVectors_, baseVector_, eqVectors_] := Module[{i = 1, Sv, equivOrbits,
  uniqueOrbits},
    uniqueOrbits = AssociationThread[SVectors -> ConstantArray[True, Length[SVectors]]];
    While[i <= Length[Keys[uniqueOrbits]],
      Sv = Keys[uniqueOrbits][[i]];
      equivOrbits = PermuteEntropyVector[Sv, baseVector, eqVectors];
      equivOrbits = Complement[equivOrbits, {Sv}];
      KeyDropFrom[uniqueOrbits, equivOrbits];
      i++;
    ];
    Keys[uniqueOrbits]
  ];

  findFirstIndex[largeList_, smallList_] := FirstPosition[largeList, _?(MemberQ[smallList, #] &)];
  findallIndices[largeList_, smallList_] := Position[largeList, _?(MemberQ[smallList, #] &)];

  SampleStates[tableaux_, sVectors_, symmetryOrbits_, uniqueOrbits_] := Module[{positions, sTableaus,
  sStab},
    positions = Flatten[findFirstIndex[sVectors, #] & /@ symmetryOrbits, {2, 1}];
    sTableaus = tableaux[[#]] & /@ positions;
    sStab = tableauToStabilizerGenerators /@ sTableaus;
    AssociationThread[uniqueOrbits -> sStab]
  ];

  SymmetryGroup[tableaux_, sVectors_, symmetryOrbits_, uniqueOrbits_] := Module[{positions,
  sTableaus},
    positions = Flatten[findallIndices[sVectors, #] & /@ symmetryOrbits, {3, 1}];
    sTableaus = tableaux[[#]] & /@ positions;
    AssociationThread[uniqueOrbits -> sTableaus]
  ];

  (* ========================================= *)
  (* Disjoint Triples & Equivalence            *)
  (* ========================================= *)

  DisjointTriples[n_]:=Module[{verts=Range[n],subs},subs=Subsets[verts,{1,n-1}];
Parallelize@Select[Subsets[subs,{3}],And@@(DisjointQ@@@Subsets[#,{2}])&]];

  EquivalentTriples[triple_, triple2_, eqQubits_] := Module[{Qexchange, eqTriples},
    Qexchange = Rest[Permutations[eqQubits]];
    eqTriples = Table[
      SortBy[#, ToExpression] & /@ (triple /. Thread[eqQubits -> Qexchange[[i]]]),
      {i, Length[Qexchange]}
    ];
    MemberQ[DeleteDuplicates[eqTriples], triple2]
  ];

  (* ========================================= *)
  (* MMI Evaluation                            *)
  (* ========================================= *)

  MMIChecker[ReducedEntropyVector_, DisjointTriples_]:=Block[{FullEntropyVector,EVecTranspose,EVecIndices,PossibleTriples,evalutationTable},

	FullEntropyVector=Join[ReducedEntropyVector,Reverse[ReducedEntropyVector],{0}];
	EVecTranspose=Flatten[Transpose[{FullEntropyVector}]];
	EVecIndices=DeleteCases[Subsets[Range[Log[2,Length[FullEntropyVector]+1]]],{}];
	
	AllTrue[
	Table[EVecTranspose[[Position[EVecIndices,Sort[Union[DisjointTriples[[k,1]],DisjointTriples[[k,2]]]]][[1,1]]]]+EVecTranspose[[Position[EVecIndices,Sort[Union[DisjointTriples[[k,1]],DisjointTriples[[k,3]]]]][[1,1]]]]+EVecTranspose[[Position[EVecIndices,Sort[Union[DisjointTriples[[k,2]],DisjointTriples[[k,3]]]]][[1,1]]]]>=EVecTranspose[[Position[EVecIndices,Sort[DisjointTriples[[k,1]]]][[1,1]]]]+ EVecTranspose[[Position[EVecIndices,Sort[DisjointTriples[[k,2]]]][[1,1]]]] + EVecTranspose[[Position[EVecIndices,Sort[DisjointTriples[[k,3]]]][[1,1]]]]+ EVecTranspose[[Position[EVecIndices,Sort[Union[DisjointTriples[[k,1]],DisjointTriples[[k,2]],DisjointTriples[[k,3]]]]][[1,1]]]],
	{k,Length[DisjointTriples]}],
	#&]
];

  MMICheckFV[FullEntropyVector_] := Module[
    {EVecTranspose, EVecIndices, PossibleTriples, DisjointTriples, evalutationTable},
    EVecTranspose = Flatten[Transpose[{FullEntropyVector}]];
    EVecIndices = DeleteCases[Subsets[Range[Log[2, Length[FullEntropyVector] + 1]]], {}];
    PossibleTriples = Subsets[DeleteCases[Subsets[Range[Log[2, Length[FullEntropyVector] + 1]], 3],
  {}], {3}];
    DisjointTriples = DeleteCases[
      Table[
        If[
          And[
            DisjointQ[PossibleTriples[[j, 1]], PossibleTriples[[j, 2]]],
            DisjointQ[PossibleTriples[[j, 2]], PossibleTriples[[j, 3]]],
            DisjointQ[PossibleTriples[[j, 1]], PossibleTriples[[j, 3]]]
          ],
          PossibleTriples[[j]],
          Null
        ],
        {j, Length[PossibleTriples]}
      ],
      Null
    ];
    evalutationTable = Table[
      Which[
        (* Saturates *)
        EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Union[DisjointTriples[[k, 1]],
  DisjointTriples[[k, 2]]]]]]] +
          EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Union[DisjointTriples[[k, 1]],
  DisjointTriples[[k, 3]]]]]]] +
          EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Union[DisjointTriples[[k, 2]],
  DisjointTriples[[k, 3]]]]]]] ==
          EVecTranspose[[First@FirstPosition[EVecIndices, Sort[DisjointTriples[[k, 1]]]]]] +
          EVecTranspose[[First@FirstPosition[EVecIndices, Sort[DisjointTriples[[k, 2]]]]]] +
          EVecTranspose[[First@FirstPosition[EVecIndices, Sort[DisjointTriples[[k, 3]]]]]] +
          EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Union @@ DisjointTriples[[k]]]]]],
        "Saturates",
        (* Satisfies *)
        EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Union[DisjointTriples[[k, 1]],
  DisjointTriples[[k, 2]]]]]]] +
          EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Union[DisjointTriples[[k, 1]],
  DisjointTriples[[k, 3]]]]]]] +
          EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Union[DisjointTriples[[k, 2]],
  DisjointTriples[[k, 3]]]]]]] >
          EVecTranspose[[First@FirstPosition[EVecIndices, Sort[DisjointTriples[[k, 1]]]]]] +
          EVecTranspose[[First@FirstPosition[EVecIndices, Sort[DisjointTriples[[k, 2]]]]]] +
          EVecTranspose[[First@FirstPosition[EVecIndices, Sort[DisjointTriples[[k, 3]]]]]] +
          EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Union @@ DisjointTriples[[k]]]]]],
        "Satisfies",
        True,
        "Fails"
      ],
      {k, Length[DisjointTriples]}
    ];
    Which[
      MemberQ[evalutationTable, "Fails"], "Fails",
      AllTrue[evalutationTable, # == "Saturates" &], "Saturates",
      True, "Satisfies"
    ]
  ];

  MMIcheckerMultiple[sVectors_] := Module[{n, disjointTriples},
    n = Log[2, 2 Length[sVectors[[1]]] + 2];
    disjointTriples = DisjointTriples[n];
    MMIChecker[#, disjointTriples] & /@ sVectors
  ];

  MMICheckerInputSelectTriples[ReducedEntropyVector_, Triples_] := Module[
    {FullEntropyVector, EVecTranspose, EVecIndices, evalutationTable},
    FullEntropyVector = Join[ReducedEntropyVector, Reverse[ReducedEntropyVector], {0}];
    EVecTranspose = Flatten[Transpose[{FullEntropyVector}]];
    EVecIndices = DeleteCases[Subsets[Range[Log[2, Length[FullEntropyVector] + 1]]], {}];
    evalutationTable = Table[
      Which[
        EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Union[Triples[[k, 1]], Triples[[k,
  2]]]]]]] +
          EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Union[Triples[[k, 1]], Triples[[k,
  3]]]]]]] +
          EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Union[Triples[[k, 2]], Triples[[k,
  3]]]]]]] ==
          EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Triples[[k, 1]]]]]] +
          EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Triples[[k, 2]]]]]] +
          EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Triples[[k, 3]]]]]] +
          EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Union @@ Triples[[k]]]]]],
        "Saturates",
        EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Union[Triples[[k, 1]], Triples[[k,
  2]]]]]]] +
          EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Union[Triples[[k, 1]], Triples[[k,
  3]]]]]]] +
          EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Union[Triples[[k, 2]], Triples[[k,
  3]]]]]]] >
          EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Triples[[k, 1]]]]]] +
          EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Triples[[k, 2]]]]]] +
          EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Triples[[k, 3]]]]]] +
          EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Union @@ Triples[[k]]]]]],
        "Satisfies",
        True,
        "Fails"
      ],
      {k, Length[Triples]}
    ];
    evalutationTable
  ];

  MMICheckerEntropyInput[ReducedEntropyVector_] := Module[
    {FullEntropyVector, EVecTranspose, EVecIndices, PossibleTriples, DisjointTriples},
    FullEntropyVector = Join[ReducedEntropyVector, Reverse[ReducedEntropyVector], {0}];
    EVecTranspose = Flatten[Transpose[{FullEntropyVector}]];
    EVecIndices = DeleteCases[Subsets[Range[Log[2, Length[FullEntropyVector] + 1]]], {}];
    PossibleTriples = Subsets[DeleteCases[Subsets[Range[Log[2, Length[FullEntropyVector] + 1]], 3],
  {}], {3}];
    DisjointTriples = DeleteCases[
      Table[
        If[
          And[
            DisjointQ[PossibleTriples[[j, 1]], PossibleTriples[[j, 2]]],
            DisjointQ[PossibleTriples[[j, 2]], PossibleTriples[[j, 3]]],
            DisjointQ[PossibleTriples[[j, 1]], PossibleTriples[[j, 3]]]
          ],
          PossibleTriples[[j]],
          Null
        ],
        {j, Length[PossibleTriples]}
      ],
      Null
    ];
    MMICheckerInputSelectTriples[ReducedEntropyVector, DisjointTriples]
  ];

  MMICheckerEntropyFullEV[FullEntropyVector_] := Module[
    {EVecTranspose, EVecIndices, PossibleTriples, DisjointTriples},
    EVecTranspose = Flatten[Transpose[{FullEntropyVector}]];
    EVecIndices = DeleteCases[Subsets[Range[Log[2, Length[FullEntropyVector] + 1]]], {}];
    PossibleTriples = Subsets[DeleteCases[Subsets[Range[Log[2, Length[FullEntropyVector] + 1]], 3],
  {}], {3}];
    DisjointTriples = DeleteCases[
      Table[
        If[
          And[
            DisjointQ[PossibleTriples[[j, 1]], PossibleTriples[[j, 2]]],
            DisjointQ[PossibleTriples[[j, 2]], PossibleTriples[[j, 3]]],
            DisjointQ[PossibleTriples[[j, 1]], PossibleTriples[[j, 3]]]
          ],
          PossibleTriples[[j]],
          Null
        ],
        {j, Length[PossibleTriples]}
      ],
      Null
    ];
    MMICheckerInputSelectTriples[
      FullEntropyVector[[;; -2]] - Reverse[FullEntropyVector[[;; -2]]]/0 + 0,
      DisjointTriples
    ]
  ];

  MMICheckerVec[ReducedEntropyVector_, EVecIndices_, DisjointTriples_] := Module[
    {FullEntropyVector, EVecTranspose},
    FullEntropyVector = Join[ReducedEntropyVector, Reverse[ReducedEntropyVector], {0}];
    EVecTranspose = Flatten[Transpose[{FullEntropyVector}]];
    Table[
      {
        EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Union[DisjointTriples[[k, 1]],
  DisjointTriples[[k, 2]]]]]]],
        EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Union[DisjointTriples[[k, 1]],
  DisjointTriples[[k, 3]]]]]]],
        EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Union[DisjointTriples[[k, 2]],
  DisjointTriples[[k, 3]]]]]]],
        EVecTranspose[[First@FirstPosition[EVecIndices, Sort[DisjointTriples[[k, 1]]]]]],
        EVecTranspose[[First@FirstPosition[EVecIndices, Sort[DisjointTriples[[k, 2]]]]]],
        EVecTranspose[[First@FirstPosition[EVecIndices, Sort[DisjointTriples[[k, 3]]]]]],
        EVecTranspose[[First@FirstPosition[EVecIndices, Sort[Union @@ DisjointTriples[[k]]]]]]
      },
      {k, Length[DisjointTriples]}
    ]
  ];

  MMIViolatingStates[tableaux_] := Module[{n, sVectors, fullSVectors, PossibleTriples,
  DisjointTriples, evals, pos},
    n = Length[tableaux[[1]]];
    sVectors = ComputeEntropyVectors[tableaux];
    fullSVectors = (Join[#, Reverse[#], {0}] &) /@ sVectors;
    PossibleTriples = Subsets[DeleteCases[Subsets[Range[Log[2, 2^n]], 3], {}], {3}];
    DisjointTriples = DeleteCases[
      Table[
        If[
          And[
            DisjointQ[PossibleTriples[[j, 1]], PossibleTriples[[j, 2]]],
            DisjointQ[PossibleTriples[[j, 2]], PossibleTriples[[j, 3]]],
            DisjointQ[PossibleTriples[[j, 1]], PossibleTriples[[j, 3]]]
          ],
          PossibleTriples[[j]],
          Null
        ],
        {j, Length[PossibleTriples]}
      ],
      Null
    ];
    evals = MMIChecker[#, DisjointTriples] & /@ sVectors;
    pos = Flatten[Position[evals, False]];
    tableaux[[pos]]
  ];

  finddMinimalNontrivialTriple[adj_, disjointTriples_] := Module[
    {triplesCenter, evalCenter, evalIntersections, nonTrivialTriples, lengths, pos, center,
  newTripleCenter, centerPos, vertices},
    vertices = Range[Length[adj]];
    triplesCenter = Table[
      evalCenter = isValidCIJ[adj, triple, vertices];
      If[MemberQ[evalCenter, True],
        centerPos = FirstPosition[evalCenter, True];
        center = triple[[centerPos]];
        newTripleCenter = DeleteDuplicates[Join[center, triple]],
        Nothing
      ],
      {triple, disjointTriples}
    ];
    evalIntersections = nonTrivialColIntersections[adj, triplesCenter];
    If[First@evalIntersections === True,
      nonTrivialTriples = evalIntersections[[2]];
      lengths = Length@Flatten[#] & /@ nonTrivialTriples;
      pos = FirstPosition[lengths, Min[lengths]];
      First@nonTrivialTriples[[pos]],
      {}
    ]
  ];

  findLargestNontrivialTriple[adj_, disjointTriples_] := Module[
    {triplesCenter, evalCenter, evalIntersections, nonTrivialTriples, lengths, pos, center,
  newTripleCenter, centerPos, vertices},
    vertices = Range[Length[adj]];
    triplesCenter = Table[
      evalCenter = isValidCIJ[adj, triple, vertices];
      If[MemberQ[evalCenter, True],
        centerPos = FirstPosition[evalCenter, True];
        center = triple[[centerPos]];
        newTripleCenter = DeleteDuplicates[Join[center, triple]],
        Nothing
      ],
      {triple, disjointTriples}
    ];
    evalIntersections = nonTrivialColIntersections[adj, triplesCenter];
    If[First@evalIntersections === True,
      nonTrivialTriples = evalIntersections[[2]];
      lengths = Length@Flatten[#] & /@ nonTrivialTriples;
      pos = FirstPosition[lengths, Max[lengths]];
      First@nonTrivialTriples[[pos]],
      {}
    ]
  ];

  (* ========================================= *)
  (* Generalized Star Detection                *)
  (* ========================================= *)

  IsGeneralizedStarQ[adj_] := Module[{n, vertices, centers, centersAndLeaves,
  isGeneralizedStarEvals},
    n = Length[adj];
    vertices = Range[n];
    centers = Subsets[vertices, {1, n - 3}];
    centersAndLeaves = {#, Subsets[Complement[vertices, #], {3}]} & /@ centers;
    isGeneralizedStarEvals = Flatten@Table[
      With[{C = centerAndL[[1]]},
        (Min[Total[adj[[C, #]], {1}]] >= 1 && Total[Flatten@adj[[#, #]]] == 0) & /@
          centerAndL[[2]]
      ],
      {centerAndL, centersAndLeaves}
    ];
    MemberQ[isGeneralizedStarEvals, True]
  ];

  IsGeneralizedStarQdisC[adj_, leaves_] := Module[{isGeneralizedStarEvals},
    isGeneralizedStarEvals = (Total[Flatten@adj[[#, #]]] == 0) & /@ leaves;
    MemberQ[isGeneralizedStarEvals, True]
  ];

  isValidCIJ[adj_, triple_, vertices_] := Module[{R1, R2, R3, validCenterQ},
    {R1, R2, R3} = triple;
    validCenterQ[Center_, subI_, subJ_] := Module[{subK},
      subK = Complement[vertices, Join[Center, subI, subJ]];
      Max@Total[adj[[Center, subI]], {2}] > 0 &&
        Max@Total[adj[[Center, subJ]], {2}] > 0 &&
        Max@Total[adj[[Center, subK]], {2}] > 0 &&
        Total[Total@adj[[subI, subJ]]] == 0 &&
        Total[Flatten@adj[[subI, subK]]] == 0 &&
        Total[Flatten@adj[[subJ, subK]]] == 0
    ];
    {
      validCenterQ[R1, R2, R3],
      validCenterQ[R2, R1, R3],
      validCenterQ[R3, R2, R1]
    }
  ];

  applyPermutations[lists_, vertices_, perms_] := Module[{rules},
    rules = Thread[vertices -> #] & /@ perms;
    DeleteDuplicates@Flatten[ParallelMap[(canSubsystems@(# /. rules)) &, lists], 1]
  ];

  canSubsystems[list_] := DeleteDuplicates@Map[Sort@Map[Sort, #] &, list];

  fourPartitions[n_Integer?Positive] :=
    Select[ResourceFunction["SetPartitions"][Range[n]], Length[#] == 4 &];

  validCenters[adj_, fourParts_] := Module[{centers, IJK, eval, isCenter},
    centers = Table[
      IJK = Complement[partition, {element}];
      eval = {adj[[IJK[[1]], IJK[[2]]]], adj[[IJK[[1]], IJK[[3]]]], adj[[IJK[[2]], IJK[[3]]]]};
      eval = Map[(Total@Flatten[#] === 0) &, eval];
      isCenter = AllTrue[eval, TrueQ];
      If[isCenter, {element, IJK}, Nothing],
      {partition, fourParts}, {element, partition}
    ];
    centers = Flatten[centers, {2}];
    If[centers === {}, {"NF"}, centers]
  ];

  findStarGraphReps[localComplementationOrbits_, fourParts_] := Module[{generalizedStarRepsCenters},
    generalizedStarRepsCenters = ParallelTable[
      Map[{#, First@validCenters[#, fourParts]} &, orbit],
      {orbit, localComplementationOrbits}
    ];
    generalizedStarRepsCenters
  ];

  InstancesCIJ[adj_, triples_, n_Integer?Positive] := Module[{evalCenter, centerPos, center,
  newTripleCenter, vertices},
    vertices = Range[n];
    Table[
      evalCenter = isValidCIJ[adj, triple, vertices];
      If[MemberQ[evalCenter, True],
        centerPos = FirstPosition[evalCenter, True];
        center = triple[[centerPos]];
        newTripleCenter = DeleteDuplicates[Join[center, triple]],
        Nothing
      ],
      {triple, triples}
    ]
  ];

  searchNonTrivialIntLC[adjs_, triples_] := Module[{orbitsLC, nonTrivialEvals, candidates,
  candidatesEdges, minimalCase, edgesOrbit},
    orbitsLC = ParallelMap[LCorbit, adjs];
    ParallelTable[
      nonTrivialEvals = (findLargestNontrivialTriple[#, triples] === {}) & /@ orbit;
      candidates = Pick[orbit, nonTrivialEvals, False];
      If[candidates =!= {},
        candidatesEdges = EdgeCount[AdjacencyGraph[#]] & /@ candidates;
        minimalCase = First@Pick[candidates, candidatesEdges, Min[candidatesEdges]];
        {minimalCase, findLargestNontrivialTriple[minimalCase, triples]},
        edgesOrbit = EdgeCount[AdjacencyGraph[#]] & /@ orbit;
        minimalCase = First@Pick[orbit, edgesOrbit, Min[edgesOrbit]];
        {minimalCase, {}}
      ],
      {orbit, orbitsLC}
    ]
  ];

  findMinimalEdgeRep[LCorbitsAndCenters_] := Module[{adjs, evals},
    ParallelTable[
      adjs = orbit[[All, 1]];
      evals = Map[EdgeCount@AdjacencyGraph[#] &, adjs];
      First@Pick[orbit, evals, Min[evals]],
      {orbit, LCorbitsAndCenters}
    ]
  ];

  findBiggestCenter[repCenters_] := Module[{centers, evals, CIJK},
    centers = repCenters[[2]][[All, 1]];
    evals = Map[Length, centers];
    CIJK = First@Pick[repCenters[[2]], evals, Max@evals];
    {repCenters[[1]], {CIJK[[1]], CIJK[[2, 1]], CIJK[[2, 2]]}}
  ];

  componentDiameters[adj_] := Module[{G},
    G = AdjacencyGraph[Unitize@adj, DirectedEdges -> False];
    GraphDiameter /@ (Subgraph[G, #] & /@ ConnectedComponents[G])
  ];

  pathTieBreaker[adj_] := Module[{diams = componentDiameters[adj], maxD, cntMax, secondD, cntSecond},
    If[diams === {}, Return[{0, 0, 0, 0}]];
    maxD = Max[diams];
    cntMax = Count[diams, maxD];
    secondD = If[Length[DeleteDuplicates[diams]] >= 2, Max@DeleteCases[diams, maxD], 0];
    cntSecond = Count[diams, secondD];
    {maxD, cntMax, secondD, cntSecond}
  ];

  (* ========================================= *)
  (* TikZ Export Utilities                     *)
  (* ========================================= *)

  fmtNum[x_?NumericQ] :=
    ToString@NumberForm[N@x, {Infinity, 3}, NumberPoint -> ".", NumberSeparator -> ""];

  fmtInt[x_?NumericQ] := IntegerString@Round[x];

  fmtVec[list_List] /; Length[list] == 3 := "(" <> StringRiffle[fmtInt /@ list, ","] <> ")";
  fmtVec[x_?NumericQ] := "(" <> fmtInt[x] <> ")";
  fmtVec[___] := "";

  isPathG[g_Graph] :=
    Module[{deg = VertexDegree[g]},
      ConnectedGraphQ[g] && VertexCount[g] >= 1 && Max[deg] <= 2 &&
        (VertexCount[g] == 1 || Count[deg, 1] == 2)
    ];

  isCycleG[g_Graph] :=
    Module[{deg = VertexDegree[g]},
      ConnectedGraphQ[g] && VertexCount[g] >= 3 && AllTrue[deg, # == 2 &]
    ];

  pad01[p_, eps_ : 0.06] := eps + (1 - 2 eps) p;

  firstOr[l_, d_] := Replace[l, {{} -> d, {h_, ___} :> h}];

  orderPath[g_Graph] :=
    Module[{V = VertexList[g], ends, adj, start, cur, prev = None, out = {}, visited = <||>, step = 0,
      max = 2 VertexCount[g] + 5},
      Which[
        Length@V == 0, Return[{}],
        Length@V == 1, Return[V],
        ! ConnectedGraphQ[g], Return[V]
      ];
      If[!(Max@VertexDegree[g] <= 2 && Count[VertexDegree[g], 1] == 2), Return[V]];
      ends = Select[V, VertexDegree[g, #] == 1 &];
      If[Length@ends < 2, Return[V]];
      adj = AssociationMap[Normal@AdjacencyList[g, #] &, V];
      start = First@ends;
      cur = start;
      While[True,
        step++;
        If[step > max, Break[]];
        If[TrueQ@KeyExistsQ[visited, cur], Break[]];
        visited[cur] = True;
        AppendTo[out, cur];
        With[{nbrs = DeleteCases[Lookup[adj, cur, {}], prev]},
          If[nbrs === {} || Length[out] == Length@V, Break[]];
          prev = cur;
          cur = First@nbrs;
        ];
      ];
      If[out === {}, V, out]
    ];

  graphEmbeddingBox[g_Graph] :=
    graphEmbeddingBox[g] =
      Module[{V = VertexList[g], assoc, emb, pts, xs, ys, xmin, xmax, ymin, ymax, s, xc, yc, norm,
        ord, n},
        If[Length@V == 0, Return[<||>]];
        Which[
          isPathG[g],
          ord = orderPath[g];
          n = Length@ord;
          assoc =
            AssociationThread[
              ord -> Transpose@{
                If[n <= 1, {0.5}, Range[0, n - 1]/Max[n - 1, 1]],
                ConstantArray[0.5, n]
              }
            ];
          assoc =
            Join[
              assoc,
              AssociationThread[
                Complement[V, Keys@assoc],
                ConstantArray[{0.5, 0.5}, Length@Complement[V, Keys@assoc]]
              ]
            ];
          assoc = Map[pad01[#, 0.06] &, assoc],
          isCycleG[g],
          assoc = AssociationThread[V, 0.5 + 0.5 CirclePoints[Length@V]];
          assoc = Map[pad01[#, 0.06] &, assoc],
          True,
          emb = Which[
            AssociationQ@GraphEmbedding[g], GraphEmbedding[g],
            MatrixQ[GraphEmbedding[g], NumericQ] && Length@GraphEmbedding[g] == Length@V,
            AssociationThread[V, GraphEmbedding[g]],
            True,
            AssociationThread[V, CirclePoints[Length@V]]
          ];
          pts = Values@emb;
          If[pts === {}, Return[<||>]];
          xs = pts[[All, 1]];
          ys = pts[[All, 2]];
          xmin = Min@xs; xmax = Max@xs;
          ymin = Min@ys; ymax = Max@ys;
          s = Max[xmax - xmin, ymax - ymin];
          If[s == 0,
            assoc = AssociationThread[V, CirclePoints[Length@V]],
            xc = (xmin + xmax)/2.;
            yc = (ymin + ymax)/2.;
            norm[{x_, y_}] := pad01[{(x - xc)/s + 0.5, (y - yc)/s + 0.5}, 0.06];
            assoc = Map[Chop[Clip[norm[#], {0., 1.}], 10^-12] &, emb]
          ]
        ];
        assoc
      ];

  toTikZ[nodes_, edges_, bbox_] :=
    StringTemplate["\\begin{tikzpicture}[scale=1]
  `bbox`
  `nodes`
  \\begin{scope}[on background layer]
  `edges`
  \\end{scope}
  \\end{tikzpicture}"][<|"nodes" -> nodes, "edges" -> edges, "bbox" -> bbox|>];

  Options[ToTikZAdjTripleColor] = {CanvasCM -> {3.5, 3.5},
    EdgeStyle -> "line width=0.9pt,line cap=round", ShowLabels -> True,
    LabelSize -> "\\scriptsize", PreserveAspect -> True,
    BBoxPadCM -> 0.10, EvalData -> None, ShowEvalLabels -> True,
    EvalLabelSize -> "\\scriptsize", EvalLabelOffsetCM -> {0, 0.18},
    EvalHeaderOffsetCM -> 0.18, ColorScheme -> "CIJ", BorderTriple -> None};

  ToTikZAdjTripleColor[adj_?MatrixQ, triple_, OptionsPattern[]] :=
    Module[{g, V, pairs, coords, wx, wy, pad, idxMap, nameOf, C = {}, I = {}, J = {}, colorOff, nodes,
      edges, bbox, tag, evalData, showEval, evAssoc, off, headerOff, graphLevelQ, header, extraTop,
      scheme, styles, borderTriple, dotSet = {}, dashSet = {}, starSet = {}},
      g = AdjacencyGraph[Normal@adj, VertexLabels -> None, GraphLayout -> "SpringElectricalEmbedding"];
      V = VertexList[g];
      pairs = List @@@ EdgeList[g];
      coords = graphEmbeddingBox[g];
      If[coords === <||>, Return["\\begin{tikzpicture}\\end{tikzpicture}"]];
      {wx, wy} =
        If[ListQ@OptionValue[CanvasCM],
          OptionValue[CanvasCM],
          {OptionValue[CanvasCM], OptionValue[CanvasCM]}
        ];
      coords = Map[(#*{wx, wy}) &, coords];
      pad = OptionValue[BBoxPadCM];
      If[ListQ@triple && Length@triple == 3,
        {C, I, J} = (Intersection[V, Flatten@List@#] &) /@ triple,
        {C, I, J} = {{}, {}, {}}
      ];
      colorOff = Total[Length /@ {C, I, J}] == 0;
      idxMap = AssociationThread[V -> Range[Length@V]];
      nameOf[v_] := "v" <> IntegerString@Lookup[idxMap, v];
      borderTriple = OptionValue[BorderTriple];
      {starSet, triSet, sqSet} =
        If[ListQ@borderTriple && Length@borderTriple == 3,
          (Intersection[V, Flatten@List@#] &) /@ borderTriple,
          {{}, {}, {}}
        ];
      scheme = OptionValue[ColorScheme];
      styles =
        Switch[scheme,
          "ABC", {"Anode", "Bnode", "CnodeABC"},
          _, {"Cnode", "Inode", "Jnode"}
        ];
      baseTag[v_] :=
        If[colorOff,
          "Knode",
          Which[MemberQ[C, v], styles[[1]], MemberQ[I, v], styles[[2]], MemberQ[J, v], styles[[3]], True,
            "Knode"]
        ];
      borderTag[v_] :=
        Which[
          MemberQ[starSet, v], "BorderStar",
          MemberQ[triSet, v], "BorderTriangle",
          MemberQ[sqSet, v], "BorderSquare",
          True, ""
        ];
      tag[v_] := StringRiffle[DeleteCases[{baseTag[v], borderTag[v]}, ""], ","];
      evalData = OptionValue[EvalData];
      showEval = TrueQ@OptionValue[ShowEvalLabels];
      off = OptionValue[EvalLabelOffsetCM];
      headerOff = OptionValue[EvalHeaderOffsetCM];
      graphLevelQ = ListQ@evalData && Length@evalData == 3 && VectorQ[evalData, NumericQ];
      evAssoc = AssociationThread[V -> ConstantArray[None, Length@V]];
      If[! graphLevelQ,
        Which[
          AssociationQ@evalData && AllTrue[Keys@evalData, MemberQ[V, #] &] &&
            AllTrue[Values@evalData, (ListQ[#] && Length@# == 3 && VectorQ[#, NumericQ]) &],
          evAssoc = KeyTake[evalData, V],
          ListQ@evalData && Length@evalData == Length@V &&
            And @@ Map[(ListQ[#] && Length@# == 3 && VectorQ[#, NumericQ]) &, evalData],
          evAssoc = AssociationThread[V, evalData]
        ]
      ];
      extraTop = If[graphLevelQ && showEval, headerOff, 0];
      bbox =
        "\\path[use as bounding box] (" <> fmtNum[-pad] <> "," <> fmtNum[-pad] <>
          ") rectangle (" <> fmtNum[wx + pad] <> "," <> fmtNum[wy + pad + extraTop] <> ");";
      header =
        If[graphLevelQ && showEval,
          "\\node[anchor=south, inner sep=0pt] at (" <> fmtNum[wx/2] <> "," <> fmtNum[wy + headerOff] <>
            ") {" <> OptionValue[EvalLabelSize] <> " " <> fmtVec[evalData] <> "};",
          ""
        ];
      nodes =
        StringRiffle[
          Table[
            With[{v = vv, nm = nameOf[vv], xy = coords[vv], sty = tag[vv],
              lab =
                If[TrueQ@OptionValue[ShowLabels],
                  "{" <> OptionValue[LabelSize] <> " " <> IntegerString@Lookup[idxMap, vv] <> "}",
                  "{}"
                ],
              ev = Lookup[evAssoc, vv, None]},
              Module[{base, evalLine = ""},
                base = "\\node[" <> sty <> "] (" <> nm <> ") at (" <> fmtNum@xy[[1]] <> "," <>
                  fmtNum@xy[[2]] <> ") " <> lab <> ";";
                If[showEval && ListQ@ev && ! graphLevelQ,
                  evalLine =
                    "\\node[anchor=south, inner sep=0pt] at (" <>
                      fmtNum@(xy[[1]] + off[[1]]) <> "," <>
                      fmtNum@(xy[[2]] + off[[2]]) <> ") {" <>
                      OptionValue[EvalLabelSize] <> " " <> fmtVec[ev] <> "};";
                ];
                StringRiffle[DeleteCases[{base, evalLine}, ""], "\n"]
              ]
            ],
            {vv, V}
          ],
          "\n"
        ];
      nodes = If[header === "", nodes, header <> "\n" <> nodes];
      edges =
        StringRiffle[
          (With[{a = nameOf@#[[1]], b = nameOf@#[[2]]},
              "\\draw[" <> OptionValue[EdgeStyle] <> "] (" <> a <> ") -- (" <> b <> ");"
            ]) & /@ pairs,
          "\n"
        ];
      toTikZ[nodes, edges, bbox]
    ];

  Options[ToTikZAdjTripleBW] =
    Options[ToTikZAdjTripleColor] /. ("LabelSize" -> _) -> ("LabelSize" -> "\\scriptsize");

  ToTikZAdjTripleBW[adj_?MatrixQ, triple_, opts : OptionsPattern[]] :=
    Module[{out = ToTikZAdjTripleColor[adj, triple, opts]},
      StringReplace[
        out,
        {
          "Anode" -> "AnodeBW",
          "Bnode" -> "BnodeBW",
          "Cnode" -> "CnodeBW",
          "Inode" -> "InodeBW",
          "Jnode" -> "JnodeBW",
          "Knode" -> "KnodeBW"
        }
      ]
    ];

  Options[TikZTableFromAdjTriples] = {Columns -> 6, Version -> "Color",
    CanvasCM -> {4.5, 4.5}, Parallel -> False, IndexWhat -> "None",
    IndexStyle -> "\\footnotesize", IndexVspace -> "-2pt",
    PanelNumbering -> True, PanelStyle -> "\\Large\\bfseries",
    PanelPosition -> "below", PanelVspace -> "4pt", FitToWidth -> True,
    FixedColumnWidth -> True, CellScale -> 0.96, BBoxPadCM -> 0.10,
    RowVspace -> "6pt", RowSep -> "0pt", MultiPage -> False,
    EvalLabelSize -> "\\scriptsize", EvalLabelOffsetCM -> {0, 0.18},
    EvalHeaderOffsetCM -> 0.18, ColorScheme -> "CIJ"};

  TikZTableFromAdjTriples[pairs_List, OptionsPattern[]] :=
    Module[{cols = OptionValue[Columns], mapf, maker, cells, rows, colspec, widthExpr, resizeTikZ,
      makeOne, wrapMinipage},
      maker =
        If[OptionValue[Version] === "BW",
          Function[{adj, tri, borders, ev},
            ToTikZAdjTripleBW[
              adj, tri,
              CanvasCM -> OptionValue[CanvasCM],
              BBoxPadCM -> OptionValue[BBoxPadCM],
              EvalData -> ev,
              EvalLabelSize -> OptionValue[EvalLabelSize],
              EvalLabelOffsetCM -> OptionValue[EvalLabelOffsetCM],
              EvalHeaderOffsetCM -> OptionValue[EvalHeaderOffsetCM],
              ColorScheme -> OptionValue[ColorScheme],
              BorderTriple -> borders
            ]
          ],
          Function[{adj, tri, borders, ev},
            ToTikZAdjTripleColor[
              adj, tri,
              CanvasCM -> OptionValue[CanvasCM],
              BBoxPadCM -> OptionValue[BBoxPadCM],
              EvalData -> ev,
              EvalLabelSize -> OptionValue[EvalLabelSize],
              EvalLabelOffsetCM -> OptionValue[EvalLabelOffsetCM],
              EvalHeaderOffsetCM -> OptionValue[EvalHeaderOffsetCM],
              ColorScheme -> OptionValue[ColorScheme],
              BorderTriple -> borders
            ]
          ]
        ];
      mapf = If[TrueQ@OptionValue[Parallel], ParallelMap, Map];
      If[TrueQ@OptionValue[FixedColumnWidth],
        widthExpr = "\\dimexpr\\linewidth/" <> ToString[cols] <> "-2\\tabcolsep\\relax";
        colspec =
          StringRiffle[
            Table[">{\\centering\\arraybackslash}m{" <> widthExpr <> "}", {cols}],
            ""
          ],
        colspec = StringJoin@Table["c", {cols}]
      ];
      resizeTikZ[pic_String] :=
        If[
          TrueQ@OptionValue[FitToWidth],
          "\\resizebox{" <>
            ToString@NumberForm[OptionValue[CellScale], {Infinity, 3}, NumberPoint -> "."] <>
            "\\linewidth}{!}{" <> pic <> "}",
          pic
        ];
      wrapMinipage[body_String] :=
        "\\begin{minipage}[t]{\\linewidth}\\centering\n" <>
          If[OptionValue[RowVspace] =!= "0pt", "\\vspace{" <> OptionValue[RowVspace] <> "}\n", ""] <>
          body <>
          If[OptionValue[RowVspace] =!= "0pt", "\n\\vspace{" <> OptionValue[RowVspace] <> "}", ""] <>
          "\n\\end{minipage}";
      makeOne[item_, idx_Integer] :=
        Module[{adj, tri, borders = None, ev = None, pic, idxs = {}, idxline = "", panel = "", top = "",
          mid = "", bot = ""},
          Switch[Length[item],
            4, {adj, tri, borders, ev} = item,
            3, {adj, tri, ev} = item,
            2, {adj, tri} = item,
            _, {adj, tri} = item[[;; 2]]
          ];
          pic = maker[adj, tri, borders, ev];
          Which[
            OptionValue[IndexWhat] === "All", idxs = Range[Length@adj],
            OptionValue[IndexWhat] === "Triple", idxs = Flatten@List@tri,
            True, idxs = {}
          ];
          If[idxs =!= {},
            idxline =
              "{" <> OptionValue[IndexStyle] <> " " <>
                StringRiffle[("(" <> ToString[#] <> ")") & /@ idxs, "\\, "] <> "}"
          ];
          If[TrueQ@OptionValue[PanelNumbering],
            panel = "{" <> OptionValue[PanelStyle] <> " (" <> ToString[idx] <> ")}"
          ];
          mid = resizeTikZ[pic];
          If[OptionValue[PanelPosition] === "above",
            If[panel =!= "", top = panel <> "\\\\[" <> OptionValue[PanelVspace] <> "]\n"],
            If[panel =!= "", bot = bot <> "\\\\[" <> OptionValue[PanelVspace] <> "]\n" <> panel]
          ];
          If[idxline =!= "", bot = bot <> "\\\\[" <> OptionValue[IndexVspace] <> "]\n" <> idxline];
          wrapMinipage[top <> mid <> bot]
        ];
      cells = MapIndexed[makeOne[#1, #2[[1]]] &, pairs];
      rows = Partition[PadRight[cells, Ceiling[Length@cells/cols]*cols, ""], cols];
      Module[{rowSep = OptionValue[RowSep], n = Length@rows, beginEnv, endEnv},
        If[TrueQ@OptionValue[MultiPage],
          beginEnv =
            "\\begin{longtable}{" <> colspec <> "}\n" <>
              "\\endfirsthead\n\\endhead\n\\endfoot\n\\endlastfoot\n";
          endEnv = "\n\\end{longtable}",
          beginEnv = "\\begin{tabular}{" <> colspec <> "}\n";
          endEnv = "\n\\end{tabular}"
        ];
        beginEnv <>
          StringRiffle[
            MapIndexed[
              Function[{r, k},
                StringRiffle[r /. "" -> " ", " & "] <> " \\\\" <>
                  If[k[[1]] < n && rowSep =!= "0pt", "[" <> rowSep <> "]", ""]
              ],
              rows
            ],
            "\n"
          ] <> endEnv
      ]
    ];

  End[];

  EndPackage[];
