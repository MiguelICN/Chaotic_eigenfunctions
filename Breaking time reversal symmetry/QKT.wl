(* ::Package:: *)

(* ::Package:: *)
(* QKT.wl \[LongDash] Quantum Kicked Top Efficient Spin Operators & Floquet Definitions *)


BeginPackage["QKT`"];

generateSpinOperators::usage = "generateSpinOperators[J] returns an association of {Sx, Sy, Sz, Sx2} spin operators for total spin J.";
Floq::usage = "Floq[J, \[Alpha], k, \[Tau]] returns the Floquet operator with twist along Sz.";
Floqn::usage = "Floqn[J, \[Alpha], k, nVec] returns the Floquet operator with twist along an arbitrary axis defined by nVec.";
generateCoherentStateCompiler::usage = 
  "generateCoherentStateCompiler[] returns a compiled function that generates spin coherent states.";

Begin["`Private`"];

ClearAll[sz, sp, sm, sx, sy, sx2];

(* Spin operators: memoized *)
sz[j_] := sz[j] = DiagonalMatrix[Table[-j + k - 1, {k, 1, 2 j + 1}]];
sp[j_] := sp[j] = DiagonalMatrix[Table[Sqrt[j (j + 1) - m (m + 1)], {m, -j, j - 1}], 1];
sm[j_] := sm[j] = DiagonalMatrix[Table[Sqrt[j (j + 1) - m (m - 1)], {m, -j + 1, j}], -1];
sx[j_] := sx[j] = SparseArray[(1/2) (sp[j] + sm[j])];
sy[j_] := sy[j] = SparseArray[(1/(2 I)) (sp[j] - sm[j])];
sx2[j_] := sx2[j] = sx[j] . sx[j];

(* Operator bundle *)
generateSpinOperators[j_] := <|
  "Sx" -> sx[j],
  "Sy" -> sy[j],
  "Sz" -> sz[j],
  "Sx2" -> sx2[j]
|>;

(* Floquet: original with twist along Sz *)
kickPart[j_, k_] := kickPart[j, k] = MatrixExp[(-I k/(2 j)) sx2[j]];
twistPart[j_, \[Alpha]_, \[Tau]_] := MatrixExp[(-I \[Alpha] \[Tau]) sz[j]];
Floq[j_, \[Alpha]_, k_, \[Tau]_] := kickPart[j, k] . twistPart[j, \[Alpha], \[Tau]];

(* Floquet: generalized twist along arbitrary axis *)
twistPartGeneral[j_, \[Alpha]_, nVec_] := twistPartGeneral[j, \[Alpha], nVec] = Module[
  {u = Normalize[nVec], gen},
  gen = u . {sx[j], sy[j], sz[j]};
  MatrixExp[-I \[Alpha] gen]
];

Floqn[j_, \[Alpha]_, k_, nVec_] := kickPart[j, k] . twistPartGeneral[j, \[Alpha], nVec];

(*Coherent states*)
generateCoherentStateCompiler[] := Module[{},
  Compile[{{J, _Integer}, {q, _Real}, {p, _Real}},
   Module[{dim, \[Alpha], result, mVals, logs, maxLog, terms, norm},
    dim = 2 J + 1;
    result = ConstantArray[0.0 + 0. I, dim];
    \[Alpha] = (q + I p)/Sqrt[4 - (q^2 + p^2)];
    If[Abs[\[Alpha]] == 0,
     result[[1]] = 1.0 + 0. I,
     mVals = Range[-J, J];
     logs = 0.5 (LogGamma[2 J + 1] - LogGamma[J + 1 + mVals] - 
         LogGamma[J + 1 - mVals]) + (J + mVals) Log[\[Alpha]];
     maxLog = Max[Re[logs]];
     terms = Exp[logs - maxLog];
     norm = Sqrt[Total[Abs[terms]^2]];
     result = terms/norm;
     ];
    result
    ],
   CompilationTarget -> "WVM",
   RuntimeOptions -> "Speed",
   Parallelization -> True,
   RuntimeAttributes -> {Listable}
   ]
];


End[];
EndPackage[];
