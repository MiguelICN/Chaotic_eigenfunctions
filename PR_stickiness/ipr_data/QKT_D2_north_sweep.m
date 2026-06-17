(* ::Package:: *)

(* ::Title::Closed:: *)
(*QKT D\:2082 Plot*)


(* ::Chapter::Closed:: *)
(*Classical Kicked Top*)


(* ::Input:: *)
(*ClassicalToSphere = Compile[{{q, _Real, 1}},*)
(*  Module[{b2 = q[[1]]^2 + q[[2]]^2},*)
(*    {q[[1]] Sqrt[1 - 0.25 b2], q[[2]] Sqrt[1 - 0.25 b2], 0.5 b2 - 1.0}*)
(*  ], CompilationTarget -> "C"];*)
(*(* disk center (r=0) -> sz=+1 = north pole; border (r=2) -> sz=-1 = south pole *)*)
(**)
(*ClassicalToDisk = Compile[{{s, _Real, 1}},*)
(*  Module[{b = Sqrt[2.0 / (1.0 - s[[3]])]},*)
(*    {s[[1]] b, s[[2]] b}*)
(*  ], CompilationTarget -> "C"];*)
(*(* inverse of ClassicalToSphere; singular at south pole (sz=-1) *)*)
(**)
(*ClassicalStep = Compile[*)
(*  {{s, _Real, 1}, {alpha, _Real}, {k, _Real}, {order, _Integer}},*)
(*  Module[{sx = s[[1]], sy = s[[2]], sz = s[[3]],*)
(*          ca = Cos[alpha], sa = Sin[alpha], sxc, sxs, r1, r2, r3},*)
(*    If[order == 0,*)
(*      (* rotation-kick *)*)
(*      r1 = ca sx - sa sy; r2 = sa sx + ca sy; r3 = sz;*)
(*      sx = r1; sy = r2; sz = r3;*)
(*      sxc = Cos[k sx]; sxs = Sin[k sx];*)
(*      r1 = sx; r2 = sxc sy - sxs sz; r3 = sxs sy + sxc sz;,*)
(*      (* kick-rotation *)*)
(*      sxc = Cos[k sx]; sxs = Sin[k sx];*)
(*      r1 = sx; r2 = sxc sy - sxs sz; r3 = sxs sy + sxc sz;*)
(*      sx = r1; sy = r2; sz = r3;*)
(*      r1 = ca sx - sa sy; r2 = sa sx + ca sy; r3 = sz;*)
(*    ];*)
(*    {r1, r2, r3}*)
(*  ], CompilationTarget -> "C"];*)
(**)
(*ClassicalMapeo = Compile[*)
(*  {{qini, _Real, 1}, {alpha, _Real}, {k, _Real},*)
(*   {n, _Integer}, {order, _Integer}},*)
(*  Module[{s, traj},*)
(*    s = ClassicalToSphere[qini];*)
(*    traj = Table[*)
(*      s = ClassicalStep[s, alpha, k, order];*)
(*      ClassicalToDisk[s], {n}];*)
(*    traj*)
(*  ],*)
(*  CompilationTarget -> "C",*)
(*  CompilationOptions -> {"InlineCompiledFunctions" -> True},*)
(*  RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"*)
(*];*)


(* ::Chapter::Closed:: *)
(*Spin Operators*)


(* ::Input:: *)
(*ClearAll[sz, sp, sm, sx, sx2];*)
(*(* Sz: diagonal {+j, j-1, ..., -j} *)*)
(*sz[j_] := sz[j] = SparseArray[Band[{1, 1}] -> Range[j, -j, -1]];*)
(*(* S+: superdiagonal; element (i, i+1) = <m_i|S+|m_{i+1}> with m_i = j-i+1 *)*)
(*sp[j_] := sp[j] = SparseArray[Band[{1, 2}] ->*)
(*    Table[Sqrt[j (j + 1) - m (m + 1)] // N, {m, j - 1, -j, -1}],*)
(*    {2 j + 1, 2 j + 1}];*)
(*(* S-: subdiagonal; transpose of S+ *)*)
(*sm[j_] := sm[j] = SparseArray[Band[{2, 1}] ->*)
(*    Table[Sqrt[j (j + 1) - m (m - 1)] // N, {m, j, -j + 1, -1}],*)
(*    {2 j + 1, 2 j + 1}];*)
(*sx[j_]  := sx[j]  = (1/2) (sp[j] + sm[j]);*)
(*sx2[j_] := sx2[j] = sx[j] . sx[j];*)


(* ::Chapter::Closed:: *)
(*Spin Coherent State Compiler *)


(* ::Input:: *)
(*Clear[generateCoherentStateCompiler];*)
(*generateCoherentStateCompiler = Compile[{{J, _Integer}, {q, _Real}, {p, _Real}},*)
(*  Module[{dim, zeta, result, mVals, logs, maxLog, terms, norm, r2},*)
(*    dim = 2 J + 1;*)
(*    result = ConstantArray[0.0 + 0. I, dim];*)
(*    r2 = q^2 + p^2;*)
(*    If[r2 >= 4.0,*)
(*      result[[1]] = 1.0 + 0. I,            (* north pole: m=+J, index 1, at border *)*)
(*      zeta = (q - I p) / Sqrt[4 - r2];    (* conjugate: south-pole param -> north at center *)*)
(*      If[Abs[zeta] == 0.0,*)
(*        result[[dim]] = 1.0 + 0. I,        (* south pole: m=-J, index dim, at center *)*)
(*        mVals = Range[J, -J, -1];          (* descending: J, J-1, ..., -J *)*)
(*        logs = 0.5 (LogGamma[2 J + 1] - LogGamma[J + 1 + mVals] -*)
(*                    LogGamma[J + 1 - mVals]) + (J + mVals) Log[zeta];*)
(*        maxLog = Max[Re[logs]];*)
(*        terms  = Exp[logs - maxLog];*)
(*        norm   = Sqrt[Total[Abs[terms]^2]];*)
(*        result = terms / norm;*)
(*      ]*)
(*    ];*)
(*    result*)
(*  ],*)
(*  CompilationTarget -> "WVM", RuntimeOptions -> "Speed",*)
(*  Parallelization -> True, RuntimeAttributes -> {Listable}*)
(*];*)


(* ::Chapter::Closed:: *)
(*Floquet Operator*)


(* ::Input:: *)
(*kickPart[j_, k_]    := MatrixExp[(-I k / (2 j)) sx2[j]];*)
(*freePart[j_, alpha_] := SparseArray[Band[{1, 1}] -> Exp[-I alpha Range[j, -j, -1]]];*)
(*Floq[j_, alpha_, k_] := kickPart[j, k] . freePart[j, alpha];*)
(*Floqsym[j_, alpha_, k_] :=freePart[j, alpha/2] . kickPart[j, k] . freePart[j, alpha/2];*)


(* ::Subsection:: *)
(*Parity Decomposition*)


(* ::Input:: *)
(*(* Sx^2 couples m to m +/- 2, so alternating array indices always decouple.*)
(*   Odd integer j: index 1 has m=j (odd), so Even-m sector starts at index 2.*)
(*   Even integer j or half-integer j: index 1 starts the first sector. *)*)
(*ParitySectorIndices[j_] :=*)
(*  Module[{d = 2 j + 1},*)
(*    If[IntegerQ[j] && OddQ[j],*)
(*      <|"Even" -> Range[2, d, 2], "Odd" -> Range[1, d, 2]|>,*)
(*      <|"Even" -> Range[1, d, 2], "Odd" -> Range[2, d, 2]|>*)
(*    ]*)
(*  ];*)


(* ::Input:: *)
(*GetQKTParityEigensystem[jParam_, alphaParam_, kParam_] :=*)
(*  Module[{mat, d, idx, matEven, matOdd, esEven, esOdd,*)
(*          phEven, phOdd, vecsEven, vecsOdd, idxE, idxO,*)
(*          fullVecsEven, fullVecsOdd},*)
(*mat = Floq[jParam, alphaParam, kParam];*)
(*d = 2 jParam + 1;*)
(*idx = ParitySectorIndices[jParam];*)
(*matEven = mat[[idx["Even"], idx["Even"]]];*)
(*matOdd  = mat[[idx["Odd"], idx["Odd"]]];*)
(*esEven = Eigensystem[matEven, Method -> "Direct"];*)
(*esOdd  = Eigensystem[matOdd, Method -> "Direct"];*)
(*phEven = Mod[Arg[esEven[[1]]], 2 Pi];*)
(*phOdd  = Mod[Arg[esOdd[[1]]], 2 Pi];*)
(*idxE = Ordering[phEven]; *)
(*idxO = Ordering[phOdd];*)
(*vecsEven = esEven[[2, idxE]];*)
(*vecsOdd  = esOdd[[2, idxO]];*)
(*    *)
(*fullVecsEven = Table[*)
(*Module[{vec = ConstantArray[0. + 0. I, d]},*)
(*        vec[[idx["Even"]]] = vecsEven[[i]]; vec],*)
(*      {i, Length[idx["Even"]]}];*)
(**)
(*    fullVecsOdd = Table[*)
(*Module[{vec = ConstantArray[0. + 0. I, d]},*)
(*        vec[[idx["Odd"]]] = vecsOdd[[i]]; vec],*)
(*      {i, Length[idx["Odd"]]}];*)
(**)
(*    <|"Even" -> <|"Quasienergies" -> esEven[[1]],*)
(*"FullVectors" -> fullVecsEven|>,*)
(*      "Odd"  -> <|"Quasienergies" -> esOdd[[1]],*)
(*"FullVectors" -> fullVecsOdd|>|>*)
(*  ];*)


(* ::Chapter::Closed:: *)
(*Parameters*)


(* ::Input:: *)
(*\[Alpha] = 0.84;*)
(*k = 3.0;*)
(*(* Quantum *)*)
(*J = 1000;*)
(*d = 2 J + 1;*)
(*(* Phase-space grid *)*)
(*radius = 2;*)
(*delta  = 0.020;*)


(* ::Input:: *)
(*(* Classical Poincare *)*)
(*kicks  = 2000;*)
(*points = 600;*)


(* ::Chapter::Closed:: *)
(*Classical Poincare Section*)


(* ::Input:: *)
(*randomregion  = RandomPoint[Disk[{0, 0}, 1.99], points];*)
(*classicalTraj = ClassicalMapeo[randomregion, \[Alpha], k, kicks, 0];*)
(*poincarePoints = ArrayReshape[classicalTraj, {points kicks, 2}];*)
(**)
(*poincare = ListPlot[poincarePoints,*)
(*  PlotStyle  -> Directive[PointSize[0.0005], Black, Opacity[0.15]],*)
(*  Frame      -> True, Axes -> True,*)
(*  GridLines  -> Automatic,*)
(*  PlotRange  -> {{-2, 2}, {-2, 2}},*)
(*  FrameLabel -> {Style["Q", 30, Black], Rotate[Style["P", 30, Black], -90 Degree]},*)
(*  LabelStyle -> Directive[30, Black],*)
(*  AspectRatio -> 1,*)
(*  ImageSize  -> 1000];*)


(* ::Chapter::Closed:: *)
(*Grid Over Disk*)


(* ::Input:: *)
(*gridPointsRaw = Flatten[Table[{x, y},*)
(*  {x, Range[-radius, radius, delta]},*)
(*  {y, Range[-radius, radius, delta]}], 1];*)
(*region = Select[gridPointsRaw, Norm[#] < radius &];*)
(*Print["Grid points: ", Length[region]];*)


(* ::Chapter::Closed:: *)
(*Coherent States*)


(* ::Input:: *)
(*LaunchKernels[12];*)


(* ::Input:: *)
(*Print["Generating coherent states over ", Length[region], " points..."];*)
(*statesCompiled = ParallelMap[*)
(*  generateCoherentStateCompiler[J, #[[1]], #[[2]]] &,*)
(*  region, Method -> "CoarsestGrained"];*)


(* ::Chapter::Closed:: *)
(*Floquet Eigensystem \[LongDash] Parity Decomposition*)


(* ::Input:: *)
(*(* Similarity transform: mat = exp(-i alpha Sz) . Floq . exp(i alpha Sz)*)
(*   = exp(-i alpha Sz) . exp(-ik/(2J) Sx^2)*)
(*   Sx^2 couples m to m+-2 only, so odd/even indices decouple in mat. *)*)
(*Print["Building Floquet matrix..."];*)
(*data=GetQKTParityEigensystem[J,\[Alpha],k];*)


(* ::Input:: *)
(*eigenvecM=data["Even","FullVectors"];*)
(*eigenvecm=data["Odd","FullVectors"];*)
(*neweigenvec=Riffle[eigenvecM,eigenvecm];*)


(* ::Chapter::Closed:: *)
(*Overlaps and D2*)


(* ::Input:: *)
(*Print["Projecting coherent states into Floquet basis..."];*)
(*overlaps  = Conjugate[neweigenvec] . Transpose[statesCompiled];*)
(**)
(*Print["Computing IPR..."];*)
(*iprValues = Total[Abs[overlaps]^4, {1}];*)
(**)
(*(* D2 = -log(IPR) / log(d),  d = 2J+1 *)*)
(*d2Values  = -Log[iprValues] / Log[d];*)


(* ::Chapter::Closed:: *)
(*D2 Density Plot + Poincare Overlay*)


(* ::Input:: *)
(*d2plot = ListDensityPlot[*)
(*  Transpose[{region[[All, 1]], region[[All, 2]], d2Values}],*)
(*  PlotRange     -> All,*)
(*  PlotTheme     -> "Detailed",*)
(*  ColorFunction -> "SunsetColors",*)
(*  PlotLabel     -> Style[*)
(*    "D\:2082, J = " <> ToString[J] <>*)
(*    ", \[Alpha] = " <> ToString[N[\[Alpha]]] <>*)
(*    ", k = " <> ToString[NumberForm[k, {4, 2}]], 25, Black],*)
(*  ImageSize     -> 1000,*)
(*  FrameLabel    -> {Style["Q", 40, Black], Rotate[Style["P", 40, Black], -90 Degree]},*)
(*  LabelStyle    -> Directive[30, Black],*)
(*  AspectRatio   -> 1,*)
(*  Frame         -> True,*)
(*  InterpolationOrder -> 3];*)
(**)
(*overlay = Show[{d2plot, poincare}];*)
(**)
(*Export[*)
(*  "north_D2_J_" <> ToString[J] <>*)
(*  "_alpha_" <> StringReplace[ToString[N[\[Alpha]]], "." -> "p"] <>*)
(*  "_k_" <> StringReplace[ToString[k], "." -> "p"] <> ".png",*)
(*  overlay, ImageResolution -> 200]*)


(* ::Title::Closed:: *)
(*IPR statistics*)


(* ::Chapter::Closed:: *)
(*Floquet Operator*)


(* ::Subsection:: *)
(*Parity Decomposition*)


(* ::Input:: *)
(*GetQKTParityEigensystemsym[jParam_, alphaParam_, kParam_] :=*)
(*  Module[{mat, d, idx, matEven, matOdd, esEven, esOdd,*)
(*          phEven, phOdd, vecsEven, vecsOdd, idxE, idxO,*)
(*          fullVecsEven, fullVecsOdd},*)
(**)
(*mat = Floqsym[jParam, alphaParam, kParam];*)
(*d = 2 jParam + 1;*)
(*idx = ParitySectorIndices[jParam];*)
(*matEven = mat[[idx["Even"], idx["Even"]]];*)
(*matOdd  = mat[[idx["Odd"], idx["Odd"]]];*)
(*esEven = Eigensystem[matEven, Method -> "Direct"];*)
(*esOdd  = Eigensystem[matOdd, Method -> "Direct"];*)
(*phEven = Mod[Arg[esEven[[1]]], 2 Pi];*)
(*phOdd  = Mod[Arg[esOdd[[1]]], 2 Pi];*)
(*idxE = Ordering[phEven]; *)
(*idxO = Ordering[phOdd];*)
(*vecsEven = esEven[[2, idxE]];*)
(*vecsOdd  = esOdd[[2, idxO]];*)
(*    *)
(*fullVecsEven = Table[*)
(*Module[{vec = ConstantArray[0. + 0. I, d]},*)
(*        vec[[idx["Even"]]] = vecsEven[[i]]; vec],*)
(*      {i, Length[idx["Even"]]}];*)
(**)
(*    fullVecsOdd = Table[*)
(*Module[{vec = ConstantArray[0. + 0. I, d]},*)
(*        vec[[idx["Odd"]]] = vecsOdd[[i]]; vec],*)
(*      {i, Length[idx["Odd"]]}];*)
(**)
(*    <|"Even" -> <|"Quasienergies" -> esEven[[1]],*)
(*"FullVectors" -> fullVecsEven|>,*)
(*      "Odd"  -> <|"Quasienergies" -> esOdd[[1]],*)
(*"FullVectors" -> fullVecsOdd|>|>*)
(*  ]*)


(* ::Chapter::Closed:: *)
(*Parameters*)


(* ::Input:: *)
(*\[Alpha] = 0.84;*)
(*k = 1.3;*)
(*(* Quantum *)*)
(*J = 1000;*)
(*d = 2 J + 1;*)
(*(* Phase-space grid *)*)
(*radius = 2;*)
(*delta  = 0.020;*)


(* ::Chapter::Closed:: *)
(*Grid Over Disk*)


(* ::Input:: *)
(*gridPointsRaw = Flatten[Table[{x, y},*)
(*  {x, Range[-radius, radius, delta]},*)
(*  {y, Range[-radius, radius, delta]}], 1];*)
(*region = Select[gridPointsRaw, Norm[#] < radius &];*)
(*Print["Grid points: ", Length[region]];*)


(* ::Chapter::Closed:: *)
(*Coherent States*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]];*)
(*LaunchKernels[12];*)


(* ::Input:: *)
(*Print["Generating coherent states over ", Length[region], " points..."];*)
(*statesCompiled = ParallelMap[*)
(*  generateCoherentStateCompiler[J, #[[1]], #[[2]]] &,*)
(*  region, Method -> "CoarsestGrained"];*)


(* ::Chapter::Closed:: *)
(*Floquet Eigensystem \[LongDash] Parity Decomposition*)


(* ::Input:: *)
(*(* Similarity transform: mat = exp(-i alpha Sz) . Floq . exp(i alpha Sz)*)
(*   = exp(-i alpha Sz) . exp(-ik/(2J) Sx^2)*)
(*   Sx^2 couples m to m+-2 only, so odd/even indices decouple in mat. *)*)
(*Print["Building Floquet matrix..."];*)
(*data=GetQKTParityEigensystem[J,\[Alpha],k];*)
(*datasym=GetQKTParityEigensystemsym[J,\[Alpha],k];*)


(* ::Input:: *)
(*eigenvecM=data["Even","FullVectors"];*)
(*eigenvecm=data["Odd","FullVectors"];*)
(*neweigenvec=Riffle[eigenvecM,eigenvecm];*)
(*eigenvecMsym=datasym["Even","FullVectors"];*)
(*eigenvecmsym=datasym["Odd","FullVectors"];*)
(*neweigenvecsym=Riffle[eigenvecMsym,eigenvecmsym];*)


(* ::Chapter::Closed:: *)
(*Overlaps and D2*)


(* ::Input:: *)
(*Print["Projecting coherent states into Floquet basis..."];*)
(*overlaps  = Conjugate[neweigenvec] . Transpose[statesCompiled];*)
(*overlapsym  = Conjugate[neweigenvecsym] . Transpose[statesCompiled];*)
(**)
(*Print["Computing IPR..."];*)
(*iprValues = Total[Abs[overlaps]^4, {1}];*)
(*iprValuesym = Total[Abs[overlapsym]^4, {1}];*)
(**)
(*(* D2 = -log(IPR) / log(d),  d = 2J+1 *)*)
(*d2Values  = -Log[iprValues] / Log[d];*)
(*d2Valuesym  = -Log[iprValuesym] / Log[d];*)


data=Import["C:\\Users\\Miguel\\Github\\Chaotic_eigenfunctions\\PR_stickiness\\ipr_J_1000_alpha_0p84_k_0p5_nonsym.mx"];
datasym=Import["C:\\Users\\Miguel\\Github\\Chaotic_eigenfunctions\\PR_stickiness\\ipr_J_1000_alpha_0p84_k_0p5_sym.mx"];


iprValues=data;
iprValuesym=datasym;


J=1000;
\[Alpha]=0.84;
k=0.5;
d=2J+1;


(* ::Input:: *)
(*Legended[Histogram[{(1/iprValues)(1/d),(1/iprValuesym)(1/d)},*)
(*50,*)
(*PlotTheme->"Detailed",*)
(*PlotLabel     -> Style[*)
(*    "PR, J = " <> ToString[J] <>*)
(*    ", \[Alpha] = " <> ToString[N[\[Alpha]]] <>*)
(*    ", k = " <> ToString[NumberForm[k, {4, 2}]], 25, Black],*)
(*FrameLabel    -> {Style["PR", 40, Black], Rotate[Style["P(PR)", 40, Black], -90 Degree]},*)
(*  LabelStyle    -> Directive[30, Black],*)
(*ImageSize->1000,PlotRange->All],*)
(*Placed[LineLegend[{ColorData[97][1],ColorData[97][2]},{"Non-symmetric","Symmetrized"}],{0.7,0.7}]]*)


(* ::Input:: *)
(*d2plotsym= ListDensityPlot[*)
(*  Transpose[{region[[All, 1]], region[[All, 2]], d2Valuesym}],*)
(*  PlotRange     -> All,*)
(*  PlotTheme     -> "Detailed",*)
(*  ColorFunction -> "SunsetColors",*)
(*  PlotLabel     -> Style[*)
(*    "Symmetrized D\:2082, J = " <> ToString[J] <>*)
(*    ", \[Alpha] = " <> ToString[N[\[Alpha]]] <>*)
(*    ", k = " <> ToString[NumberForm[k, {4, 2}]], 25, Black],*)
(*  ImageSize     -> 1000,*)
(*  FrameLabel    -> {Style["Q", 40, Black], Rotate[Style["P", 40, Black], -90 Degree]},*)
(*  LabelStyle    -> Directive[30, Black],*)
(*  AspectRatio   -> 1,*)
(*  Frame         -> True,*)
(*  InterpolationOrder -> 3];*)


(* ::Input:: *)
(*Export[*)
(*  "north_sym_D2_J_" <> ToString[J] <>*)
(*  "_alpha_" <> StringReplace[ToString[N[\[Alpha]]], "." -> "p"] <>*)
(*  "_k_" <> StringReplace[ToString[k], "." -> "p"] <> ".png",*)
(*  d2plotsym, ImageResolution -> 300]*)


(* ::Title:: *)
(*K-Sweep*)


(* ::Chapter:: *)
(*Parameters*)


(* ::Input:: *)
(*\[Alpha]KS  = 0.84;*)
(*JKS   = 1000;*)
(*dKS   = 2 JKS + 1;*)
(*kListKS  = Range[0.05,30.0,0.05];*)
(*radiusKS = 2;*)
(*deltaKS  = 0.015;*)


kListKS//Length


(* ::Chapter:: *)
(*Precomputation*)


(* ::Input:: *)
(*regionKS = Select[Flatten[Table[{x, y},*)
(*    {x, Range[-radiusKS, radiusKS, deltaKS]},*)
(*    {y, Range[-radiusKS, radiusKS, deltaKS]}], 1], Norm[#] < radiusKS &];*)
(*Export["region_J_" <> ToString[JKS] <> "_delta_" <>*)
(*  StringReplace[ToString[N[deltaKS]], "." -> "p"] <> ".mx", regionKS];*)
(*regionKS//Length*)


SetDirectory[NotebookDirectory[]];
LaunchKernels[12];


(* ::Input:: *)
(*statesKS = ParallelMap[generateCoherentStateCompiler[JKS, #[[1]], #[[2]]] &,*)
(*  regionKS, Method -> "CoarsestGrained"];*)


(* ::Input:: *)
(*idxKS        = ParitySectorIndices[JKS];*)
(*statesEvenKS = Transpose[statesKS][[idxKS["Even"], All]];*)
(*statesOddKS  = Transpose[statesKS][[idxKS["Odd"],  All]];*)
(*Clear[statesKS];*)


(* ::Input:: *)
(*{sx2EvalsEKS, sx2EvecsEKS} = Eigensystem[N @ Normal @ sx2[JKS][[idxKS["Even"], idxKS["Even"]]]];*)
(*{sx2EvalsOKS, sx2EvecsOKS} = Eigensystem[N @ Normal @ sx2[JKS][[idxKS["Odd"],  idxKS["Odd"]]]];*)
(*sx2EvecsEHKS = ConjugateTranspose[sx2EvecsEKS];*)
(*sx2EvecsOHKS = ConjugateTranspose[sx2EvecsOKS];*)


(* ::Input:: *)
(*freeVecKS  = Exp[-I \[Alpha]KS       Range[JKS, -JKS, -1]];*)
(*freeVecHKS = Exp[-I (\[Alpha]KS / 2) Range[JKS, -JKS, -1]];*)
(*freeEKS  = freeVecKS[[idxKS["Even"]]];  freeHEKS = freeVecHKS[[idxKS["Even"]]];*)
(*freeOKS  = freeVecKS[[idxKS["Odd"]]];   freeHOKS = freeVecHKS[[idxKS["Odd"]]];*)
(*Clear[freeVecKS, freeVecHKS];*)


(* ::Chapter::Closed:: *)
(*K-Sweep Loop*)


(* ::Input:: *)
(*Do[*)
(*  kickEKS = sx2EvecsEHKS . (Exp[(-I kVal / (2 JKS)) sx2EvalsEKS] * sx2EvecsEKS);*)
(*  kickOKS = sx2EvecsOHKS . (Exp[(-I kVal / (2 JKS)) sx2EvalsOKS] * sx2EvecsOKS);*)
(**)
(*  floqEKS    = Transpose[freeEKS  * Transpose[kickEKS]];*)
(*  floqOKS    = Transpose[freeOKS  * Transpose[kickOKS]];*)
(*  floqEsymKS = Transpose[freeHEKS * Transpose[freeHEKS * kickEKS]];*)
(*  floqOsymKS = Transpose[freeHOKS * Transpose[freeHOKS * kickOKS]];*)
(*  Clear[kickEKS, kickOKS];*)
(**)
(*  {esEkns, esOkns, esEksym, esOksym} = Parallelize[{*)
(*    Eigensystem[floqEKS,    Method -> "Direct"][[2]],*)
(*    Eigensystem[floqOKS,    Method -> "Direct"][[2]],*)
(*    Eigensystem[floqEsymKS, Method -> "Direct"][[2]],*)
(*    Eigensystem[floqOsymKS, Method -> "Direct"][[2]]*)
(*  }];*)
(*  Clear[floqEKS, floqOKS, floqEsymKS, floqOsymKS];*)
(**)
(*  iprEkns    = Total[Abs[Conjugate[esEkns]  . statesEvenKS]^4, {1}]; Clear[esEkns];*)
(*  iprOkns    = Total[Abs[Conjugate[esOkns]  . statesOddKS]^4,  {1}]; Clear[esOkns];*)
(*  iprEksym   = Total[Abs[Conjugate[esEksym] . statesEvenKS]^4, {1}]; Clear[esEksym];*)
(*  iprOksym   = Total[Abs[Conjugate[esOksym] . statesOddKS]^4,  {1}]; Clear[esOksym];*)
(**)
(*  iprKS    = iprEkns + iprOkns;   Clear[iprEkns, iprOkns];*)
(*  iprSymKS = iprEksym + iprOksym; Clear[iprEksym, iprOksym];*)
(**)
(*  nameKS = "_J_" <> ToString[JKS] <>*)
(*           "_alpha_" <> StringReplace[ToString[N[\[Alpha]KS]], "." -> "p"] <>*)
(*           "_k_"     <> StringReplace[ToString[N[kVal]],       "." -> "p"];*)
(*  Export["ipr" <> nameKS <> "_nonsym.mx", iprKS];*)
(*  Export["ipr" <> nameKS <> "_sym.mx",    iprSymKS];*)
(*  Clear[iprKS, iprSymKS];*)
(*  Print["k = ", kVal, " done."];*)
(**)
(*, {kVal, kListKS}]*)


(* ::Chapter:: *)
(*Plotting, dq = 0.020*)


\[Alpha]AV     = 0.84;
JAV     = 1000;
dAV     = 2 JAV + 1;
kListAV = Range[0.05, 30.0, 0.05];
deltaAV = 0.020;


nameAV[kVal_] := "_J_" <> ToString[JAV] <>
                 "_alpha_" <> StringReplace[ToString[N[\[Alpha]AV]], "." -> "p"] <>
                 "_k_"     <> StringReplace[ToString[N[kVal]], "." -> "p"];


regionAV = Import["region_J_" <> ToString[JAV] <> "_delta_" <>
  StringReplace[ToString[N[deltaAV]], "." -> "p"] <> ".mx"];


iprDataNS  = Association[# -> Import["ipr" <> nameAV[#] <> "_nonsym.mx"] & /@ kListAV];
iprDataSym = Association[# -> Import["ipr" <> nameAV[#] <> "_sym.mx"]    & /@ kListAV];


Histogram[{1/iprDataNS[[50]],1/iprDataSym[[50]]},100,PlotTheme->"Detailed"]


(* ::Input:: *)
(*meanD2NS  = {#, Mean[-Log[iprDataNS[#]]  / Log[dAV]]} & /@ kListAV;*)
(*meanD2Sym = {#, Mean[-Log[iprDataSym[#]] / Log[dAV]]} & /@ kListAV;*)


(* ::Input:: *)
(*meanD2NS  = {#, (1/dAV)Mean[1/iprDataNS[#]]} & /@ kListAV;*)
(*meanD2Sym = {#, (1/dAV)Mean[1/iprDataSym[#]]} & /@ kListAV;*)


(* ::Input:: *)
(*avgPlot = ListLinePlot[{meanD2NS},*)
(*  PlotMarkers -> {Automatic, 10},*)
(*  PlotStyle   -> {Directive[Blue, Thick], Directive[Red, Thick]},*)
(*  PlotLegends -> Placed[LineLegend[{"Non-sym", "Sym"},*)
(*    LegendMarkerSize -> 15, LegendFunction -> "Frame"], {0.9, 0.8}],*)
(*  PlotLabel   -> Style["\!\(\*OverscriptBox[\(D\), \(_\)]\)\:2082 vs k,  J = " <>*)
(*    ToString[JAV] <> ", \[Alpha] = " <> ToString[N[\[Alpha]AV]], 22, Black],*)
(*  FrameLabel  -> {Style["k", 30, Black], Rotate[Style["\!\(\*OverscriptBox[\(D\), \(_\)]\)\:2082", 30, Black], 90 Degree]},*)
(*  LabelStyle  -> Directive[20, Black],*)
(*Frame -> True, GridLines -> Automatic,*)
(*ImageSize -> 800,*)
(*PlotRange->All]*)


(* ::Input:: *)
(*avgPlot = ListLinePlot[{meanD2Sym},*)
(*  PlotMarkers -> {Automatic, 10},*)
(*  PlotStyle   -> {Directive[Red, Thick]},*)
(*  PlotLegends -> Placed[LineLegend[{"Non-sym", "Sym"},*)
(*    LegendMarkerSize -> 15, LegendFunction -> "Frame"], {0.9, 0.8}],*)
(*  PlotLabel   -> Style["\!\(\*OverscriptBox[\(D\), \(_\)]\)\:2082 vs k,  J = " <>*)
(*    ToString[JAV] <> ", \[Alpha] = " <> ToString[N[\[Alpha]AV]], 22, Black],*)
(*  FrameLabel  -> {Style["k", 30, Black], Rotate[Style["\!\(\*OverscriptBox[\(D\), \(_\)]\)\:2082", 30, Black], 90 Degree]},*)
(*  LabelStyle  -> Directive[20, Black],*)
(*Frame -> True, GridLines -> Automatic,*)
(*ImageSize -> 800,*)
(*PlotRange->All]*)


(* ::Input:: *)
(*avgPlot = ListLinePlot[{meanD2NS,meanD2Sym},*)
(*  PlotMarkers -> {Automatic, 6},*)
(*  PlotStyle   -> {Directive[Blue,Thick], Directive[Brown,Thin]},*)
(*  PlotLegends -> Placed[LineLegend[{"Non-sym", "Sym"},*)
(*    LegendMarkerSize -> 15, LegendFunction -> "Frame"], {0.9, 0.8}],*)
(*  PlotLabel   -> Style["\!\(\*OverscriptBox[\(D\), \(_\)]\)\:2082 vs k,  J = " <>*)
(*    ToString[JAV] <> ", \[Alpha] = " <> ToString[N[\[Alpha]AV]], 22, Black],*)
(*  FrameLabel  -> {Style["k", 30, Black], Rotate[Style["\!\(\*OverscriptBox[\(D\), \(_\)]\)\:2082", 30, Black],-90 Degree]},*)
(*  LabelStyle  -> Directive[20, Black],*)
(*Frame -> True, GridLines -> Automatic,*)
(*ImageSize -> 800,*)
(*PlotRange->{{0,16},{0,0.5}}]*)


(* ::Input:: *)
(*avgPlot = ListLinePlot[{meanD2NS,meanD2Sym},*)
(*  PlotMarkers -> {Automatic, 6},*)
(*  PlotStyle   -> {Directive[Blue,Thick], Directive[Brown,Thin]},*)
(*  PlotLegends -> Placed[LineLegend[{"Non-sym", "Sym"},*)
(*    LegendMarkerSize -> 15, LegendFunction -> "Frame"], {0.9, 0.9}],*)
(*  PlotLabel   -> Style["\!\(\*OverscriptBox[\(D\), \(_\)]\)\:2082 vs k,  J = " <>*)
(*    ToString[JAV] <> ", \[Alpha] = " <> ToString[N[\[Alpha]AV]], 22, Black],*)
(*  FrameLabel  -> {Style["k", 30, Black], Rotate[Style["\!\(\*OverscriptBox[\(D\), \(_\)]\)\:2082", 30, Black],-90 Degree]},*)
(*  LabelStyle  -> Directive[20, Black],*)
(*Frame -> True, GridLines -> Automatic,*)
(*ImageSize -> 1000,*)
(*PlotRange->{{10,20},{0.45,0.49}}]*)


(* ::Input:: *)
(*Export["meanD2_J_" <> ToString[JAV] <>*)
(*  "_alpha_" <> StringReplace[ToString[N[\[Alpha]AV]], "." -> "p"] <> ".png",*)
(*  avgPlot, ImageResolution -> 200]*)


avgPlot = ListLinePlot[{meanD2NS,meanD2Sym},
  PlotMarkers -> {Automatic, 6},
  PlotStyle   -> {Directive[Blue,Thick], Directive[Brown,Thin]},
  PlotLegends -> Placed[LineLegend[{"Non-sym", "Sym"},
    LegendMarkerSize -> 15, LegendFunction -> "Frame"], {0.9, 0.8}],
  PlotLabel   -> Style["PR vs k,  J = " <>
    ToString[JAV] <> ", \[Alpha] = " <> ToString[N[\[Alpha]AV]], 22, Black],
  FrameLabel  -> {Style["k", 30, Black], Rotate[Style["PR", 30, Black],-90 Degree]},
  LabelStyle  -> Directive[20, Black],
Frame -> True, GridLines -> Automatic,
ImageSize -> 800,
PlotRange->{{0,16},{0,0.5}}]
