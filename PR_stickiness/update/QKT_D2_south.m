(* ::Package:: *)

(* ::Title:: *)
(*QKT D\:2082 Plot*)


(* ::Chapter:: *)
(*Classical Kicked Top*)


ClassicalToSphere = Compile[{{q, _Real, 1}},
  Module[{b2 = q[[1]]^2 + q[[2]]^2},
    {q[[1]] Sqrt[1 - 0.25 b2], q[[2]] Sqrt[1 - 0.25 b2], 1.0 - 0.5 b2}
  ], CompilationTarget -> "C"];
(* disk center (r=0) -> sz=+1 = north pole; border (r=2) -> sz=-1 = south pole *)

ClassicalToDisk = Compile[{{s, _Real, 1}},
  Module[{b = Sqrt[2.0 / (1.0 + s[[3]])]},
    {s[[1]] b, s[[2]] b}
  ], CompilationTarget -> "C"];
(* inverse of ClassicalToSphere; singular at south pole (sz=-1) *)

ClassicalStep = Compile[
  {{s, _Real, 1}, {alpha, _Real}, {k, _Real}, {order, _Integer}},
  Module[{sx = s[[1]], sy = s[[2]], sz = s[[3]],
          ca = Cos[alpha], sa = Sin[alpha], sxc, sxs, r1, r2, r3},
    If[order == 0,
      (* rotation-kick *)
      r1 = ca sx - sa sy; r2 = sa sx + ca sy; r3 = sz;
      sx = r1; sy = r2; sz = r3;
      sxc = Cos[k sx]; sxs = Sin[k sx];
      r1 = sx; r2 = sxc sy - sxs sz; r3 = sxs sy + sxc sz;,
      (* kick-rotation *)
      sxc = Cos[k sx]; sxs = Sin[k sx];
      r1 = sx; r2 = sxc sy - sxs sz; r3 = sxs sy + sxc sz;
      sx = r1; sy = r2; sz = r3;
      r1 = ca sx - sa sy; r2 = sa sx + ca sy; r3 = sz;
    ];
    {r1, r2, r3}
  ], CompilationTarget -> "C"];

ClassicalMapeo = Compile[
  {{qini, _Real, 1}, {alpha, _Real}, {k, _Real},
   {n, _Integer}, {order, _Integer}},
  Module[{s, traj},
    s = ClassicalToSphere[qini];
    traj = Table[
      s = ClassicalStep[s, alpha, k, order];
      ClassicalToDisk[s], {n}];
    traj
  ],
  CompilationTarget -> "C",
  CompilationOptions -> {"InlineCompiledFunctions" -> True},
  RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"
];


(* ::Chapter:: *)
(*Spin Operators*)


ClearAll[sz, sp, sm, sx, sx2];

(* Sz: diagonal {+j, j-1, ..., -j} *)
sz[j_] := sz[j] = SparseArray[Band[{1, 1}] -> Range[j, -j, -1]];

(* S+: superdiagonal; element (i, i+1) = <m_i|S+|m_{i+1}> with m_i = j-i+1 *)
sp[j_] := sp[j] = SparseArray[Band[{1, 2}] ->
    Table[Sqrt[j (j + 1) - m (m + 1)] // N, {m, j - 1, -j, -1}],
    {2 j + 1, 2 j + 1}];

(* S-: subdiagonal; transpose of S+ *)
sm[j_] := sm[j] = SparseArray[Band[{2, 1}] ->
    Table[Sqrt[j (j + 1) - m (m - 1)] // N, {m, j, -j + 1, -1}],
    {2 j + 1, 2 j + 1}];

sx[j_]  := sx[j]  = (1/2) (sp[j] + sm[j]);
sx2[j_] := sx2[j] = sx[j] . sx[j];


(* ::Chapter:: *)
(*Spin Coherent State Compiler *)


generateCoherentStateCompiler = Compile[{{J, _Integer}, {q, _Real}, {p, _Real}},
  Module[{dim, zeta, result, mVals, logs, maxLog, terms, norm, r2},
    dim = 2 J + 1;
    result = ConstantArray[0.0 + 0. I, dim];
    r2 = q^2 + p^2;
    If[r2 >= 4.0,
      result[[dim]] = 1.0 + 0. I,          (* south pole: m=-J, index dim *)
      zeta = (q + I p) / Sqrt[4 - r2];
      If[Abs[zeta] == 0.0,
        result[[1]] = 1.0 + 0. I,          (* north pole: m=+J, index 1  *)
        mVals = Range[J, -J, -1];          (* descending: J, J-1, ..., -J *)
        logs = 0.5 (LogGamma[2 J + 1] - LogGamma[J + 1 + mVals] -
                    LogGamma[J + 1 - mVals]) + (J - mVals) Log[zeta];
        maxLog = Max[Re[logs]];
        terms  = Exp[logs - maxLog];
        norm   = Sqrt[Total[Abs[terms]^2]];
        result = terms / norm;
      ]
    ];
    result
  ],
  CompilationTarget -> "C", RuntimeOptions -> "Speed",
  Parallelization -> True, RuntimeAttributes -> {Listable}
];


(* ::Chapter:: *)
(*Floquet Operator *)


kickPart[j_, k_]    := MatrixExp[(-I k / (2 j)) sx2[j]];
freePart[j_, alpha_] := SparseArray[Band[{1, 1}] -> Exp[-I alpha Range[j, -j, -1]]];

Floq[j_, alpha_, k_] := kickPart[j, k] . freePart[j, alpha];
Floqsym[j_, alpha_, k_] := freePart[j, alpha/2] . kickPart[j, k] . freePart[j, alpha/2];


(* ================================================================ *)
(* 2. PARAMETERS                                                     *)
(* ================================================================ *)


(* ::Chapter:: *)
(*Parameters*)


\[Alpha] = 0.84;
k = 3.0;

(* Classical Poincare *)
kicks  = 2000;
points = 600;

(* Quantum *)
J = 1000;
d = 2 J + 1;

(* Phase-space grid *)
radius = 2;
delta  = 0.020;


(* ================================================================ *)
(* 3. EXECUTION                                                      *)
(* ================================================================ *)


(* ::Chapter::Closed:: *)
(*Classical Poincare Section*)


randomregion  = RandomPoint[Disk[{0, 0}, 1.99], points];
classicalTraj = ClassicalMapeo[randomregion, \[Alpha], k, kicks, 0];
poincarePoints = ArrayReshape[classicalTraj, {points kicks, 2}];

poincare = ListPlot[poincarePoints,
  PlotStyle  -> Directive[PointSize[0.0005], Black, Opacity[0.15]],
  Frame      -> True, Axes -> True,
  GridLines  -> Automatic,
  PlotRange  -> {{-2, 2}, {-2, 2}},
  FrameLabel -> {Style["Q", 30, Black], Rotate[Style["P", 30, Black], -90 Degree]},
  LabelStyle -> Directive[30, Black],
  AspectRatio -> 1,
  ImageSize  -> 1000];


(* ::Chapter::Closed:: *)
(*Grid Over Disk*)


gridPointsRaw = Flatten[Table[{x, y},
  {x, Range[-radius, radius, delta]},
  {y, Range[-radius, radius, delta]}], 1];
region = Select[gridPointsRaw, Norm[#] < radius &];
Print["Grid points: ", Length[region]];


(* ::Chapter:: *)
(*Coherent States*)


LaunchKernels[12];


Print["Generating coherent states over ", Length[region], " points..."];
statesCompiled = ParallelMap[
  generateCoherentStateCompiler[J, #[[1]], #[[2]]] &,
  region, Method -> "CoarsestGrained"];


(* ::Chapter:: *)
(*Floquet Eigensystem \[LongDash] Parity Decomposition*)


(* Similarity transform: mat = exp(-i alpha Sz) . Floq . exp(i alpha Sz)
   = exp(-i alpha Sz) . exp(-ik/(2J) Sx^2)
   Sx^2 couples m to m+-2 only, so odd/even indices decouple in mat. *)

Print["Building Floquet matrix..."];
mat = Floq[J, \[Alpha], k];

matM = mat[[1 ;; ;; 2, 1 ;; ;; 2]];    (* indices 1,3,...,2J+1 : size J+1 *)
matm = mat[[2 ;; ;; 2, 2 ;; ;; 2]];    (* indices 2,4,...,2J   : size J   *)

Print["Diagonalizing parity sectors..."];
{eigenvalM, eigenvecM} = Eigensystem[matM, Method -> "Direct"];
{eigenvalm, eigenvecm} = Eigensystem[matm, Method -> "Direct"];

(* Embed sector eigenvectors back into the full d-dimensional basis *)
zerosM = ConstantArray[0., J];
zerosm = ConstantArray[0., J + 1];
neweigenvecM = Table[Riffle[eigenvecM[[n]], zerosM], {n, J + 1}];
neweigenvecm = Table[Riffle[zerosm, eigenvecm[[n]]], {n, J}];
neweigenvec  = Riffle[neweigenvecM, neweigenvecm];


(* ::Chapter:: *)
(*Overlaps and D2*)


Print["Projecting coherent states into Floquet basis..."];
overlaps  = Conjugate[neweigenvec] . Transpose[statesCompiled];

Print["Computing IPR..."];
iprValues = Total[Abs[overlaps]^4, {1}];

(* D2 = -log(IPR) / log(d),  d = 2J+1 *)
d2Values  = -Log[iprValues] / Log[d];


(* ================================================================ *)
(* 4. PLOTTING                                                       *)
(* ================================================================ *)


(* ::Chapter:: *)
(*D2 Density Plot + Poincare Overlay*)


d2plot = ListDensityPlot[
  Transpose[{region[[All, 1]], region[[All, 2]], d2Values}],
  PlotRange     -> All,
  PlotTheme     -> "Detailed",
  ColorFunction -> "SunsetColors",
  PlotLabel     -> Style[
    "D\:2082, J = " <> ToString[J] <>
    ", \[Alpha] = " <> ToString[N[\[Alpha]]] <>
    ", k = " <> ToString[NumberForm[k, {4, 2}]], 25, Black],
  ImageSize     -> 1000,
  FrameLabel    -> {Style["Q", 40, Black], Rotate[Style["P", 40, Black], -90 Degree]},
  LabelStyle    -> Directive[30, Black],
  AspectRatio   -> 1,
  Frame         -> True,
  InterpolationOrder -> 3];

overlay = Show[{d2plot, poincare}];

Export[
  "D2_J_" <> ToString[J] <>
  "_alpha_" <> StringReplace[ToString[N[\[Alpha]]], "." -> "p"] <>
  "_k_" <> StringReplace[ToString[k], "." -> "p"] <> ".png",
  overlay, ImageResolution -> 300]
