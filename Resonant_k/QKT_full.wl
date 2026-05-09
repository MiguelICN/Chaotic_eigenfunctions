(* ::Package:: *)

(* QKT_full.wl \[LongDash] Quantum & Classical Kicked Top                     *)
(* Operators, Parity Decomposition, Spectral Statistics,             *)
(* Coherent States, and Phase-Space Tools                            *)
If[Length[PacletFind["ForScience"]] == 0,
    PacletInstall[FileNameJoin[{DirectoryName[$InputFileName], "ForScience-0.88.45.paclet"}]]
];


(* ::Section::Closed:: *)
(*Begin package*)


BeginPackage["QKT`"];


<<ForScience`;


(* ::Section::Closed:: *)
(*Usage definitions*)


(* ::Subsection::Closed:: *)
(*Quantum Kicked Top*)


generateSpinOperators::usage = "generateSpinOperators[j] returns <|\"Sx\"->..., \"Sy\"->..., \"Sz\"->..., \"Sx2\"->...|> for total spin j. Descending-m basis.";


Floq::usage = "Floq[j, \[Alpha], k] returns the Floquet operator exp(-ik Sx^2/(2j)) exp(-i\[Alpha] Sz).";


Floqsym::usage = "Floqsym[j, \[Alpha], k] returns the symmetrized Floquet operator exp(-i\[Alpha] Sz/2) exp(-ik Sx^2/(2j)) exp(-i\[Alpha] Sz/2).";


Floqn::usage = "Floqn[j, \[Alpha], k, nVec] returns the Floquet operator with twist along Sx and rotation along nVec.";


Floqb::usage = "Floqb[j, \[Alpha], k, \[Gamma]] returns the Floquet operator breaking one time-reversal symmetry.";


Floqbn::usage = "Floqbn[j, \[Alpha], k, nVec, \[Gamma]] returns the Floquet operator breaking both time-reversal symmetries.";


generateCoherentStateCompiler::usage = "generateCoherentStateCompiler[] returns a compiled function f[J, q, p] for spin-j coherent states on the stereographic disk.";


MeanLevelSpacingRatio::usage = "MeanLevelSpacingRatio[eigenvalues] returns <r>, the mean consecutive spacing ratio.";


(* ::Subsection::Closed:: *)
(*Parity Decomposition*)


ParitySectorIndices::usage = "ParitySectorIndices[j] returns <|\"Even\"->idxList, \"Odd\"->idxList|> for the two decoupled parity sectors. Works for integer and half-integer j.";


GetQKTParitySpectra::usage = "GetQKTParitySpectra[j, \[Alpha], k] returns <|\"Even\"->specData, \"Odd\"->specData|> with unfolded spectra for both parity sectors.";


GetQKTParityEigensystem::usage = "GetQKTParityEigensystem[j, \[Alpha], k] returns eigenvalues, eigenvectors, and full-basis reconstructed eigenvectors for both parity sectors.";


GenerateQKTEnsemble::usage = "GenerateQKTEnsemble[j, \[Alpha], kList] generates QKT spectra over kList (parallelized). Returns <|\"Even\"->list, \"Odd\"->list|>.";


(* ::Subsection::Closed:: *)
(*Spectral Statistics*)


UnfoldLinear::usage = "UnfoldLinear[phases] linearly unfolds sorted phases on [0,2\[Pi]) to mean spacing 1.";
ValidateUnfolding::usage = "ValidateUnfolding[specData] returns staircase residuals and KS statistic.";
NNSpacings::usage = "NNSpacings[specData] returns nearest-neighbor spacings with wrap-around.";
PoolSpacings::usage = "PoolSpacings[spectraList] pools and sorts all NN spacings from an ensemble.";
PoolSpacingRatios::usage = "PoolSpacingRatios[spectraList] pools all spacing ratios from an ensemble.";
GetFluctuations::usage = "GetFluctuations[specData] returns the fluctuating part of the staircase, centered on zero.";
LengthSpectrum::usage = "LengthSpectrum[fluc] returns the DFT of the fluctuating density.";
FluctuationPowerSpectrum::usage = "FluctuationPowerSpectrum[fluc] returns the power spectrum of staircase fluctuations.";
IPR::usage = "IPR[vec] returns the inverse participation ratio.";


(* ::Subsection::Closed:: *)
(*Analytical RMT Predictions*)


WignerCOE::usage = "WignerCOE[s] returns the COE Wigner surmise P(s).";
IntWignerCOE::usage = "IntWignerCOE[s] returns the integrated COE Wigner surmise I(s).";
PrCOE::usage = "PrCOE[r] returns the COE spacing ratio distribution (Atas et al. 2013).";
PrCUE::usage = "PrCUE[r] returns the CUE spacing ratio distribution.";
PrPoisson::usage = "PrPoisson[r] returns the Poisson spacing ratio distribution.";
COEAsymptoticSigma2::usage = "COEAsymptoticSigma2[L] returns the asymptotic COE number variance.";
COESFFAnalytical::usage = "COESFFAnalytical[\[Tau]] returns the analytical COE SFF.";


(* ::Subsection::Closed:: *)
(*Number Variance and SFF*)


CompileLevelCounter::usage = "CompileLevelCounter[extSpectrum] returns a compiled binary-search counting function.";
ComputeSigma2Sliding::usage = "ComputeSigma2Sliding[counterFunc, nLevels, lValues] computes Sigma2 via sliding windows.";
EvaluateSpectrumSigma2::usage = "EvaluateSpectrumSigma2[specData, lValues] computes Sigma2 for one spectrum.";
EnsembleSigma2::usage = "EnsembleSigma2[spectraList, lValues] returns ensemble-averaged Sigma2 with errors.";
SingleSpectrumSigma2::usage = "SingleSpectrumSigma2[specData, lValues] returns Sigma2 from spectral averaging only.";
CompileSFF::usage = "CompileSFF[spectrum, tauValues] returns K(\[Tau]).";
EnsembleSFF::usage = "EnsembleSFF[spectraList, tauMax, nPts] returns ensemble-averaged SFF with errors.";
SingleSFF::usage = "SingleSFF[specData, tauMax, nPts] computes K(\[Tau]) for a single spectrum.";


(* ::Subsection::Closed:: *)
(*COE Random Matrices*)


RandomCOESpectrum::usage = "RandomCOESpectrum[n] generates one COE(n) unfolded spectrum.";
GenerateCOEEnsemble::usage = "GenerateCOEEnsemble[n, nReal] generates nReal COE(n) unfolded spectra (parallelized).";


(* ::Subsection::Closed:: *)
(*Classical Kicked Top*)


ClassicalToSphere::usage = "ClassicalToSphere[{q1,q2}] maps a stereographic disk point to S^2.";
ClassicalToDisk::usage = "ClassicalToDisk[{Sx,Sy,Sz}] maps S^2 back to the disk.";
ClassicalStep::usage = "ClassicalStep[{Sx,Sy,Sz}, \[Alpha], k, order] performs one stroboscopic step.";
ClassicalMapeo::usage = "ClassicalMapeo[{q1,q2}, \[Alpha], k, n, order] generates a trajectory. Compiled, listable.";


(* ::Section::Closed:: *)
(*Beginning of Package*)


Begin["`Private`"];


(* ::Section:: *)
(*Routine definitions*)


(* ::Subsection::Closed:: *)
(*Quantum definitions*)


ClearAll[sz, sp, sm, sx, sy, sx2];


(* 1. Sz in the standard descending basis: {j, j-1, ..., -j} *)
sz[j_] := sz[j] = SparseArray[Band[{1, 1}] -> Range[j, -j, -1]];

(* 2. S+ is on the superdiagonal. Iterating ket m from j-1 down to -j *)
sp[j_] := sp[j] = SparseArray[Band[{1, 2}] -> 
    Table[Sqrt[j (j + 1) - m (m + 1)] // N, {m, j - 1, -j, -1}], 
    {2 j + 1, 2 j + 1}];

(* 3. S- is on the subdiagonal. Iterating ket m from j down to -j+1 *)
sm[j_] := sm[j] = SparseArray[Band[{2, 1}] -> 
    Table[Sqrt[j (j + 1) - m (m - 1)] // N, {m, j, -j + 1, -1}], 
    {2 j + 1, 2 j + 1}];

(* 4. Sx and Sy (Already sparse since sp and sm are sparse) *)
sx[j_] := sx[j] = (1/2) (sp[j] + sm[j]);
sy[j_] := sy[j] = (1/(2 I)) (sp[j] - sm[j]);

(* 5. Matrix multiplication of sparse arrays *)
sx2[j_] := sx2[j] = sx[j] . sx[j];


generateSpinOperators[j_] := <|"Sx" -> sx[j], "Sy" -> sy[j], "Sz" -> sz[j], "Sx2" -> sx2[j]|>;


(* ------------------------------------------------------------------ *)
(* Coherent State Compiler (Descending Basis, Boundary-Protected)     *)
(* ------------------------------------------------------------------ *)
generateCoherentStateCompiler[] := Compile[{{J, _Integer}, {q, _Real}, {p, _Real}},
   Module[{dim, \[Alpha], result, mVals, logs, maxLog, terms, norm, r2},
    dim = 2 J + 1;
    result = ConstantArray[0.0 + 0. I, dim];
    r2 = q^2 + p^2;
    If[r2 >= 4.0,
     result[[dim]] = 1.0 + 0. I,
     \[Alpha] = (q + I p)/Sqrt[4 - r2];
     If[Abs[\[Alpha]] == 0.0,
      result[[1]] = 1.0 + 0. I,
      mVals = Range[J, -J, -1];
      logs = 0.5 (LogGamma[2 J + 1] - LogGamma[J + 1 + mVals] - LogGamma[J + 1 - mVals]) + (J - mVals) Log[\[Alpha]];
      maxLog = Max[Re[logs]];
      terms = Exp[logs - maxLog];
      norm = Sqrt[Total[Abs[terms]^2]];
      result = terms/norm;
      ]
     ];
    result
    ],
   CompilationTarget -> "WVM", RuntimeOptions -> "Speed", Parallelization -> True, RuntimeAttributes -> {Listable}
];


twistPart[j_, k_] := MatrixExp[(-I k/(2 j)) sx2[j]];
freePart[j_, \[Alpha]_] := SparseArray[Band[{1, 1}] -> Exp[-I \[Alpha] Range[j, -j, -1]]];


Floq[j_, \[Alpha]_, k_] := twistPart[j, k] . freePart[j, \[Alpha]];
Floqsym[j_, \[Alpha]_, k_] := freePart[j, \[Alpha]/2] . twistPart[j, k] . freePart[j, \[Alpha]/2];


twistPartGeneral[j_, \[Alpha]_, nVec_] := twistPartGeneral[j, \[Alpha], nVec] = Module[{u = Normalize[nVec], gen},
  gen = u . {sx[j], sy[j], sz[j]};
  MatrixExp[-I \[Alpha] gen]
];


Floqn[j_, \[Alpha]_, k_, nVec_] := twistPart[j, k] . twistPartGeneral[j, \[Alpha], nVec];


MeanLevelSpacingRatio[eigenvalues_] := Mean[Min /@ Transpose[{#, 1/#}] &[Ratios[Differences[Sort[eigenvalues]]]]];


(* ::Subsection:: *)
(*Parity Decomposition*)


(* Sx^2 couples m to m +/- 2, so alternating array indices always decouple.
   Odd integer j: index 1 has m=j (odd), so Even-m sector starts at index 2.
   Even integer j or half-integer j: index 1 starts the first sector. *)

ParitySectorIndices[j_] :=
  Module[{d = 2 j + 1},
    If[IntegerQ[j] && OddQ[j],
      <|"Even" -> Range[2, d, 2], "Odd" -> Range[1, d, 2]|>,
      <|"Even" -> Range[1, d, 2], "Odd" -> Range[2, d, 2]|>
    ]
  ];


UnfoldLinear[phases_List] :=
  Module[{n = Length[phases]},
    <|"N" -> n, "Phases" -> phases, "Spectrum" -> (n / (2 Pi)) * phases|>
  ];


GetQKTParitySpectra[jParam_, alphaParam_, kParam_] :=
  Module[{mat, d, idx, matEven, matOdd, eigsEven, eigsOdd, phEven, phOdd},
    mat = Floq[jParam, alphaParam, kParam];
    d = 2 jParam + 1;
    idx = ParitySectorIndices[jParam];
    matEven = mat[[idx["Even"], idx["Even"]]];
    matOdd  = mat[[idx["Odd"], idx["Odd"]]];
    eigsEven = Eigenvalues[matEven, Method -> "Direct"];
    eigsOdd  = Eigenvalues[matOdd, Method -> "Direct"];
    phEven = Sort[Mod[Arg[eigsEven], 2 Pi]];
    phOdd  = Sort[Mod[Arg[eigsOdd], 2 Pi]];
    <|"Even" -> UnfoldLinear[phEven], "Odd" -> UnfoldLinear[phOdd]|>
  ];


GetQKTParityEigensystem[jParam_, alphaParam_, kParam_] :=
  Module[{mat, d, idx, matEven, matOdd, esEven, esOdd,
          phEven, phOdd, vecsEven, vecsOdd, idxE, idxO,
          fullVecsEven, fullVecsOdd},
    mat = Floq[jParam, alphaParam, kParam];
    d = 2 jParam + 1;
    idx = ParitySectorIndices[jParam];
    matEven = mat[[idx["Even"], idx["Even"]]];
    matOdd  = mat[[idx["Odd"], idx["Odd"]]];
    esEven = Eigensystem[matEven, Method -> "Direct"];
    esOdd  = Eigensystem[matOdd, Method -> "Direct"];
    phEven = Mod[Arg[esEven[[1]]], 2 Pi];
    phOdd  = Mod[Arg[esOdd[[1]]], 2 Pi];
    idxE = Ordering[phEven]; idxO = Ordering[phOdd];
    vecsEven = esEven[[2, idxE]];
    vecsOdd  = esOdd[[2, idxO]];
    fullVecsEven = Table[
      Module[{vec = ConstantArray[0. + 0. I, d]},
        vec[[idx["Even"]]] = vecsEven[[i]]; vec],
      {i, Length[idx["Even"]]}];
    fullVecsOdd = Table[
      Module[{vec = ConstantArray[0. + 0. I, d]},
        vec[[idx["Odd"]]] = vecsOdd[[i]]; vec],
      {i, Length[idx["Odd"]]}];
    <|"Even" -> <|"Quasienergies" -> esEven[[1]],
                  "Phases" -> phEven,
                  "Vectors" -> vecsEven,
                  "FullVectors" -> fullVecsEven|>,
      "Odd"  -> <|"Quasienergies" -> esOdd[[1]], 
                  "Phases" -> phOdd,
                  "Vectors" -> vecsOdd,
                  "FullVectors" -> fullVecsOdd|>|>
  ];


GenerateQKTEnsemble[jParam_, alphaParam_, kList_List] :=
  Module[{allSpectra},
    DistributeDefinitions[Floq, GetQKTParitySpectra, UnfoldLinear,
      ParitySectorIndices, twistPart, freePart, sx, sy, sz, sp, sm, sx2];
    allSpectra = ParallelTable[
      GetQKTParitySpectra[jParam, alphaParam, k],
      {k, kList}, Method -> "CoarsestGrained"
    ];
    <|"Even" -> allSpectra[[All, "Even"]],
      "Odd"  -> allSpectra[[All, "Odd"]]|>
  ];


(* ::Subsection::Closed:: *)
(*Spectral Statistics*)


ValidateUnfolding[specData_Association] :=
  Module[{phases, n, empiricalCDF, linearCDF, residuals, ks},
    phases = specData["Phases"]; n = specData["N"];
    empiricalCDF = Range[n] / n;
    linearCDF = phases / (2 Pi);
    residuals = n (empiricalCDF - linearCDF);
    ks = Max[Abs[residuals]] / Sqrt[n];
    <|"Phases" -> phases, "EmpiricalCDF" -> empiricalCDF,
      "LinearCDF" -> linearCDF, "Residuals" -> residuals,
      "KS" -> ks, "N" -> n|>
  ];


NNSpacings[specData_Association] :=
  Module[{s = specData["Spectrum"], n = specData["N"], diffs},
    diffs = Differences[s];
    Append[diffs, n - Last[s] + First[s]]
  ];


PoolSpacings[spectraList_List] := Sort[Flatten[Map[NNSpacings, spectraList]]];


PoolSpacingRatios[spectraList_List] :=
  Flatten[Map[
    Function[{specData},
      Module[{sp},
        sp = Differences[specData["Spectrum"]];
        Min /@ Transpose[{Most[sp] / Rest[sp], Rest[sp] / Most[sp]}]
      ]
    ], spectraList]];


GetFluctuations[specData_Association] :=
  Module[{x = specData["Spectrum"], n = specData["N"], raw},
    raw = Range[n] - x;
    <|"x" -> x, "delta" -> raw - Mean[raw], "N" -> n|>
  ];


LengthSpectrum[fluc_Association] :=
  Module[{n, nHalf, ft},
    n = fluc["N"]; nHalf = Floor[n / 2];
    ft = Fourier[fluc["delta"], FourierParameters -> {1, -1}];
    <|"Period" -> Range[0, nHalf - 1], "Amplitude" -> Re[ft[[1 ;; nHalf]]]|>
  ];


FluctuationPowerSpectrum[fluc_Association] :=
  Module[{n, nHalf, ft},
    n = fluc["N"]; nHalf = Floor[n / 2];
    ft = Abs[Fourier[fluc["delta"]]]^2;
    <|"Freq" -> Range[0, n - 1][[1 ;; nHalf]] / n, "Power" -> ft[[1 ;; nHalf]]|>
  ];


IPR[vec_] := Total[Abs[vec]^4] / Total[Abs[vec]^2]^2;


(* ::Subsection::Closed:: *)
(*Analytical RMT*)


WignerCOE[s_] := (Pi / 2) s Exp[-Pi s^2 / 4];
IntWignerCOE[s_] := 1 - Exp[-Pi s^2 / 4];
PrCOE[r_]     := (27 / 4) (r + r^2) / (1 + r + r^2)^(5/2);
PrCUE[r_]     := (81 Sqrt[3] / (4 Pi)) (r + r^2)^2 / (1 + r + r^2)^4;
PrPoisson[r_] := 2 / (1 + r)^2;
COEAsymptoticSigma2[L_] := (2 / Pi^2) (Log[2 Pi L] + 1 + EulerGamma - Pi^2 / 8);
COESFFAnalytical[tau_] := Piecewise[{
    {2 tau - tau Log[1 + 2 tau], 0 < tau < 1},
    {2 - tau Log[(2 tau + 1) / (2 tau - 1)], tau >= 1}}];


(* ::Subsection::Closed:: *)
(*Number Variance*)


CompileLevelCounter[extendedSpectrum_List] :=
  Module[{packed, m},
    packed = Developer`ToPackedArray[N[extendedSpectrum]]; m = Length[packed];
    Compile[{{z, _Real}},
      Module[{lo = 1, hi = m, mid = 0},
        While[lo <= hi, mid = Quotient[lo + hi, 2];
          If[packed[[mid]] <= z, lo = mid + 1, hi = mid - 1]]; hi],
      CompilationTarget -> "WVM", RuntimeOptions -> "Speed",
      RuntimeAttributes -> {Listable}, Parallelization -> False]];


ComputeSigma2Sliding[counterFunc_, nLevels_Integer, lValues_List] :=
  Module[{centers = N[Range[0, nLevels - 1]], nC},
    nC = Length[centers];
    Map[Function[{L}, Module[{counts, s2t, s2s},
        counts = counterFunc[centers + L] - counterFunc[centers];
        s2t = Total[(counts - L)^2] / nC; s2s = Variance[counts];
        {L, s2t, s2s}]], N[lValues]]];


EvaluateSpectrumSigma2[specData_Association, lValues_List] :=
  Module[{ext, counter},
    ext = Join[specData["Spectrum"] - specData["N"], specData["Spectrum"],
               specData["Spectrum"] + specData["N"]];
    counter = CompileLevelCounter[ext];
    ComputeSigma2Sliding[counter, specData["N"], lValues]];


EnsembleSigma2[spectraList_List, lValues_List] :=
  Module[{nMat, allR, eL, tv, sv},
    nMat = Length[spectraList];
    DistributeDefinitions[CompileLevelCounter, ComputeSigma2Sliding, EvaluateSpectrumSigma2];
    allR = ParallelTable[EvaluateSpectrumSigma2[spectraList[[i]], lValues],
      {i, nMat}, Method -> "CoarsestGrained"];
    eL = allR[[1, All, 1]]; tv = Map[#[[All, 2]] &, allR]; sv = Map[#[[All, 3]] &, allR];
    <|"LValues" -> eL,
      "Sigma2TheorMean" -> Mean[tv], "SETheorMean" -> StandardDeviation[tv] / Sqrt[nMat * 1.],
      "Sigma2SampleMean" -> Mean[sv], "SESampleMean" -> StandardDeviation[sv] / Sqrt[nMat * 1.]|>];


SingleSpectrumSigma2[specData_Association, lValues_List] :=
  Module[{r = EvaluateSpectrumSigma2[specData, lValues]},
    <|"LValues" -> r[[All, 1]], "Sigma2Theor" -> r[[All, 2]], "Sigma2Sample" -> r[[All, 3]]|>];


(* ::Subsection::Closed:: *)
(*Spectral Form Factor*)


CompileSFF = Compile[{{spec, _Real, 1}, {taus, _Real, 1}},
    Module[{n = Length[spec], cs, sn},
      Table[cs = Total[Cos[2.0 Pi spec t]]; sn = Total[Sin[2.0 Pi spec t]];
        (cs^2 + sn^2) / n, {t, taus}]],
    CompilationTarget -> "WVM", RuntimeOptions -> "Speed"];


EnsembleSFF[spectraList_List, tauMax_Real, nPts_Integer] :=
  Module[{nMat, taus, specArray, allK},
    nMat = Length[spectraList];
    taus = N[Range[1, nPts] (tauMax / nPts)];
    specArray = Map[#["Spectrum"] &, spectraList];
    DistributeDefinitions[CompileSFF, taus, specArray];
    allK = ParallelTable[CompileSFF[specArray[[i]], taus],
      {i, nMat}, Method -> "CoarsestGrained"];
    <|"Tau" -> taus, "MeanK" -> Mean[allK],
      "SE" -> StandardDeviation[allK] / Sqrt[nMat * 1.]|>];


SingleSFF[specData_Association, tauMax_Real, nPts_Integer] :=
  Module[{taus = N[Range[1, nPts] (tauMax / nPts)]},
    <|"Tau" -> taus, "K" -> CompileSFF[specData["Spectrum"], taus]|>];


(* ::Subsection::Closed:: *)
(*COE Random Matrices*)


RandomCOESpectrum[n_Integer] :=
  Module[{z, q, r, diag, cue, sym, eigs, phases},
    z = RandomVariate[NormalDistribution[], {n, n}] + I * RandomVariate[NormalDistribution[], {n, n}];
    {q, r} = QRDecomposition[z]; diag = DiagonalMatrix[Sign[Diagonal[r]]];
    cue = diag . q; sym = cue . Transpose[cue]; eigs = Eigenvalues[sym];
    phases = Sort[Mod[Arg[eigs], 2 Pi]]; UnfoldLinear[phases]];


GenerateCOEEnsemble[n_Integer, nReal_Integer] :=
  Module[{}, DistributeDefinitions[RandomCOESpectrum, UnfoldLinear];
    ParallelTable[RandomCOESpectrum[n], {nReal}, Method -> "CoarsestGrained"]];


(* ::Subsection::Closed:: *)
(*Classical definitions*)


ClassicalToSphere = Compile[{{q, _Real, 1}},
   Module[{b2 = q[[1]]^2 + q[[2]]^2},
    {q[[1]]*Sqrt[1 - 0.25*b2], q[[2]]*Sqrt[1 - 0.25*b2], 0.5*b2 - 1.0}
   ], CompilationTarget -> "WVM"];


ClassicalToDisk = Compile[{{s, _Real, 1}},
   Module[{b = Sqrt[2.0/(1.0 - s[[3]])]},
    {s[[1]]*b, s[[2]]*b}
   ], CompilationTarget -> "WVM"];


ClassicalStep = Compile[{{s, _Real, 1}, {alpha, _Real}, {k, _Real}, {order, _Integer}},
   Module[{sx = s[[1]], sy = s[[2]], sz = s[[3]], ca = Cos[alpha], sa = Sin[alpha], sxc, sxs, r1, r2, r3},
    If[order == 0,
     r1 = ca*sx - sa*sy; r2 = sa*sx + ca*sy; r3 = sz;
     sx = r1; sy = r2; sz = r3;
     sxc = Cos[k*sx]; sxs = Sin[k*sx];
     r1 = sx; r2 = sxc*sy - sxs*sz; r3 = sxs*sy + sxc*sz;
     ,
     sxc = Cos[k*sx]; sxs = Sin[k*sx];
     r1 = sx; r2 = sxc*sy - sxs*sz; r3 = sxs*sy + sxc*sz;
     sx = r1; sy = r2; sz = r3;
     r1 = ca*sx - sa*sy; r2 = sa*sx + ca*sy; r3 = sz;
    ];
    {r1, r2, r3}
   ], CompilationTarget -> "WVM"];


ClassicalMapeo = Compile[{{qini, _Real, 1}, {alpha, _Real}, {k, _Real}, {n, _Integer}, {order, _Integer}},
   Module[{s, traj},
    s = ClassicalToSphere[qini];
    traj = Table[
      s = ClassicalStep[s, alpha, k, order];
      ClassicalToDisk[s]
     , {n}];
    traj
   ],
   CompilationTarget -> "WVM",
   CompilationOptions -> {"InlineCompiledFunctions" -> True},
   RuntimeAttributes -> {Listable},
   RuntimeOptions -> "Speed"
];


(* ::Section:: *)
(*End of Package*)


End[];


EndPackage[];
