(* ::Package:: *)

(* QKT.wl \[LongDash] Quantum & Classical Kicked Top Efficient Operators *)

(* If ForScience paclet not installed, install it.
   See https://github.com/MMA-ForScience/ForScience *)
If[Length[PacletFind["ForScience"]] == 0,
    PacletInstall[FileNameJoin[{DirectoryName[$InputFileName], "ForScience-0.88.45.paclet"}]]
];


(* ::Section:: *)
(*Begin package*)


BeginPackage["QKT`"];


(* For nice formatting of usage messages, see https://github.com/MMA-ForScience/ForScience *)
<<ForScience`;


(* ::Section:: *)
(*Usage definitions*)


(* ::Subsection::Closed:: *)
(*Quantum Kicked Top*)


generateSpinOperators::usage = "generateSpinOperators[J] returns an association of {Sx, Sy, Sz, Sx2} spin operators for total spin J.";


Floq::usage = "Floq[J, \[Alpha], k] returns the Floquet operator with rotation along Sz and twist around Sx.";


Floqsym::usage = "Floqsym[J, \[Alpha], k] returns the Symmetric Floquet operator with rotation along Sz and twist around Sx.";


Floqn::usage = "Floqn[J, \[Alpha], k, nVec] returns the Floquet operator with twist along an arbitrary axis defined by nVec.";


Floqb::usage = "Floqb[J,\[Alpha],k,\[Gamma]] returns the Floquet operator which breaks 1 time reversal symmetries, COE.";


Floqbn::usage = "Floqbn[J,\[Alpha],k,nVec,\[Gamma]] returns the Floquet operator which can break both time reversal symmetries, CUE.";


generateCoherentStateCompiler::usage = "generateCoherentStateCompiler[] returns a compiled function that generates spin coherent states.";


MeanLevelSpacingRatio::usage = "MeanLevelSpacingRatio[eigenvalues] returns the average level spacing ratio <r>.";


(* ::Subsection::Closed:: *)
(*Spectral Statistics*)


UnfoldLinear::usage = "UnfoldLinear[phases] linearly unfolds sorted phases on [0,2Pi) to mean spacing 1.";
ValidateUnfolding::usage = "ValidateUnfolding[specData] returns staircase residuals and KS statistic.";
NNSpacings::usage = "NNSpacings[specData] returns nearest-neighbor spacings with wrap-around.";
GetQKTParitySpectra::usage = "GetQKTParitySpectra[j, alpha, k] returns <|Even->specData, Odd->specData|>.";
GenerateQKTEnsemble::usage = "GenerateQKTEnsemble[j, alpha, kList] generates QKT spectra for both parity sectors.";
RandomCOESpectrum::usage = "RandomCOESpectrum[n] returns one COE(n) unfolded spectrum.";
GenerateCOEEnsemble::usage = "GenerateCOEEnsemble[n, nReal] generates nReal COE(n) unfolded spectra.";
CompileLevelCounter::usage = "CompileLevelCounter[extSpectrum] returns a compiled binary-search counting function.";
ComputeSigma2Sliding::usage = "ComputeSigma2Sliding[counterFunc, nLevels, lValues] computes Sigma2 via sliding windows with dual estimator.";
EvaluateSpectrumSigma2::usage = "EvaluateSpectrumSigma2[specData, lValues] computes Sigma2 for one spectrum.";
EnsembleSigma2::usage = "EnsembleSigma2[spectraList, lValues] returns ensemble-averaged Sigma2 with standard errors.";
SingleSpectrumSigma2::usage = "SingleSpectrumSigma2[specData, lValues] returns Sigma2 association for a single spectrum (spectral averaging only).";
WignerCOE::usage = "WignerCOE[s] returns the COE (beta=1) Wigner surmise P(s).";
IntWignerCOE::usage = "IntWignerCOE[s] returns the integrated COE Wigner surmise.";
COEAsymptoticSigma2::usage = "COEAsymptoticSigma2[L] returns the asymptotic COE number variance.";
CompileSFF::usage = "CompileSFF is a compiled SFF function.";
EnsembleSFF::usage = "EnsembleSFF[spectraList, tauMax, nPts] returns ensemble-averaged SFF.";
COESFFPrediction::usage = "COESFFPrediction[tau] returns the COE SFF prediction.";
COESFFAnalytical::usage = "COESFFAnalytical[tau] returns the analytical COE SFF (large-N).";
SingleSFF::usage = "SingleSFF[specData, tauMax, nPts] computes SFF for a single spectrum.";
PoolSpacings::usage = "PoolSpacings[spectraList] returns sorted pooled NN spacings from an ensemble.";
PoolSpacingRatios::usage = "PoolSpacingRatios[spectraList] returns pooled spacing ratios from an ensemble.";
PrCOE::usage = "PrCOE[r] returns the COE spacing ratio distribution (Atas et al 2013).";
PrCUE::usage = "PrCUE[r] returns the CUE spacing ratio distribution.";
PrPoisson::usage = "PrPoisson[r] returns the Poisson spacing ratio distribution.";
IntPrCOE::usage = "IntPrCOE[r] returns the integrated COE spacing ratio CDF.";
KLDivergence::usage = "KLDivergence[pData, pTheory] returns the KL divergence D_KL(data||theory) from histograms.";
GetFluctuations::usage = "GetFluctuations[specData] returns the fluctuating part of the spectral staircase.";
FluctuationPowerSpectrum::usage = "FluctuationPowerSpectrum[fluc] returns the power spectrum |delta_k|^2.";
LengthSpectrum::usage = "LengthSpectrum[fluc] returns the length spectrum (FT of fluctuating density).";
GetQKTParityEigensystem::usage = "GetQKTParityEigensystem[j, alpha, k] returns eigenvalues, eigenvectors, and unfolded spectra for both sectors.";
IPR::usage = "IPR[vec] returns the inverse participation ratio of a vector.";


(* ::Subsection::Closed:: *)
(*Classical Kicked Top*)


ClassicalToSphere::usage = "ClassicalToSphere[qvec] projects a 2D disk point to a 3D sphere coordinate.";


ClassicalToDisk::usage = "ClassicalToDisk[svec] projects a 3D sphere coordinate back to the 2D disk.";


ClassicalStep::usage = "ClassicalStep[svec, alpha, k, order] performs one iteration of the Kicked Top map.";


ClassicalMapeo::usage = "ClassicalMapeo[qini, alpha, k, n, order] generates a full trajectory of n points starting from qini.";


(* ::Section:: *)
(*Beginning of Package*)


Begin["`Private`"];


(* ::Section:: *)
(*Routine definitions*)


(* ::Subsection::Closed:: *)
(*Quantum definitions*)


ClearAll[sz, sp, sm, sx, sy, sx2];


(* Spin operators: memoized for efficiency *)
sz[j_] := sz[j] = DiagonalMatrix[Table[-j + k - 1, {k, 1, 2 j + 1}]];
sp[j_] := sp[j] = DiagonalMatrix[Table[Sqrt[j (j + 1) - m (m + 1)], {m, -j, j - 1}], 1];
sm[j_] := sm[j] = DiagonalMatrix[Table[Sqrt[j (j + 1) - m (m - 1)], {m, -j + 1, j}], -1];
sx[j_] := sx[j] = SparseArray[(1/2) (sp[j] + sm[j])];
sy[j_] := sy[j] = SparseArray[(1/(2 I)) (sp[j] - sm[j])];
sx2[j_] := sx2[j] = sx[j] . sx[j];


generateSpinOperators[j_] := <| "Sx" -> sx[j], "Sy" -> sy[j], "Sz" -> sz[j], "Sx2" -> sx2[j] |>;


(* Coherent state generator *)
generateCoherentStateCompiler[] := Compile[{{J, _Integer}, {q, _Real}, {p, _Real}},
   Module[{dim, \[Alpha], result, mVals, logs, maxLog, terms, norm},
    dim = 2 J + 1;
    result = ConstantArray[0.0 + 0. I, dim];
    \[Alpha] = (q + I p)/Sqrt[4 - (q^2 + p^2)];
    If[Abs[\[Alpha]] == 0,
     result[[1]] = 1.0 + 0. I,
     mVals = Range[-J, J];
     logs = 0.5 (LogGamma[2 J + 1] - LogGamma[J + 1 + mVals] - LogGamma[J + 1 - mVals]) + (J + mVals) Log[\[Alpha]];
     maxLog = Max[Re[logs]];
     terms = Exp[logs - maxLog];
     norm = Sqrt[Total[Abs[terms]^2]];
     result = terms/norm;
     ];
    result
    ],
   CompilationTarget -> "WVM", RuntimeOptions -> "Speed", Parallelization -> True, RuntimeAttributes -> {Listable}
];


twistPart[j_, k_] := MatrixExp[(-I k/(2 j)) sx2[j]];
freePart[j_, \[Alpha]_] := MatrixExp[(-I \[Alpha]) sz[j]];


Floq[j_, \[Alpha]_, k_] := twistPart[j, k] . freePart[j, \[Alpha]];
Floqsym[j_, \[Alpha]_, k_] := freePart[j, \[Alpha]/2] . twistPart[j, k] . freePart[j, \[Alpha]/2];


twistPartGeneral[j_, \[Alpha]_, nVec_] := twistPartGeneral[j, \[Alpha], nVec] = Module[{u = Normalize[nVec], gen},
  gen = u . {sx[j], sy[j], sz[j]};
  MatrixExp[-I \[Alpha] gen]
];


Floqn[j_, \[Alpha]_, k_, nVec_] := twistPart[j, k] . twistPartGeneral[j, \[Alpha], nVec];


MeanLevelSpacingRatio[eigenvalues_] := Mean[Min /@ Transpose[{#, 1/#}] &[Ratios[Differences[Sort[eigenvalues]]]]];


(* ::Subsection:: *)
(*Spectral Statistics definitions*)


UnfoldLinear[phases_List] :=
  Module[{n = Length[phases]},
    <|"N" -> n, "Phases" -> phases, "Spectrum" -> (n / (2 Pi)) * phases|>
  ];

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

WignerCOE[s_] := (Pi / 2) s Exp[-Pi s^2 / 4];
IntWignerCOE[s_] := 1 - Exp[-Pi s^2 / 4];   (* exact closed form *)
COEAsymptoticSigma2[L_] := (2 / Pi^2) (Log[2 Pi L] + 1 + EulerGamma - Pi^2 / 8);
COESFFPrediction[tau_] := Min[2 tau, 1];

GetQKTParitySpectra[jParam_, alphaParam_, kParam_] :=
  Module[{mat, matEven, matOdd, eigsEven, eigsOdd, phEven, phOdd},
    mat = Floq[jParam, alphaParam, kParam];
    matEven = mat[[1 ;; ;; 2, 1 ;; ;; 2]];
    matOdd  = mat[[2 ;; ;; 2, 2 ;; ;; 2]];
    eigsEven = Eigenvalues[matEven, Method -> "Direct"];
    eigsOdd  = Eigenvalues[matOdd, Method -> "Direct"];
    phEven = Sort[Mod[Arg[eigsEven], 2 Pi]];
    phOdd  = Sort[Mod[Arg[eigsOdd], 2 Pi]];
    <|"Even" -> UnfoldLinear[phEven], "Odd" -> UnfoldLinear[phOdd]|>
  ];

GenerateQKTEnsemble[jParam_, alphaParam_, kList_List] :=
  Module[{allSpectra},
    DistributeDefinitions[Floq, GetQKTParitySpectra, UnfoldLinear,
      twistPart, freePart, sx, sy, sz, sp, sm, sx2];
    allSpectra = ParallelTable[
      GetQKTParitySpectra[jParam, alphaParam, k],
      {k, kList}, Method -> "CoarsestGrained"
    ];
    <|"Even" -> allSpectra[[All, "Even"]], "Odd" -> allSpectra[[All, "Odd"]]|>
  ];

RandomCOESpectrum[n_Integer] :=
  Module[{z, q, r, diag, cue, sym, eigs, phases},
    z = RandomVariate[NormalDistribution[], {n, n}] +
        I * RandomVariate[NormalDistribution[], {n, n}];
    {q, r} = QRDecomposition[z];
    diag = DiagonalMatrix[Sign[Diagonal[r]]];
    cue = diag . q;
    sym = cue . Transpose[cue];
    eigs = Eigenvalues[sym];
    phases = Sort[Mod[Arg[eigs], 2 Pi]];
    UnfoldLinear[phases]
  ];

GenerateCOEEnsemble[n_Integer, nReal_Integer] :=
  Module[{},
    DistributeDefinitions[RandomCOESpectrum, UnfoldLinear];
    ParallelTable[RandomCOESpectrum[n], {nReal}, Method -> "CoarsestGrained"]
  ];

CompileLevelCounter[extendedSpectrum_List] :=
  Module[{packed, m},
    packed = Developer`ToPackedArray[N[extendedSpectrum]];
    m = Length[packed];
    Compile[{{z, _Real}},
      Module[{lo = 1, hi = m, mid = 0},
        While[lo <= hi,
          mid = Quotient[lo + hi, 2];
          If[packed[[mid]] <= z, lo = mid + 1, hi = mid - 1]
        ]; hi],
      CompilationTarget -> "WVM", RuntimeOptions -> "Speed",
      RuntimeAttributes -> {Listable}, Parallelization -> False
    ]
  ];

(* No hard L cap \[LongDash] user controls range via lValues *)
ComputeSigma2Sliding[counterFunc_, nLevels_Integer, lValues_List] :=
  Module[{centers, nC},
    centers = N[Range[0, nLevels - 1]];
    nC = Length[centers];
    Map[
      Function[{L},
        Module[{counts, s2theor, s2sample},
          counts = counterFunc[centers + L] - counterFunc[centers];
          s2theor  = Total[(counts - L)^2] / nC;
          s2sample = Variance[counts];
          {L, s2theor, s2sample}
        ]
      ],
      N[lValues]
    ]
  ];

EvaluateSpectrumSigma2[specData_Association, lValues_List] :=
  Module[{ext, counter},
    ext = Join[
      specData["Spectrum"] - specData["N"],
      specData["Spectrum"],
      specData["Spectrum"] + specData["N"]
    ];
    counter = CompileLevelCounter[ext];
    ComputeSigma2Sliding[counter, specData["N"], lValues]
  ];

EnsembleSigma2[spectraList_List, lValues_List] :=
  Module[{nMat, allResults, effectiveL, theorVals, sampleVals,
          meanTheor, seTheor, meanSample, seSample},
    nMat = Length[spectraList];
    DistributeDefinitions[CompileLevelCounter, ComputeSigma2Sliding,
      EvaluateSpectrumSigma2];
    allResults = ParallelTable[
      EvaluateSpectrumSigma2[spectraList[[i]], lValues],
      {i, nMat}, Method -> "CoarsestGrained"
    ];
    effectiveL  = allResults[[1, All, 1]];
    theorVals   = Map[#[[All, 2]] &, allResults];
    sampleVals  = Map[#[[All, 3]] &, allResults];
    meanTheor   = Mean[theorVals];
    seTheor     = StandardDeviation[theorVals] / Sqrt[nMat * 1.0];
    meanSample  = Mean[sampleVals];
    seSample    = StandardDeviation[sampleVals] / Sqrt[nMat * 1.0];
    <|"LValues" -> effectiveL,
      "Sigma2TheorMean" -> meanTheor, "SETheorMean" -> seTheor,
      "Sigma2SampleMean" -> meanSample, "SESampleMean" -> seSample|>
  ];

(* Single spectrum: spectral averaging only, no ensemble *)
SingleSpectrumSigma2[specData_Association, lValues_List] :=
  Module[{results},
    results = EvaluateSpectrumSigma2[specData, lValues];
    <|"LValues" -> results[[All, 1]],
      "Sigma2Theor" -> results[[All, 2]],
      "Sigma2Sample" -> results[[All, 3]]|>
  ];

CompileSFF =
  Compile[{{spec, _Real, 1}, {taus, _Real, 1}},
    Module[{n = Length[spec], cs, sn},
      Table[cs = Total[Cos[2.0 Pi spec t]];
        sn = Total[Sin[2.0 Pi spec t]];
        (cs^2 + sn^2) / n, {t, taus}]
    ], CompilationTarget -> "WVM", RuntimeOptions -> "Speed"
  ];

EnsembleSFF[spectraList_List, tauMax_Real, nPts_Integer] :=
  Module[{nMat, taus, specArray, allK, meanK, seK},
    nMat = Length[spectraList];
    taus = N[Range[1, nPts] (tauMax / nPts)];
    specArray = Map[#["Spectrum"] &, spectraList];
    DistributeDefinitions[CompileSFF, taus, specArray];
    allK = ParallelTable[
      CompileSFF[specArray[[i]], taus],
      {i, nMat}, Method -> "CoarsestGrained"
    ];
    meanK = Mean[allK];
    seK = StandardDeviation[allK] / Sqrt[nMat * 1.0];
    <|"Tau" -> taus, "MeanK" -> meanK, "SE" -> seK|>
  ];

COESFFAnalytical[tau_] :=
  Piecewise[{
    {2 tau - tau Log[1 + 2 tau], 0 < tau < 1},
    {2 - tau Log[(2 tau + 1) / (2 tau - 1)], tau >= 1}
  }];

SingleSFF[specData_Association, tauMax_Real, nPts_Integer] :=
  Module[{taus},
    taus = N[Range[1, nPts] (tauMax / nPts)];
    <|"Tau" -> taus, "K" -> CompileSFF[specData["Spectrum"], taus]|>
  ];

(* --- Pooled statistics (fast, no per-element Map) --- *)
PoolSpacings[spectraList_List] :=
  Sort[Flatten[Map[NNSpacings, spectraList]]];

PoolSpacingRatios[spectraList_List] :=
  Flatten[Map[
    Function[{specData},
      Module[{sp},
        sp = Differences[specData["Spectrum"]];
        Min /@ Transpose[{Most[sp] / Rest[sp], Rest[sp] / Most[sp]}]
      ]
    ],
    spectraList
  ]];

(* --- Spacing ratio distributions (Atas et al 2013) --- *)
PrCOE[r_]     := (27 / 4) (r + r^2) / (1 + r + r^2)^(5/2);
PrCUE[r_]     := (81 Sqrt[3] / (4 Pi)) (r + r^2)^2 / (1 + r + r^2)^4;
PrPoisson[r_] := 2 / (1 + r)^2;

IntPrCOE[rMax_?NumericQ] := NIntegrate[PrCOE[r], {r, 0, rMax}];

(* --- KL divergence from binned histograms --- *)
(* pData and pTheory are probability densities at same bin centers *)
KLDivergence[dataVals_List, theoryFunc_, binSpec_: {0, 4, 0.05}] :=
  Module[{edges, centers, hData, pD, pT, mask, dx},
    edges = Range[Sequence @@ binSpec];
    dx = binSpec[[3]];
    centers = MovingAverage[edges, 2];
    hData = HistogramList[dataVals, edges, "ProbabilityDensity"][[2]];
    pD = hData * dx;  (* bin probabilities *)
    pT = Map[theoryFunc, centers] * dx;
    mask = Flatten[Position[pD, _?(# > 0 &)]];
    Total[(pD[[mask]]) * Log[pD[[mask]] / pT[[mask]]]]
  ];


GetFluctuations[specData_Association] :=
  Module[{x = specData["Spectrum"], n = specData["N"], raw},
    raw = Range[n] - x;
    <|"x" -> x, "delta" -> raw - Mean[raw], "N" -> n|>
  ];

FluctuationPowerSpectrum[fluc_Association] :=
  Module[{ft, n, nHalf, freqs},
    n = fluc["N"]; nHalf = Floor[n / 2];
    ft = Abs[Fourier[fluc["delta"]]]^2;
    freqs = Range[0, n - 1] / n;
    <|"Freq" -> freqs[[1 ;; nHalf]], "Power" -> ft[[1 ;; nHalf]]|>
  ];

LengthSpectrum[fluc_Association] :=
  Module[{n, nHalf, ft, periods},
    n = fluc["N"]; nHalf = Floor[n / 2];
    ft = Fourier[fluc["delta"], FourierParameters -> {1, -1}];
    periods = Range[0, nHalf - 1];
    <|"Period" -> periods, "Amplitude" -> Re[ft[[1 ;; nHalf]]]|>
  ];

(* --- Eigensystem with both parity sectors --- *)
GetQKTParityEigensystem[jParam_, alphaParam_, kParam_] :=
  Module[{mat, matEven, matOdd, esEven, esOdd, phEven, phOdd,
          vecsEven, vecsOdd, idxE, idxO},
    mat = Floq[jParam, alphaParam, kParam];
    matEven = mat[[1 ;; ;; 2, 1 ;; ;; 2]];
    matOdd  = mat[[2 ;; ;; 2, 2 ;; ;; 2]];
    esEven = Eigensystem[matEven, Method -> "Direct"];
    esOdd  = Eigensystem[matOdd, Method -> "Direct"];
    phEven = Mod[Arg[esEven[[1]]], 2 Pi];
    phOdd  = Mod[Arg[esOdd[[1]]], 2 Pi];
    idxE = Ordering[phEven]; idxO = Ordering[phOdd];
    vecsEven = esEven[[2, idxE]];
    vecsOdd  = esOdd[[2, idxO]];
    <|"Even" -> <|"Spec" -> UnfoldLinear[phEven[[idxE]]],
                  "Vectors" -> vecsEven|>,
      "Odd"  -> <|"Spec" -> UnfoldLinear[phOdd[[idxO]]],
                  "Vectors" -> vecsOdd|>|>
  ];

IPR[vec_] := Total[Abs[vec]^4] / Total[Abs[vec]^2]^2;


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
     (* Rotation along Z then Kick along X *)
     r1 = ca*sx - sa*sy; r2 = sa*sx + ca*sy; r3 = sz;
     sx = r1; sy = r2; sz = r3;
     sxc = Cos[k*sx]; sxs = Sin[k*sx];
     r1 = sx; r2 = sxc*sy - sxs*sz; r3 = sxs*sy + sxc*sz;
     ,
     (* Kick along X then Rotation along Z *)
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
