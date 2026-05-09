(* ================================================================ *)
(* FLOQUET NUMBER VARIANCE: COMPLETE DIAGNOSTIC WORKFLOW            *)
(* Systematic sliding windows + ensemble averaging                  *)
(* ================================================================ *)

(* ================================================================ *)
(* SECTION 0: DEFINITIONS                                           *)
(* ================================================================ *)

ClearAll[
  GetRawPhases, UnfoldLinear, ValidateUnfolding,
  CompileLevelCounter, ComputeSigma2Sliding,
  EvaluateSpectrumSigma2, EnsembleSigma2,
  CUESigma2Exact, PoissonSigma2,
  CompileSFF, EnsembleSFF,
  NNSpacings
];

(* --- 0a. Raw phase extraction from Floquet operator --- *)
GetRawPhases[jParam_, alphaParam_, kParam_] :=
  Module[{mat, matM, eigs},
    mat = Floq[jParam, alphaParam, kParam];
    matM = mat[[2 ;; ;; 2, 2 ;; ;; 2]];  (* parity sector *)
    eigs = Eigenvalues[matM, Method -> "Direct"];
    Sort[Mod[Arg[eigs], 2 Pi]]
  ];

(* --- 0b. Linear unfolding --- *)
UnfoldLinear[phases_List] :=
  Module[{n = Length[phases]},
    <|"N" -> n,
      "Phases" -> phases,
      "Spectrum" -> (n / (2 Pi)) * phases|>
  ];

(* --- 0c. Unfolding validation diagnostics --- *)
(* Compares empirical staircase against linear prediction *)
ValidateUnfolding[specData_Association] :=
  Module[{phases, n, empiricalCDF, linearCDF, residuals, ks},
    phases = specData["Phases"];
    n = specData["N"];
    empiricalCDF = Range[n] / n;           (* i/N *)
    linearCDF = phases / (2 Pi);           (* theta_i / 2pi *)
    residuals = n * (empiricalCDF - linearCDF); (* in units of mean spacing *)
    ks = Max[Abs[residuals]] / Sqrt[n];    (* KS-like statistic *)
    <|"Phases" -> phases,
      "EmpiricalCDF" -> empiricalCDF,
      "LinearCDF" -> linearCDF,
      "Residuals" -> residuals,
      "KS" -> ks,
      "N" -> n|>
  ];

(* --- 0d. Nearest-neighbor spacings --- *)
NNSpacings[specData_Association] :=
  Module[{s = specData["Spectrum"], n = specData["N"], diffs},
    diffs = Differences[s];
    (* include wrap-around spacing *)
    Append[diffs, n - Last[s] + First[s]]
  ];

(* --- 0e. Compiled binary-search level counter --- *)
CompileLevelCounter[extendedSpectrum_List] :=
  Module[{packed, m},
    packed = Developer`ToPackedArray[N[extendedSpectrum]];
    m = Length[packed];
    Compile[{{z, _Real}},
      Module[{lo = 1, hi = m, mid = 0},
        While[lo <= hi,
          mid = Quotient[lo + hi, 2];
          If[packed[[mid]] <= z, lo = mid + 1, hi = mid - 1]
        ];
        hi
      ],
      CompilationTarget -> "WVM", RuntimeOptions -> "Speed",
      RuntimeAttributes -> {Listable}, Parallelization -> False
    ]
  ];

(* --- 0f. Sigma2 via systematic sliding windows --- *)
(* Returns {L, Sigma2_theorMean, Sigma2_sampleMean} for each L *)
ComputeSigma2Sliding[counterFunc_, nLevels_Integer, lValues_List] :=
  Module[{validL, centers, nC},
    validL = Select[N[lValues], 0 < # <= (nLevels / 2.0) &];
    centers = N[Range[0, nLevels - 1]];  (* one center per level *)
    nC = Length[centers];
    Map[
      Function[{L},
        Module[{counts, s2theor, s2sample},
          counts = counterFunc[centers + L] - counterFunc[centers];
          s2theor  = Total[(counts - L)^2] / nC;  (* theoretical mean = L *)
          s2sample = Variance[counts];             (* sample mean subtracted *)
          {L, s2theor, s2sample}
        ]
      ],
      validL
    ]
  ];

(* --- 0g. Per-spectrum Sigma2 evaluation --- *)
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

(* --- 0h. Ensemble Sigma2 with dual estimator --- *)
EnsembleSigma2[spectraList_List, lValues_List] :=
  Module[{nMat, allResults, theorVals, sampleVals,
          meanTheor, seTheor, meanSample, seSample, effectiveL},
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
    <|"LValues"         -> effectiveL,
      "Sigma2TheorMean" -> meanTheor,
      "SETheorMean"     -> seTheor,
      "Sigma2SampleMean"-> meanSample,
      "SESampleMean"    -> seSample|>
  ];

(* --- 0i. Exact finite-N CUE number variance --- *)
CUESigma2Exact[L_?NumericQ, nDim_Integer] :=
  L - 2 NIntegrate[
    (L - r) (Sin[Pi r] / (nDim Sin[Pi r / nDim]))^2,
    {r, 0, L},
    Method -> "GaussKronrod", AccuracyGoal -> 10
  ];

PoissonSigma2[L_] := L;

(* --- 0j. Spectral form factor --- *)
CompileSFF =
  Compile[{{spec, _Real, 1}, {taus, _Real, 1}},
    Module[{n = Length[spec], cs, sn},
      Table[
        cs = Total[Cos[2.0 Pi spec t]];
        sn = Total[Sin[2.0 Pi spec t]];
        (cs^2 + sn^2) / n,
        {t, taus}
      ]
    ],
    CompilationTarget -> "WVM", RuntimeOptions -> "Speed"
  ];

EnsembleSFF[spectraList_List, tauMax_Real, nPts_Integer] :=
  Module[{nMat, taus, allK, meanK, seK},
    nMat = Length[spectraList];
    taus = N[Range[1, nPts] (tauMax / nPts)];
    DistributeDefinitions[CompileSFF];
    allK = ParallelTable[
      CompileSFF[spectraList[[i]]["Spectrum"], taus],
      {i, nMat}, Method -> "CoarsestGrained"
    ];
    meanK = Mean[allK];
    seK   = StandardDeviation[allK] / Sqrt[nMat * 1.0];
    <|"Tau" -> taus, "MeanK" -> meanK, "SE" -> seK|>
  ];

(* --- 0k. CUE SFF prediction (user normalization) --- *)
(* K_CUE(tau) = min(tau, 1) for tau > 0 *)
CUESFFPrediction[tau_] := Min[tau, 1];


(* ================================================================ *)
(* SECTION 1: PARAMETER SETUP & ENSEMBLE DIAGONALIZATION            *)
(* ================================================================ *)

J = 500;
alpha = 0.84;
delta = 0.01;
kCenter = 20.;

SeedRandom[123];
kList = RandomReal[{kCenter - delta, kCenter + delta}, 100];

Print["Diagonalizing ", Length[kList], " realizations..."];
AbsoluteTiming[
  rawPhases = ParallelTable[
    GetRawPhases[J, alpha, k], {k, kList},
    Method -> "CoarsestGrained"
  ];
]
ensembleSpectra = Map[UnfoldLinear, rawPhases];
effectiveDim = ensembleSpectra[[1]]["N"];
Print["Effective Hilbert space dimension: ", effectiveDim];


(* ================================================================ *)
(* SECTION 2: UNFOLDING VALIDATION                                  *)
(* ================================================================ *)

Print["\n=== STEP 2: UNFOLDING VALIDATION ==="];

(* --- 2a. Staircase vs linear prediction for 3 sample realizations --- *)
validationData = Map[ValidateUnfolding, ensembleSpectra[[{1, 50, 100}]]];

pltStaircase = Show[
  Table[
    ListPlot[
      Transpose[{validationData[[i]]["Phases"], validationData[[i]]["EmpiricalCDF"]}],
      PlotStyle -> {Opacity[0.6], {Blue, Red, Darker[Green]}[[i]]},
      Joined -> True,
      PlotRange -> All
    ],
    {i, 3}
  ],
  Plot[x / (2 Pi), {x, 0, 2 Pi},
    PlotStyle -> {Black, Dashed, Thick}
  ],
  FrameLabel -> {"\[Theta]", "N(\[Theta]) / N"},
  PlotLabel -> "Staircase vs Linear Prediction",
  Frame -> True, ImageSize -> 500
];
Print[pltStaircase];

(* --- 2b. Residuals in units of mean spacing --- *)
pltResiduals = Show[
  Table[
    ListLinePlot[
      Transpose[{validationData[[i]]["Phases"], validationData[[i]]["Residuals"]}],
      PlotStyle -> {{Opacity[0.7], {Blue, Red, Darker[Green]}[[i]]}},
      PlotRange -> All
    ],
    {i, 3}
  ],
  FrameLabel -> {"\[Theta]", "Residual (mean spacings)"},
  PlotLabel -> "Unfolding Residuals",
  Frame -> True, ImageSize -> 500
];
Print[pltResiduals];

(* --- 2c. KS statistic distribution across full ensemble --- *)
allKS = Map[ValidateUnfolding[#]["KS"] &, ensembleSpectra];
pltKS = Histogram[allKS, 20,
  FrameLabel -> {"KS / \!\(\*SqrtBox[\(N\)]\)", "Count"},
  PlotLabel -> "KS Statistic Distribution (ensemble)",
  Frame -> True, ImageSize -> 400
];
Print[pltKS];
Print["KS statistics: mean = ", Mean[allKS] // NumberForm[#, 4] &,
  ", max = ", Max[allKS] // NumberForm[#, 4] &];


(* ================================================================ *)
(* SECTION 3: NEAREST-NEIGHBOR SPACING DISTRIBUTION                 *)
(* ================================================================ *)

Print["\n=== STEP 3: NNSD CHECK ==="];

(* Pool spacings from all realizations *)
allSpacings = Flatten[Map[NNSpacings, ensembleSpectra]];

(* CUE Wigner surmise (beta=2) *)
wignerCUE[s_] := (32 / Pi^2) s^2 Exp[-4 s^2 / Pi];

pltNNSD = Show[
  Histogram[allSpacings, {0, 4, 0.05}, "ProbabilityDensity",
    PlotRange -> {{0, 4}, {0, 1.2}},
    ChartStyle -> Directive[LightBlue, EdgeForm[Gray]]
  ],
  Plot[wignerCUE[s], {s, 0, 4},
    PlotStyle -> {Red, Thick},
    PlotLegends -> {"CUE surmise"}
  ],
  Plot[Exp[-s], {s, 0, 4},
    PlotStyle -> {Black, Dashed},
    PlotLegends -> {"Poisson"}
  ],
  FrameLabel -> {"s", "P(s)"},
  PlotLabel -> "Nearest-Neighbor Spacing Distribution",
  Frame -> True, ImageSize -> 500
];
Print[pltNNSD];

(* Quantitative: mean and variance of spacings *)
Print["<s> = ", Mean[allSpacings] // NumberForm[#, 5] &,
  " (expect 1.0)"];
Print["Var(s) = ", Variance[allSpacings] // NumberForm[#, 5] &,
  " (CUE: ", N[1 - 8 / Pi^2 + 16 / Pi^2 - 1] // NumberForm[#, 4] &, ")"];
(* Exact CUE variance of spacing: (4 - 128/(9*Pi^2) + ...) ~ 0.178 *)


(* ================================================================ *)
(* SECTION 4: NUMBER VARIANCE - COMPUTATION & DUAL ESTIMATOR CHECK  *)
(* ================================================================ *)

Print["\n=== STEP 4: NUMBER VARIANCE (DUAL ESTIMATOR) ==="];

lDomain = Range[0.5, 100.0, 0.1];

AbsoluteTiming[
  sigma2Result = EnsembleSigma2[ensembleSpectra, lDomain];
]

(* --- 4a. Dual estimator comparison --- *)
(* If unfolding is correct, these should agree *)
pltDualEstimator = Show[
  ListLinePlot[
    {Transpose[{sigma2Result["LValues"], sigma2Result["Sigma2TheorMean"]}],
     Transpose[{sigma2Result["LValues"], sigma2Result["Sigma2SampleMean"]}]},
    PlotStyle -> {{Blue, Thick}, {Red, Dashed}},
    PlotLegends -> {"Theor mean (L)", "Sample mean"},
    PlotRange -> All
  ],
  FrameLabel -> {"L", "\!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\)(L)"},
  PlotLabel -> "Dual Estimator Diagnostic",
  Frame -> True, ImageSize -> 500
];
Print[pltDualEstimator];

(* Relative difference *)
relDiff = Abs[sigma2Result["Sigma2TheorMean"] - sigma2Result["Sigma2SampleMean"]] /
  (sigma2Result["Sigma2TheorMean"] + 10^-10);
pltRelDiff = ListLinePlot[
  Transpose[{sigma2Result["LValues"], relDiff}],
  PlotRange -> {0, All},
  FrameLabel -> {"L", "| \[CapitalDelta]\!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\) | / \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\)"},
  PlotLabel -> "Relative Difference: Theor vs Sample Mean",
  Frame -> True, Filling -> Axis,
  FillingStyle -> Directive[Opacity[0.3], Blue],
  ImageSize -> 500
];
Print[pltRelDiff];
Print["Max relative difference: ", Max[relDiff] // NumberForm[#, 4] &];


(* ================================================================ *)
(* SECTION 5: SIGMA2 vs CUE vs POISSON                             *)
(* ================================================================ *)

Print["\n=== STEP 5: NUMBER VARIANCE vs ANALYTICAL PREDICTIONS ==="];

(* Compute CUE prediction at a subset of L values (slow for large L) *)
lSubset = Range[0.5, 100.0, 0.5];
Print["Computing exact CUE Sigma2 for N = ", effectiveDim, "..."];
AbsoluteTiming[
  cueSigma2Values = Map[CUESigma2Exact[#, effectiveDim] &, lSubset];
]

pltSigma2 = Show[
  (* Numerical data with error band *)
  ListLinePlot[
    {Transpose[{sigma2Result["LValues"],
      sigma2Result["Sigma2TheorMean"] + sigma2Result["SETheorMean"]}],
     Transpose[{sigma2Result["LValues"],
      sigma2Result["Sigma2TheorMean"] - sigma2Result["SETheorMean"]}]},
    PlotStyle -> None,
    Filling -> {1 -> {2}},
    FillingStyle -> Directive[Opacity[0.2], Blue]
  ],
  ListLinePlot[
    Transpose[{sigma2Result["LValues"], sigma2Result["Sigma2TheorMean"]}],
    PlotStyle -> {Blue, Thick},
    PlotLegends -> {"Floquet data"}
  ],
  (* CUE exact *)
  ListLinePlot[
    Transpose[{lSubset, cueSigma2Values}],
    PlotStyle -> {Red, Thick},
    PlotLegends -> {"CUE(N) exact"}
  ],
  (* Poisson *)
  Plot[L, {L, 0, 100},
    PlotStyle -> {Black, Dashed, Thin},
    PlotLegends -> {"Poisson"}
  ],
  FrameLabel -> {"L", "\!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\)(L)"},
  PlotLabel -> "Number Variance",
  Frame -> True, PlotRange -> {0, All},
  ImageSize -> 600
];
Print[pltSigma2];

(* --- 5b. Residual: data minus CUE --- *)
(* Interpolate CUE for comparison *)
cueSigma2Interp = Interpolation[Transpose[{lSubset, cueSigma2Values}]];
lCommon = Select[sigma2Result["LValues"],
  lSubset[[1]] <= # <= lSubset[[-1]] &];
idxCommon = Flatten[Position[sigma2Result["LValues"], #] & /@ lCommon];
residualVsCUE = sigma2Result["Sigma2TheorMean"][[idxCommon]] -
  Map[cueSigma2Interp, lCommon];
seCommon = sigma2Result["SETheorMean"][[idxCommon]];

pltResidualCUE = Show[
  ListLinePlot[
    {Transpose[{lCommon, residualVsCUE + seCommon}],
     Transpose[{lCommon, residualVsCUE - seCommon}]},
    PlotStyle -> None, Filling -> {1 -> {2}},
    FillingStyle -> Directive[Opacity[0.3], Blue]
  ],
  ListLinePlot[
    Transpose[{lCommon, residualVsCUE}],
    PlotStyle -> {Blue}, PlotRange -> All
  ],
  Plot[0, {x, 0, 100}, PlotStyle -> {Red, Dashed}],
  FrameLabel -> {"L",
    "\!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\)(data) - \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\)(CUE)"},
  PlotLabel -> "Deviation from CUE",
  Frame -> True, ImageSize -> 500
];
Print[pltResidualCUE];


(* ================================================================ *)
(* SECTION 6: L-RANGE STABILITY CHECK                              *)
(* ================================================================ *)

Print["\n=== STEP 6: L-RANGE FINITE SIZE CHECK ==="];

(* Sigma2/Sigma2_CUE ratio — should be ~1 if CUE, deviating at large L *)
ratioToCUE = sigma2Result["Sigma2TheorMean"][[idxCommon]] /
  (Map[cueSigma2Interp, lCommon] + 10^-10);

pltRatio = Show[
  ListLinePlot[
    Transpose[{lCommon, ratioToCUE}],
    PlotStyle -> {Blue, Thick},
    PlotRange -> {0.8, 1.2}
  ],
  Plot[1, {x, 0, 100}, PlotStyle -> {Red, Dashed}],
  (* Mark L = N/4 and L = N/2 boundaries *)
  Graphics[{Gray, Dashed,
    Line[{{effectiveDim / 4, 0}, {effectiveDim / 4, 2}}],
    Line[{{effectiveDim / 2, 0}, {effectiveDim / 2, 2}}],
    Text["N/4", {effectiveDim / 4, 1.18}],
    Text["N/2", {effectiveDim / 2, 1.18}]
  }],
  FrameLabel -> {"L",
    "\!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\)(data) / \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\)(CUE)"},
  PlotLabel -> "Ratio to CUE (finite-size diagnostic)",
  Frame -> True, ImageSize -> 500
];
Print[pltRatio];


(* ================================================================ *)
(* SECTION 7: SFF CROSS-VALIDATION                                 *)
(* ================================================================ *)

Print["\n=== STEP 7: SPECTRAL FORM FACTOR ==="];

AbsoluteTiming[
  sffResult = EnsembleSFF[ensembleSpectra, 3.0, 1000];
]

pltSFF = Show[
  (* Data with error band *)
  ListLinePlot[
    {Transpose[{sffResult["Tau"], sffResult["MeanK"] + sffResult["SE"]}],
     Transpose[{sffResult["Tau"], sffResult["MeanK"] - sffResult["SE"]}]},
    PlotStyle -> None, Filling -> {1 -> {2}},
    FillingStyle -> Directive[Opacity[0.2], Blue]
  ],
  ListLinePlot[
    Transpose[{sffResult["Tau"], sffResult["MeanK"]}],
    PlotStyle -> {Blue, Thick},
    PlotLegends -> {"Floquet data"}
  ],
  (* CUE prediction: min(tau, 1) *)
  Plot[Min[tau, 1], {tau, 0, 3},
    PlotStyle -> {Red, Thick, Dashed},
    PlotLegends -> {"CUE: min(\[Tau],1)"}
  ],
  (* Poisson: K = 1 *)
  Plot[1, {tau, 0, 3},
    PlotStyle -> {Black, Dotted},
    PlotLegends -> {"Poisson"}
  ],
  FrameLabel -> {"\[Tau]", "K(\[Tau])"},
  PlotLabel -> "Spectral Form Factor",
  Frame -> True,
  PlotRange -> {{0, 3}, {0, 1.5}},
  ImageSize -> 600
];
Print[pltSFF];

(* --- 7b. SFF residual vs CUE --- *)
cueSFFVals = Map[CUESFFPrediction, sffResult["Tau"]];
sffResidual = sffResult["MeanK"] - cueSFFVals;

pltSFFResidual = Show[
  ListLinePlot[
    Transpose[{sffResult["Tau"], sffResidual}],
    PlotStyle -> {Blue},
    Filling -> Axis,
    FillingStyle -> Directive[Opacity[0.3], Blue]
  ],
  Plot[0, {x, 0, 3}, PlotStyle -> {Red, Dashed}],
  FrameLabel -> {"\[Tau]", "K(data) - K(CUE)"},
  PlotLabel -> "SFF Deviation from CUE",
  Frame -> True, ImageSize -> 500
];
Print[pltSFFResidual];


(* ================================================================ *)
(* SUMMARY                                                          *)
(* ================================================================ *)

Print["\n=== DIAGNOSTIC SUMMARY ==="];
Print["Dimension N = ", effectiveDim];
Print["Ensemble size = ", Length[ensembleSpectra]];
Print["L range = [", lDomain[[1]], ", ", lDomain[[-1]], "]"];
Print["L/N max = ", lDomain[[-1]] / effectiveDim // NumberForm[#, 3] &];
Print["Mean KS = ", Mean[allKS] // NumberForm[#, 4] &];
Print["Dual estimator max relative diff = ", Max[relDiff] // NumberForm[#, 4] &];
