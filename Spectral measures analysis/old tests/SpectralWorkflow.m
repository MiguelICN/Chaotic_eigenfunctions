(* ================================================================ *)
(* SPECTRAL DIAGNOSTICS: COMPLETE WORKFLOW                          *)
(* Path-independent. Place this file and QKT_full.wl in same folder *)
(* ================================================================ *)


(* ================================================================ *)
(* MODULE 0: SETUP                                                  *)
(* ================================================================ *)

baseDir = NotebookDirectory[];
Get[FileNameJoin[{baseDir, "QKT_full.wl"}]];
LaunchKernels[];


(* ================================================================ *)
(* MODULE 1: GLOBAL PARAMETERS                                     *)
(* ================================================================ *)

J       = 500;
alpha   = 0.84;
kCenter = 27.5;
delta   = 0.001;
nQKT    = 500;
nCOE    = 500;
Lmax    = 250;
lDomain = Join[Range[0.1, 5, 0.1], Range[5.5, Lmax, 0.5]];
tauMax  = 3.0;
nPtsSFF = 1000;

(* --- Directories --- *)
EnsureDir[dir_] := If[!DirectoryQ[dir],
  CreateDirectory[dir, CreateIntermediateDirectories -> True]];

runLabel = "J" <> ToString[J] <> "_k" <> ToString[N[kCenter]];
resultsDir = FileNameJoin[{baseDir, "Results", runLabel}];
evenDir    = FileNameJoin[{resultsDir, "Even"}];
oddDir     = FileNameJoin[{resultsDir, "Odd"}];
plotDir    = FileNameJoin[{resultsDir, "plots"}];
EnsureDir[evenDir]; EnsureDir[oddDir]; EnsureDir[plotDir];

(* --- Plot style --- *)
$pOpts = Sequence[ImageSize -> 800, Frame -> True,
  FrameStyle -> Directive[30, Black],
  LabelStyle -> Directive[30, Black]];
$legStyle = Directive[Black, 18];


(* ================================================================ *)
(* MODULE 2: HELPER PLOT FUNCTIONS                                  *)
(* ================================================================ *)

(* --- Staircase + residuals --- *)
PlotStaircase[specData_, label_] :=
  Module[{vd = ValidateUnfolding[specData]},
    {Show[
      ListLinePlot[Transpose[{vd["Phases"], vd["EmpiricalCDF"]}],
        PlotStyle -> {Blue, Thick}],
      Plot[x / (2 Pi), {x, 0, 2 Pi}, PlotStyle -> {Black, Dashed, Thick}],
      $pOpts, FrameLabel -> {"\[Theta]", "\[ScriptCapitalN](\[Theta])/N"},
      PlotLabel -> label <> ": Staircase"],
     Show[
      ListLinePlot[Transpose[{vd["Phases"], vd["Residuals"]}],
        PlotStyle -> {Blue, AbsoluteThickness[1]}],
      $pOpts, FrameLabel -> {"\[Theta]", "Residual (mean spacings)"},
      PlotLabel -> label <> ": Residuals",
      GridLines -> {None, {0}}, GridLinesStyle -> Directive[Red, Dashed]],
     vd["KS"]}
  ];

(* --- NNSD histogram --- *)
PlotNNSD[spacings_, label_] :=
  Legended[Show[
    Histogram[spacings, {0, 4, 0.05}, "ProbabilityDensity",
      ChartStyle -> Directive[LightBlue, EdgeForm[Gray]],
      PlotRange -> {{0, 3.5}, {0, 1.1}}],
    Plot[WignerCOE[s], {s, 0, 3.5}, PlotStyle -> {Red, Thick}],
    Plot[Exp[-s], {s, 0, 3.5}, PlotStyle -> {Black, Dashed}],
    $pOpts, FrameLabel -> {"s", "P(s)"}, PlotLabel -> label <> ": P(s)"],
  Placed[LineLegend[
    {Directive[Red, Thick], Directive[Black, Dashed]},
    {"COE surmise", "Poisson"}, LabelStyle -> $legStyle
  ], {Right, Bottom}]];

(* --- Cumulative I(s) with KL divergence --- *)
PlotCumulativeS[spacings_, label_] :=
  Module[{sorted, nSp, empCDF, coeCDF, ks, dkl},
    sorted = Sort[spacings]; nSp = Length[sorted];
    empCDF = Range[nSp] / nSp;
    coeCDF = IntWignerCOE /@ sorted;
    ks = Max[Abs[empCDF - coeCDF]];
    dkl = KLDivergence[spacings, WignerCOE];
    Legended[Show[
      ListLinePlot[Transpose[{sorted, empCDF}],
        PlotStyle -> {Blue, Thick}, PlotRange -> {{0, 3.5}, {0, 1}}],
      Plot[IntWignerCOE[s], {s, 0, 3.5}, PlotStyle -> {Red, Thick}],
      Plot[1 - Exp[-s], {s, 0, 3.5}, PlotStyle -> {Black, Dashed}],
      $pOpts, FrameLabel -> {"s", "I(s)"},
      PlotLabel -> label <> ": I(s)\nKS=" <>
        ToString[NumberForm[ks, 4]] <> ", D\[LetterSpace]KL=" <>
        ToString[NumberForm[dkl, 4]]],
    Placed[LineLegend[
      {Directive[Blue, Thick], Directive[Red, Thick], Directive[Black, Dashed]},
      {"Data", "COE", "Poisson"}, LabelStyle -> $legStyle
    ], {Right, Bottom}]]
  ];

(* --- P(r) histogram --- *)
PlotPrHist[ratios_, label_] :=
  Legended[Show[
    Histogram[ratios, {0, 1, 0.02}, "ProbabilityDensity",
      ChartStyle -> Directive[LightBlue, EdgeForm[Gray]],
      PlotRange -> {{0, 1}, {0, 2.2}}],
    Plot[PrCOE[r], {r, 0, 1}, PlotStyle -> {Red, Thick}],
    Plot[PrCUE[r], {r, 0, 1}, PlotStyle -> {Blue, Thick, Dashed}],
    Plot[PrPoisson[r], {r, 0, 1}, PlotStyle -> {Black, Dashed}],
    $pOpts, FrameLabel -> {"r", "P(r)"}, PlotLabel -> label <> ": P(r)"],
  Placed[LineLegend[
    {Directive[Red, Thick], Directive[Blue, Thick, Dashed], Directive[Black, Dashed]},
    {"COE (\[Beta]=1)", "CUE (\[Beta]=2)", "Poisson"}, LabelStyle -> $legStyle
  ], {Right, Bottom}]];

(* --- Cumulative I(r) with KL divergence --- *)
PlotCumulativeR[ratios_, label_] :=
  Module[{sorted, nR, empCDF, coeCDF, ks, dkl},
    sorted = Sort[ratios]; nR = Length[sorted];
    empCDF = Range[nR] / nR;
    (* Subsample for speed *)
    Module[{step = Max[1, Floor[nR / 1000]], idx, sR, eC, cC},
      idx = Range[1, nR, step];
      sR = sorted[[idx]]; eC = empCDF[[idx]];
      cC = IntPrCOE /@ sR;
      ks = Max[Abs[empCDF[[idx]] - cC]];
      dkl = KLDivergence[ratios, PrCOE, {0, 1, 0.02}];
      Legended[Show[
        ListLinePlot[Transpose[{sR, eC}],
          PlotStyle -> {Blue, Thick}, PlotRange -> {{0, 1}, {0, 1}}],
        ListLinePlot[Transpose[{sR, cC}], PlotStyle -> {Red, Thick}],
        $pOpts, FrameLabel -> {"r", "I(r)"},
        PlotLabel -> label <> ": I(r)\nKS=" <>
          ToString[NumberForm[ks, 4]] <> ", D\[LetterSpace]KL=" <>
          ToString[NumberForm[dkl, 4]]],
      Placed[LineLegend[
        {Directive[Blue, Thick], Directive[Red, Thick]},
        {"Data", "COE"}, LabelStyle -> $legStyle
      ], {Right, Bottom}]]
    ]
  ];

(* --- Sigma2 comparison: Even, Odd, COE --- *)
PlotSigma2[s2Even_, s2Odd_, s2COE_, label_] :=
  Legended[Show[
    (* COE band *)
    ListLinePlot[
      {Transpose[{s2COE["LValues"],
        s2COE["Sigma2TheorMean"] + s2COE["SETheorMean"]}],
       Transpose[{s2COE["LValues"],
        s2COE["Sigma2TheorMean"] - s2COE["SETheorMean"]}]},
      PlotStyle -> None, Filling -> {1 -> {2}},
      FillingStyle -> Directive[Opacity[0.15], Red]],
    ListLinePlot[Transpose[{s2COE["LValues"], s2COE["Sigma2TheorMean"]}],
      PlotStyle -> {Red, Thick}],
    (* Even band *)
    ListLinePlot[
      {Transpose[{s2Even["LValues"],
        s2Even["Sigma2TheorMean"] + s2Even["SETheorMean"]}],
       Transpose[{s2Even["LValues"],
        s2Even["Sigma2TheorMean"] - s2Even["SETheorMean"]}]},
      PlotStyle -> None, Filling -> {1 -> {2}},
      FillingStyle -> Directive[Opacity[0.15], Blue]],
    ListLinePlot[Transpose[{s2Even["LValues"], s2Even["Sigma2TheorMean"]}],
      PlotStyle -> {Blue, Thick}],
    (* Odd band *)
    ListLinePlot[
      {Transpose[{s2Odd["LValues"],
        s2Odd["Sigma2TheorMean"] + s2Odd["SETheorMean"]}],
       Transpose[{s2Odd["LValues"],
        s2Odd["Sigma2TheorMean"] - s2Odd["SETheorMean"]}]},
      PlotStyle -> None, Filling -> {1 -> {2}},
      FillingStyle -> Directive[Opacity[0.15], Darker[Green]]],
    ListLinePlot[Transpose[{s2Odd["LValues"], s2Odd["Sigma2TheorMean"]}],
      PlotStyle -> {Darker[Green], Thick}],
    (* Asymptotic + Poisson *)
    Plot[COEAsymptoticSigma2[L], {L, 0.5, Lmax},
      PlotStyle -> {Orange, Thick, Dashed}],
    Plot[L, {L, 0.5, Lmax}, PlotStyle -> {Black, Dashed, Thin}],
    $pOpts,
    FrameLabel -> {"L", "\!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\)(L)"},
    PlotLabel -> label, PlotRange -> {All, {0, Automatic}}],
  Placed[LineLegend[
    {Directive[Blue, Thick], Directive[Darker[Green], Thick],
     Directive[Red, Thick], Directive[Orange, Thick, Dashed],
     Directive[Black, Dashed, Thin]},
    {"Even", "Odd", "COE(" <> ToString[nCOE] <> ")",
     "COE asymptotic", "Poisson"}, LabelStyle -> $legStyle
  ], {Right, Bottom}]];

(* --- Residual: Even and Odd minus COE --- *)
PlotResidual[s2Even_, s2Odd_, s2COE_, label_] :=
  Module[{nMin, lR, resE, resO, seE, seO},
    nMin = Min[Length[s2Even["LValues"]], Length[s2Odd["LValues"]],
               Length[s2COE["LValues"]]];
    lR = s2COE["LValues"][[1 ;; nMin]];
    resE = s2Even["Sigma2TheorMean"][[1 ;; nMin]] -
      s2COE["Sigma2TheorMean"][[1 ;; nMin]];
    resO = s2Odd["Sigma2TheorMean"][[1 ;; nMin]] -
      s2COE["Sigma2TheorMean"][[1 ;; nMin]];
    seE = Sqrt[s2Even["SETheorMean"][[1 ;; nMin]]^2 +
      s2COE["SETheorMean"][[1 ;; nMin]]^2];
    seO = Sqrt[s2Odd["SETheorMean"][[1 ;; nMin]]^2 +
      s2COE["SETheorMean"][[1 ;; nMin]]^2];
    Legended[Show[
      ListLinePlot[{Transpose[{lR, resE + seE}], Transpose[{lR, resE - seE}]},
        PlotStyle -> None, Filling -> {1 -> {2}},
        FillingStyle -> Directive[Opacity[0.2], Blue]],
      ListLinePlot[Transpose[{lR, resE}], PlotStyle -> {Blue, Thick}],
      ListLinePlot[{Transpose[{lR, resO + seO}], Transpose[{lR, resO - seO}]},
        PlotStyle -> None, Filling -> {1 -> {2}},
        FillingStyle -> Directive[Opacity[0.2], Darker[Green]]],
      ListLinePlot[Transpose[{lR, resO}], PlotStyle -> {Darker[Green], Thick}],
      Plot[0, {x, 0, Max[lR]}, PlotStyle -> {Red, Dashed}],
      $pOpts,
      FrameLabel -> {"L",
        "\!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\)(QKT) - \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\)(COE)"},
      PlotLabel -> label],
    Placed[LineLegend[
      {Directive[Blue, Thick], Directive[Darker[Green], Thick]},
      {"Even - COE", "Odd - COE"}, LabelStyle -> $legStyle
    ], {Right, Bottom}]]
  ];

(* --- Dual estimator consistency --- *)
PlotDualEstimator[s2_, label_] :=
  Module[{rd},
    rd = Abs[s2["Sigma2TheorMean"] - s2["Sigma2SampleMean"]] /
      (s2["Sigma2TheorMean"] + 10^-10);
    ListLinePlot[Transpose[{s2["LValues"], rd}],
      Filling -> Axis, FillingStyle -> Directive[Opacity[0.3], Blue],
      $pOpts, PlotRange -> {0, All},
      FrameLabel -> {"L", "Relative difference"},
      PlotLabel -> label <> ": Dual Estimator"]
  ];

(* --- Fluctuating staircase --- *)
PlotFluctuations[flucEven_, flucOdd_, label_] :=
  Legended[Show[
    ListLinePlot[Transpose[{flucEven["x"], flucEven["delta"]}],
      PlotStyle -> {Blue, AbsoluteThickness[0.5]}],
    ListLinePlot[Transpose[{flucOdd["x"], flucOdd["delta"]}],
      PlotStyle -> {Darker[Green], AbsoluteThickness[0.5]}],
    $pOpts,
    FrameLabel -> {"Unfolded level index x",
      "\!\(\*OverscriptBox[\(N\), \(~\)]\)(x)"},
    PlotLabel -> label <> ": Fluctuating Staircase",
    GridLines -> {None, {0}}, GridLinesStyle -> Directive[Red, Dashed]],
  Placed[LineLegend[
    {Directive[Blue, AbsoluteThickness[1.5]],
     Directive[Darker[Green], AbsoluteThickness[1.5]]},
    {"Even", "Odd"}, LabelStyle -> $legStyle
  ], {Right, Bottom}]];

(* --- Length spectrum --- *)
PlotLengthSpectrum[lsEven_, lsOdd_, label_, nShow_: 50] :=
  Legended[Show[
    ListLinePlot[
      Transpose[{lsEven["Period"][[1 ;; nShow]],
                  lsEven["Amplitude"][[1 ;; nShow]]}],
      PlotStyle -> {Blue, AbsoluteThickness[1.5]}],
    ListLinePlot[
      Transpose[{lsOdd["Period"][[1 ;; nShow]],
                  lsOdd["Amplitude"][[1 ;; nShow]]}],
      PlotStyle -> {Red, AbsoluteThickness[1.5]}],
    $pOpts, FrameLabel -> {"Period (kicks)", "Length Spectrum"},
    PlotLabel -> label <> ": Length Spectrum",
    GridLines -> {None, {0}}, GridLinesStyle -> Directive[Gray, Dashed]],
  Placed[LineLegend[
    {Directive[Blue, AbsoluteThickness[1.5]],
     Directive[Red, AbsoluteThickness[1.5]]},
    {"Even", "Odd"}, LabelStyle -> $legStyle
  ], {Right, Bottom}]];

(* --- SFF --- *)
PlotSFF[sffEvenEns_, sffOddEns_, sffCOE_, label_] :=
  Legended[Show[
    (* COE band *)
    ListLinePlot[
      {Transpose[{sffCOE["Tau"], sffCOE["MeanK"] + sffCOE["SE"]}],
       Transpose[{sffCOE["Tau"], sffCOE["MeanK"] - sffCOE["SE"]}]},
      PlotStyle -> None, Filling -> {1 -> {2}},
      FillingStyle -> Directive[Opacity[0.15], Red]],
    ListLinePlot[Transpose[{sffCOE["Tau"], sffCOE["MeanK"]}],
      PlotStyle -> {Red, Thick}],
    (* Even band *)
    ListLinePlot[
      {Transpose[{sffEvenEns["Tau"], sffEvenEns["MeanK"] + sffEvenEns["SE"]}],
       Transpose[{sffEvenEns["Tau"], sffEvenEns["MeanK"] - sffEvenEns["SE"]}]},
      PlotStyle -> None, Filling -> {1 -> {2}},
      FillingStyle -> Directive[Opacity[0.15], Blue]],
    ListLinePlot[Transpose[{sffEvenEns["Tau"], sffEvenEns["MeanK"]}],
      PlotStyle -> {Blue, Thick}],
    (* Odd band *)
    ListLinePlot[
      {Transpose[{sffOddEns["Tau"], sffOddEns["MeanK"] + sffOddEns["SE"]}],
       Transpose[{sffOddEns["Tau"], sffOddEns["MeanK"] - sffOddEns["SE"]}]},
      PlotStyle -> None, Filling -> {1 -> {2}},
      FillingStyle -> Directive[Opacity[0.15], Darker[Green]]],
    ListLinePlot[Transpose[{sffOddEns["Tau"], sffOddEns["MeanK"]}],
      PlotStyle -> {Darker[Green], Thick}],
    (* Analytical *)
    Plot[COESFFAnalytical[tau], {tau, 0.001, tauMax},
      PlotStyle -> {Black, Thick, Dashed}],
    $pOpts, FrameLabel -> {"\[Tau]", "K(\[Tau])"},
    PlotLabel -> label <> ": SFF",
    PlotRange -> {{0, tauMax}, {0, 1.5}}],
  Placed[LineLegend[
    {Directive[Blue, Thick], Directive[Darker[Green], Thick],
     Directive[Red, Thick], Directive[Black, Thick, Dashed]},
    {"Even", "Odd", "COE(" <> ToString[nCOE] <> ")", "COE analytical"},
    LabelStyle -> $legStyle
  ], {Right, Bottom}]];

(* --- IPR distribution --- *)
PlotIPR[iprsEven_, iprsOdd_, nDim_, label_] :=
  Legended[Show[
    Histogram[{iprsEven, iprsOdd}, Automatic, "ProbabilityDensity",
      ChartStyle -> {Directive[Blue, Opacity[0.5]],
                     Directive[Darker[Green], Opacity[0.5]]}],
    (* COE prediction: IPR ~ 3/N *)
    Graphics[{Red, Dashed, AbsoluteThickness[2],
      Line[{{3.0 / nDim, 0}, {3.0 / nDim, 100}}]}],
    $pOpts, FrameLabel -> {"IPR", "Density"},
    PlotLabel -> label <> ": IPR Distribution",
    PlotRange -> {All, All}],
  Placed[LineLegend[
    {Directive[Blue, Opacity[0.5]], Directive[Darker[Green], Opacity[0.5]],
     Directive[Red, Dashed]},
    {"Even", "Odd", "COE: 3/N"}, LabelStyle -> $legStyle
  ], {Right, Bottom}]];


(* === Export helper === *)
SavePlot[plt_, name_] := (
  Export[FileNameJoin[{plotDir, name <> ".pdf"}], plt];
  Print["  Saved: ", name, ".pdf"];
);


(* ================================================================ *)
(* MODULE 3: COE REFERENCE                                          *)
(* ================================================================ *)

Print["=== COE REFERENCE ==="];
dimCOE = J + 1;

Print["Generating ", nCOE, " COE(", dimCOE, ") matrices..."];
{tCOE, coeSpectra} = AbsoluteTiming[GenerateCOEEnsemble[dimCOE, nCOE]];
Print["  Done in ", Round[tCOE, 0.1], "s"];

(* --- COE unfolding validation --- *)
{pltSt, pltRes, ksVal} = PlotStaircase[coeSpectra[[1]], "COE"];
Print["  KS = ", NumberForm[ksVal, 4]];
SavePlot[pltSt, "coe_staircase"];
SavePlot[pltRes, "coe_residuals"];

(* --- COE NNSD --- *)
spacingsCOE = PoolSpacings[coeSpectra];
SavePlot[PlotNNSD[spacingsCOE, "COE"], "coe_nnsd"];
SavePlot[PlotCumulativeS[spacingsCOE, "COE"], "coe_cumulative_s"];

(* --- COE P(r) --- *)
ratiosCOE = PoolSpacingRatios[coeSpectra];
SavePlot[PlotPrHist[ratiosCOE, "COE"], "coe_pr"];
SavePlot[PlotCumulativeR[ratiosCOE, "COE"], "coe_cumulative_r"];
Print["  COE <r> = ", NumberForm[Mean[ratiosCOE], 5]];

(* --- COE Sigma2 --- *)
Print["Computing COE Sigma2..."];
{tS2, sigma2COE} = AbsoluteTiming[EnsembleSigma2[coeSpectra, lDomain]];
Print["  Done in ", Round[tS2, 0.1], "s"];

SavePlot[PlotDualEstimator[sigma2COE, "COE"], "coe_dual_estimator"];
relDiffCOE = Max[Abs[sigma2COE["Sigma2TheorMean"] - sigma2COE["Sigma2SampleMean"]] /
  (sigma2COE["Sigma2TheorMean"] + 10^-10)];
Print["  COE max dual diff = ", NumberForm[relDiffCOE, 4]];

(* --- COE SFF --- *)
Print["Computing COE SFF..."];
{tSFF, sffCOE} = AbsoluteTiming[EnsembleSFF[coeSpectra, tauMax, nPtsSFF]];
Print["  Done in ", Round[tSFF, 0.1], "s"];

(* Free COE spectra, keep summary data *)
coeSpectra = .;
Print["COE REFERENCE COMPLETE\n"];


(* ================================================================ *)
(* MODULE 4: QKT ENSEMBLE                                          *)
(* ================================================================ *)

Print["=== QKT ENSEMBLE ==="];
SeedRandom[123];
kList = RandomReal[{kCenter - delta, kCenter + delta}, nQKT];

Print["Diagonalizing ", nQKT, " QKT matrices (J=", J, ")..."];
{tDiag, qktEnsemble} = AbsoluteTiming[
  GenerateQKTEnsemble[J, alpha, kList]];
Print["  Done in ", Round[tDiag, 0.1], "s"];

dimE = qktEnsemble["Even"][[1]]["N"];
dimO = qktEnsemble["Odd"][[1]]["N"];
Print["  Even N = ", dimE, ", Odd N = ", dimO];

(* --- Export spectra --- *)
Export[FileNameJoin[{evenDir, "spectra.wdx"}], qktEnsemble["Even"], "WDX"];
Export[FileNameJoin[{oddDir, "spectra.wdx"}], qktEnsemble["Odd"], "WDX"];
Print["  Spectra exported"];


(* ================================================================ *)
(* MODULE 5: UNFOLDING VALIDATION                                   *)
(* ================================================================ *)

Print["\n=== UNFOLDING VALIDATION ==="];
Do[
  specList = qktEnsemble[sec];
  nR = Length[specList];
  sIdx = Round[Subdivide[1, nR, 2]];
  vd = Map[ValidateUnfolding, specList[[sIdx]]];
  allKS = Map[ValidateUnfolding[#]["KS"] &, specList];

  Print["--- ", sec, " ---"];
  Print["  KS: mean=", NumberForm[Mean[allKS], 4],
    ", max=", NumberForm[Max[allKS], 4]];

  SavePlot[Show[
    Table[ListLinePlot[
      Transpose[{vd[[i]]["Phases"], vd[[i]]["EmpiricalCDF"]}],
      PlotStyle -> {{Blue, Red, Darker[Green]}[[i]], Opacity[0.6]}
    ], {i, 3}],
    Plot[x / (2 Pi), {x, 0, 2 Pi}, PlotStyle -> {Black, Dashed, Thick}],
    $pOpts, FrameLabel -> {"\[Theta]", "\[ScriptCapitalN](\[Theta])/N"},
    PlotLabel -> sec <> ": Staircase"
  ], sec <> "_staircase"];

  SavePlot[Show[
    Table[ListLinePlot[
      Transpose[{vd[[i]]["Phases"], vd[[i]]["Residuals"]}],
      PlotStyle -> {{Blue, Red, Darker[Green]}[[i]], Opacity[0.7]}
    ], {i, 3}],
    $pOpts, FrameLabel -> {"\[Theta]", "Residual (mean spacings)"},
    PlotLabel -> sec <> ": Residuals",
    GridLines -> {None, {0}}, GridLinesStyle -> Directive[Red, Dashed]
  ], sec <> "_residuals"];
, {sec, {"Even", "Odd"}}];


(* ================================================================ *)
(* MODULE 6: NNSD + P(r)                                            *)
(* ================================================================ *)

Print["\n=== SPACING STATISTICS ==="];
Do[
  specList = qktEnsemble[sec];
  spacings = PoolSpacings[specList];
  ratios   = PoolSpacingRatios[specList];

  Print["--- ", sec, " ---"];
  Print["  <s>=", NumberForm[Mean[spacings], 5],
    ", Var(s)=", NumberForm[Variance[spacings], 5],
    "  (COE: ", NumberForm[N[4/Pi - 1], 4], ")"];
  Print["  <r>=", NumberForm[Mean[ratios], 5],
    "  (COE: 0.5359)"];

  SavePlot[PlotNNSD[spacings, sec], sec <> "_nnsd"];
  SavePlot[PlotCumulativeS[spacings, sec], sec <> "_cumulative_s"];
  SavePlot[PlotPrHist[ratios, sec], sec <> "_pr"];
  SavePlot[PlotCumulativeR[ratios, sec], sec <> "_cumulative_r"];

  Export[FileNameJoin[{If[sec == "Even", evenDir, oddDir],
    "spacings.wdx"}], spacings, "WDX"];
  Export[FileNameJoin[{If[sec == "Even", evenDir, oddDir],
    "ratios.wdx"}], ratios, "WDX"];
, {sec, {"Even", "Odd"}}];


(* ================================================================ *)
(* MODULE 7: NUMBER VARIANCE                                        *)
(* ================================================================ *)

Print["\n=== NUMBER VARIANCE ==="];

Print["Computing Even Sigma2..."];
{t1, sigma2Even} = AbsoluteTiming[EnsembleSigma2[qktEnsemble["Even"], lDomain]];
Print["  Done in ", Round[t1, 0.1], "s"];

Print["Computing Odd Sigma2..."];
{t2, sigma2Odd} = AbsoluteTiming[EnsembleSigma2[qktEnsemble["Odd"], lDomain]];
Print["  Done in ", Round[t2, 0.1], "s"];

Export[FileNameJoin[{evenDir, "sigma2.wdx"}], sigma2Even, "WDX"];
Export[FileNameJoin[{oddDir, "sigma2.wdx"}], sigma2Odd, "WDX"];

mainLabel = "\!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\): J=" <>
  ToString[J] <> ", k=" <> ToString[N[kCenter]];
SavePlot[PlotSigma2[sigma2Even, sigma2Odd, sigma2COE, mainLabel], "sigma2"];
SavePlot[PlotResidual[sigma2Even, sigma2Odd, sigma2COE,
  "Residual: J=" <> ToString[J] <> ", k=" <> ToString[N[kCenter]]], "residual"];

Do[
  s2 = If[sec == "Even", sigma2Even, sigma2Odd];
  SavePlot[PlotDualEstimator[s2, sec], sec <> "_dual_estimator"];
, {sec, {"Even", "Odd"}}];

(* Parity-averaged residual *)
Module[{nMin, lR, avgS2, avgCOE, resAvg},
  nMin = Min[Length[sigma2Even["LValues"]], Length[sigma2Odd["LValues"]],
             Length[sigma2COE["LValues"]]];
  lR = sigma2COE["LValues"][[1 ;; nMin]];
  avgS2 = (dimE sigma2Even["Sigma2TheorMean"][[1 ;; nMin]] +
           dimO sigma2Odd["Sigma2TheorMean"][[1 ;; nMin]]) / (dimE + dimO);
  resAvg = avgS2 - sigma2COE["Sigma2TheorMean"][[1 ;; nMin]];
  SavePlot[Legended[Show[
    ListLinePlot[Transpose[{lR, resAvg}],
      PlotStyle -> {Purple, Thick},
      Filling -> Axis, FillingStyle -> Directive[Opacity[0.2], Purple]],
    Plot[0, {x, 0, Max[lR]}, PlotStyle -> {Red, Dashed}],
    $pOpts,
    FrameLabel -> {"L",
      "\!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\)(avg) - \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\)(COE)"},
    PlotLabel -> "Parity-Averaged Residual"],
  Placed[LineLegend[{Directive[Purple, Thick]},
    {"Weighted avg"}, LabelStyle -> $legStyle], {Right, Bottom}]
  ], "residual_parity_avg"];
  Print["  Max |avg residual| = ", NumberForm[Max[Abs[resAvg]], 4]];
];


(* ================================================================ *)
(* MODULE 8: SPECTRAL FORM FACTOR                                   *)
(* ================================================================ *)

Print["\n=== SFF ==="];

Print["Computing QKT SFF..."];
{t1, sffEvenEns} = AbsoluteTiming[EnsembleSFF[qktEnsemble["Even"], tauMax, nPtsSFF]];
{t2, sffOddEns}  = AbsoluteTiming[EnsembleSFF[qktEnsemble["Odd"], tauMax, nPtsSFF]];
Print["  Done in ", Round[t1 + t2, 0.1], "s"];

SavePlot[PlotSFF[sffEvenEns, sffOddEns, sffCOE,
  "J=" <> ToString[J] <> ", k=" <> ToString[N[kCenter]]], "sff"];


(* ================================================================ *)
(* MODULE 9: FLUCTUATING STAIRCASE + LENGTH SPECTRUM                *)
(* ================================================================ *)

Print["\n=== FLUCTUATION ANALYSIS ==="];

(* Use first realization *)
flucEven = GetFluctuations[qktEnsemble["Even"][[1]]];
flucOdd  = GetFluctuations[qktEnsemble["Odd"][[1]]];

SavePlot[PlotFluctuations[flucEven, flucOdd,
  "J=" <> ToString[J] <> ", k=" <> ToString[N[kCenter]]], "fluctuations"];

lsEven = LengthSpectrum[flucEven];
lsOdd  = LengthSpectrum[flucOdd];

SavePlot[PlotLengthSpectrum[lsEven, lsOdd,
  "J=" <> ToString[J] <> ", k=" <> ToString[N[kCenter]], 50], "length_spectrum"];


(* ================================================================ *)
(* MODULE 10: EIGENVECTOR STATISTICS                                *)
(* ================================================================ *)

Print["\n=== EIGENVECTOR ANALYSIS ==="];

(* Get eigenvectors for one realization *)
Print["Computing eigensystem..."];
{tEig, eigData} = AbsoluteTiming[
  GetQKTParityEigensystem[J, alpha, kCenter]];
Print["  Done in ", Round[tEig, 0.1], "s"];

iprsEven = Map[IPR, eigData["Even"]["Vectors"]];
iprsOdd  = Map[IPR, eigData["Odd"]["Vectors"]];

Print["  Even: <IPR>=", NumberForm[Mean[iprsEven], 4],
  " (COE: ", NumberForm[3.0 / dimE, 4], ")"];
Print["  Odd:  <IPR>=", NumberForm[Mean[iprsOdd], 4],
  " (COE: ", NumberForm[3.0 / dimO, 4], ")"];

SavePlot[PlotIPR[iprsEven, iprsOdd, dimE,
  "J=" <> ToString[J] <> ", k=" <> ToString[N[kCenter]]], "ipr"];

Export[FileNameJoin[{evenDir, "ipr.wdx"}], iprsEven, "WDX"];
Export[FileNameJoin[{oddDir, "ipr.wdx"}], iprsOdd, "WDX"];


(* ================================================================ *)
(* MODULE 11: SUMMARY                                               *)
(* ================================================================ *)

Print["\n=== SUMMARY ==="];
Print["J = ", J, ", alpha = ", alpha, ", k = ", kCenter, " +/- ", delta];
Print["QKT ensemble = ", nQKT, ", COE ensemble = ", nCOE];
Print["Even N = ", dimE, ", Odd N = ", dimO, ", COE N = ", dimCOE];
Print["L range = [", lDomain[[1]], ", ", lDomain[[-1]], "]"];
Print["Results in: ", resultsDir];
Print["\n=== ALL DONE ==="];
