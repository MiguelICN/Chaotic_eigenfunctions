(* ::Package:: *)

(* ================================================================ *)
(* SpectralAnalysis.m                                                *)
(* Complete spectral diagnostics workflow for the QKT                *)
(* Modular: each section runs independently                          *)
(* Linux/Windows compatible paths via FileNameJoin                   *)
(* Place this file and QKT_full.wl in the same folder               *)
(* ================================================================ *)


(* ================================================================ *)
(* MODULE 0: SETUP                                                  *)
(* ================================================================ *)

baseDir = NotebookDirectory[];
Get[FileNameJoin[{baseDir, "QKT_full.wl"}]];
LaunchKernels[];
$HistoryLength = 0;

EnsureDir[dir_] := If[!DirectoryQ[dir],
  CreateDirectory[dir, CreateIntermediateDirectories -> True]];

kStr[k_] := StringReplace[
  If[StringEndsQ[#, "."], # <> "0", #] &[ToString[N[k]]], "." -> "_"];


(* ================================================================ *)
(* MODULE 1: PARAMETERS                                             *)
(* ================================================================ *)

J       = 500;
alpha   = 0.84;
kCenter = 27.5;
delta   = 0.01;
nQKT    = 500;
nCOE    = 500;
Lmax    = 250;
lDomain = Join[Range[0.1, 5, 0.1], Range[5.5, Lmax, 0.5]];
tauMax  = 3.0;
nPtsSFF = 1000;

(* --- Directories --- *)
runLabel = "J" <> ToString[J] <> "_k" <> kStr[kCenter];
resultsDir = FileNameJoin[{baseDir, "Results", runLabel}];
evenDir    = FileNameJoin[{resultsDir, "Even"}];
oddDir     = FileNameJoin[{resultsDir, "Odd"}];
plotDir    = FileNameJoin[{resultsDir, "plots"}];
collageDir = FileNameJoin[{resultsDir, "individual"}];
EnsureDir[evenDir]; EnsureDir[oddDir];
EnsureDir[plotDir]; EnsureDir[collageDir];

(* --- Plot style --- *)
$pOpts = Sequence[ImageSize -> 800, Frame -> True,
  FrameStyle -> Directive[Black, 25],
  LabelStyle -> Directive[Black, 25]];
$cOpts = Sequence[Frame -> True,
  FrameStyle -> Directive[Black, 16],
  LabelStyle -> Directive[Black, 16]];
$cSize = 380;
$legStyle = Directive[Black, 18];

SavePlot[plt_, name_] := (
  Export[FileNameJoin[{plotDir, name <> ".png"}],
    Rasterize[plt, ImageResolution -> 150, Background -> White]];
  Print["  Saved: ", name, ".png"]);


(* ================================================================ *)
(* MODULE 2: SPECTRUM GENERATION / IMPORT                           *)
(* Set generateNew = True to compute, False to import .wdx files    *)
(* ================================================================ *)

generateNew = True;

If[generateNew,
  (* --- Generate --- *)
  Print["=== GENERATING QKT ENSEMBLE ==="];
  SeedRandom[42];
  kList = RandomReal[{kCenter - delta, kCenter + delta}, nQKT];
  {tGen, qktEnsemble} = AbsoluteTiming[
    GenerateQKTEnsemble[J, alpha, kList]];
  Print["  ", nQKT, " matrices done in ", Round[tGen, 0.1], "s"];
  Export[FileNameJoin[{evenDir, "spectra.wdx"}], qktEnsemble["Even"], "WDX"];
  Export[FileNameJoin[{oddDir, "spectra.wdx"}], qktEnsemble["Odd"], "WDX"];
  Print["  Spectra exported"];
  ,
  (* --- Import --- *)
  Print["=== IMPORTING SPECTRA ==="];
  qktEnsemble = <|
    "Even" -> Import[FileNameJoin[{evenDir, "spectra.wdx"}]],
    "Odd"  -> Import[FileNameJoin[{oddDir, "spectra.wdx"}]]|>;
  If[FailureQ[qktEnsemble["Even"]] || FailureQ[qktEnsemble["Odd"]],
    Print["ERROR: Import failed. Check paths."]; Abort[]];
  Print["  Imported successfully"];
];

dimE = qktEnsemble["Even"][[1]]["N"];
dimO = qktEnsemble["Odd"][[1]]["N"];
nReal = Length[qktEnsemble["Even"]];
Print["  Even N=", dimE, ", Odd N=", dimO, ", nReal=", nReal];


(* ================================================================ *)
(* MODULE 3: LOCAL HELPERS                                          *)
(* IntPrCOEFast interpolator + GoF metrics + collage builders      *)
(* ================================================================ *)

Print["\n=== BUILDING HELPERS ==="];
IntPrCOEFast = Interpolation[
  Table[{r, NIntegrate[PrCOE[x], {x, 0, r}]}, {r, 0, 1, 0.001}]];
Print["  IntPrCOEFast[0.5] = ", IntPrCOEFast[0.5]];

(* --- GoF metrics from sorted data --- *)
KSFromSorted[empCDF_List, theoCDFvals_List] := Max[Abs[empCDF - theoCDFvals]];

CvMFromSorted[n_Integer, theoCDFvals_List] :=
  Module[{midEmp = N[(2 Range[n] - 1) / (2 n)]},
    Total[(midEmp - theoCDFvals)^2] / n + 1 / (12 n)];

MSDFromSorted[empCDF_List, theoCDFvals_List] := Mean[(empCDF - theoCDFvals)^2];

(* --- GoF for one spectrum (pure numerical, parallelizable) --- *)
ComputeGoFSingle[specData_Association] :=
  Module[{spacings, sp, ratios, sortS, sortR, nS, nR,
          empS, empR, theoS, theoR},
    spacings = NNSpacings[specData]; sortS = Sort[spacings]; nS = Length[sortS];
    empS = N[Range[nS] / nS]; theoS = N[IntWignerCOE /@ sortS];
    sp = Differences[specData["Spectrum"]];
    ratios = Min /@ Transpose[{Most[sp] / Rest[sp], Rest[sp] / Most[sp]}];
    sortR = Sort[ratios]; nR = Length[sortR];
    empR = N[Range[nR] / nR]; theoR = N[IntPrCOEFast /@ sortR];
    {KSFromSorted[empS, theoS], CvMFromSorted[nS, theoS], MSDFromSorted[empS, theoS],
     KSFromSorted[empR, theoR], CvMFromSorted[nR, theoR], MSDFromSorted[empR, theoR]}
  ];

(* --- Extract spacings + ratios --- *)
ExtractStats[specData_Association] :=
  Module[{spacings, sp, ratios},
    spacings = NNSpacings[specData];
    sp = Differences[specData["Spectrum"]];
    ratios = Min /@ Transpose[{Most[sp] / Rest[sp], Rest[sp] / Most[sp]}];
    <|"SortedS" -> Sort[spacings], "SortedR" -> Sort[ratios]|>];

(* --- Collage row for one sector --- *)
MakeRowPlots[sortedS_List, sortedR_List, flucData_Association, label_String] :=
  Module[{nS, nR, empS, empR, theoS, theoR, ksS, ksR, diffS, diffR},
    nS = Length[sortedS]; nR = Length[sortedR];
    empS = N[Range[nS] / nS]; empR = N[Range[nR] / nR];
    theoS = N[IntWignerCOE /@ sortedS];
    theoR = N[IntPrCOEFast /@ sortedR];
    ksS = Max[Abs[empS - theoS]]; ksR = Max[Abs[empR - theoR]];
    diffS = (empS - theoS)^2; diffR = (empR - theoR)^2;
    {(* P(s) *)
     Show[Histogram[sortedS, {0, 4, 0.05}, "ProbabilityDensity",
       ChartStyle -> Directive[LightGray, EdgeForm[Gray]], PlotRange -> {{0, 3.5}, {0, 1.1}}],
      Plot[WignerCOE[x], {x, 0, 3.5}, PlotStyle -> {Red, Dashed, Thick}],
      Plot[Exp[-x], {x, 0, 3.5}, PlotStyle -> {Black, Dashed}],
      $cOpts, FrameLabel -> {"s", "P(s)"}, PlotLabel -> label <> " P(s)", ImageSize -> $cSize],
     (* P(r) *)
     Show[Histogram[sortedR, {0, 1, 0.02}, "ProbabilityDensity",
       ChartStyle -> Directive[LightGray, EdgeForm[Gray]], PlotRange -> {{0, 1}, {0, 2.2}}],
      Plot[PrCOE[x], {x, 0, 1}, PlotStyle -> {Red, Dashed, Thick}],
      Plot[PrPoisson[x], {x, 0, 1}, PlotStyle -> {Black, Dashed}],
      $cOpts, FrameLabel -> {"r", "P(r)"}, PlotLabel -> label <> " P(r)", ImageSize -> $cSize],
     (* I(s) *)
     Show[ListLinePlot[Transpose[{sortedS, empS}], PlotStyle -> {Blue, Thick},
       PlotRange -> {{0, 3.5}, {0, 1}}],
      Plot[IntWignerCOE[x], {x, 0, 3.5}, PlotStyle -> {Red, Dashed}],
      Plot[1 - Exp[-x], {x, 0, 3.5}, PlotStyle -> {Black, Dashed}],
      $cOpts, FrameLabel -> {"s", "I(s)"},
      PlotLabel -> label <> " I(s) KS=" <> ToString[NumberForm[ksS, 4]], ImageSize -> $cSize],
     (* I(r) *)
     Show[ListLinePlot[Transpose[{sortedR, empR}], PlotStyle -> {Blue, Thick},
       PlotRange -> {{0, 1}, {0, 1}}],
      ListLinePlot[Transpose[{sortedR, theoR}], PlotStyle -> {Red, Dashed}],
      $cOpts, FrameLabel -> {"r", "I(r)"},
      PlotLabel -> label <> " I(r) KS=" <> ToString[NumberForm[ksR, 4]], ImageSize -> $cSize],
     (* CDF diff (s) *)
     ListLinePlot[Transpose[{sortedS, diffS}], PlotStyle -> {Blue, Thick},
      Filling -> Axis, FillingStyle -> Directive[Opacity[0.3], Blue],
      $cOpts, FrameLabel -> {"s", "(I-I_COE)^2"}, PlotLabel -> label <> " CDF diff (s)",
      PlotRange -> All, ImageSize -> $cSize],
     (* CDF diff (r) *)
     ListLinePlot[Transpose[{sortedR, diffR}], PlotStyle -> {Blue, Thick},
      Filling -> Axis, FillingStyle -> Directive[Opacity[0.3], Blue],
      $cOpts, FrameLabel -> {"r", "(I-I_COE)^2"}, PlotLabel -> label <> " CDF diff (r)",
      PlotRange -> All, ImageSize -> $cSize],
     (* Fluctuation *)
     ListLinePlot[Transpose[{flucData["x"], flucData["delta"]}],
      PlotStyle -> {Blue, AbsoluteThickness[0.5]}, $cOpts,
      FrameLabel -> {"x", "\!\(\*OverscriptBox[\(N\), \(~\)]\)(x)"},
      PlotLabel -> label <> " Fluctuation", ImageSize -> $cSize,
      GridLines -> {None, {0}}, GridLinesStyle -> Directive[Red, Dashed]]}
  ];


(* ================================================================ *)
(* MODULE 4: INDIVIDUAL REALIZATION COLLAGES                        *)
(* ================================================================ *)

Print["\n=== INDIVIDUAL COLLAGES ==="];

Do[
  Module[{stE, stO, flE, flO, colE, colO, collage},
    stE = ExtractStats[qktEnsemble["Even"][[i]]];
    stO = ExtractStats[qktEnsemble["Odd"][[i]]];
    flE = GetFluctuations[qktEnsemble["Even"][[i]]];
    flO = GetFluctuations[qktEnsemble["Odd"][[i]]];
    colE = MakeRowPlots[stE["SortedS"], stE["SortedR"], flE, "Even"];
    colO = MakeRowPlots[stO["SortedS"], stO["SortedR"], flO, "Odd"];
    collage = Grid[
      Prepend[Transpose[{colE, colO}],
        {Style["Even", Bold, 20], Style["Odd", Bold, 20]}],
      Spacings -> {1, 0.5}, Frame -> All, FrameStyle -> GrayLevel[0.85]];
    Export[FileNameJoin[{collageDir,
      "realization_" <> StringPadLeft[ToString[i], 4, "0"] <> ".png"}],
      Rasterize[collage, ImageResolution -> 150, Background -> White]];
    If[Mod[i, 50] == 0 || i == 1, Print["  Exported ", i, "/", nReal]];
  ];
, {i, nReal}];
Print["Individual collages complete"];


(* ================================================================ *)
(* MODULE 5: ENSEMBLE SHORT-RANGE STATISTICS                        *)
(* Pooled P(s), P(r), I(s), I(r) + ensemble collage                *)
(* ================================================================ *)

Print["\n=== ENSEMBLE SHORT-RANGE ==="];

Module[{spE, spO, rE, rO, flE, flO, colE, colO, ensCollage},
  spE = Sort[Flatten[Map[NNSpacings, qktEnsemble["Even"]]]];
  spO = Sort[Flatten[Map[NNSpacings, qktEnsemble["Odd"]]]];
  rE  = Sort[PoolSpacingRatios[qktEnsemble["Even"]]];
  rO  = Sort[PoolSpacingRatios[qktEnsemble["Odd"]]];
  Print["  Even: <s>=", NumberForm[Mean[spE], 5], ", <r>=", NumberForm[Mean[rE], 5]];
  Print["  Odd:  <s>=", NumberForm[Mean[spO], 5], ", <r>=", NumberForm[Mean[rO], 5]];
  flE = GetFluctuations[qktEnsemble["Even"][[1]]];
  flO = GetFluctuations[qktEnsemble["Odd"][[1]]];
  colE = MakeRowPlots[spE, rE, flE, "Even (ens)"];
  colO = MakeRowPlots[spO, rO, flO, "Odd (ens)"];
  ensCollage = Grid[
    Prepend[Transpose[{colE, colO}],
      {Style["Even (pooled " <> ToString[nReal] <> ")", Bold, 20],
       Style["Odd (pooled " <> ToString[nReal] <> ")", Bold, 20]}],
    Spacings -> {1, 0.5}, Frame -> All, FrameStyle -> GrayLevel[0.85]];
  Export[FileNameJoin[{plotDir, "ensemble_collage.png"}],
    Rasterize[ensCollage, ImageResolution -> 150, Background -> White]];
  Print["  Ensemble collage exported"];
];


(* ================================================================ *)
(* MODULE 6: COE REFERENCE + NUMBER VARIANCE + SFF                  *)
(* ================================================================ *)

Print["\n=== COE REFERENCE ==="];
dimCOE = dimE;
{tCOE, coeSpectra} = AbsoluteTiming[GenerateCOEEnsemble[dimCOE, nCOE]];
Print["  ", nCOE, " COE(", dimCOE, ") done in ", Round[tCOE, 0.1], "s"];

Print["Computing COE Sigma2..."];
{tS, sigma2COE} = AbsoluteTiming[EnsembleSigma2[coeSpectra, lDomain]];
Print["  Done in ", Round[tS, 0.1], "s"];

Print["Computing COE SFF..."];
{tF, sffCOE} = AbsoluteTiming[EnsembleSFF[coeSpectra, tauMax, nPtsSFF]];
Print["  Done in ", Round[tF, 0.1], "s"];
coeSpectra = .;

Print["\n=== QKT NUMBER VARIANCE ==="];
{t1, sigma2Even} = AbsoluteTiming[EnsembleSigma2[qktEnsemble["Even"], lDomain]];
{t2, sigma2Odd}  = AbsoluteTiming[EnsembleSigma2[qktEnsemble["Odd"], lDomain]];
Print["  Even: ", Round[t1, 0.1], "s, Odd: ", Round[t2, 0.1], "s"];

SavePlot[Legended[Show[
  ListLinePlot[Transpose[{sigma2COE["LValues"], sigma2COE["Sigma2TheorMean"]}], PlotStyle -> {Red, Thick}],
  ListLinePlot[Transpose[{sigma2Even["LValues"], sigma2Even["Sigma2TheorMean"]}], PlotStyle -> {Blue, Thick}],
  ListLinePlot[Transpose[{sigma2Odd["LValues"], sigma2Odd["Sigma2TheorMean"]}], PlotStyle -> {Darker[Green], Thick}],
  Plot[COEAsymptoticSigma2[L], {L, 0.5, Lmax}, PlotStyle -> {Orange, Thick, Dashed}],
  Plot[L, {L, 0.5, Lmax}, PlotStyle -> {Black, Dashed, Thin}],
  $pOpts, FrameLabel -> {"L", "\!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\)(L)"},
  PlotLabel -> "Number Variance: J=" <> ToString[J],
  PlotRange -> {All, {0, Automatic}}],
  Placed[LineLegend[{Directive[Blue, Thick], Directive[Darker[Green], Thick],
    Directive[Red, Thick], Directive[Orange, Thick, Dashed]},
    {"Even", "Odd", "COE", "COE asym."}, LabelStyle -> $legStyle], {Right, Bottom}]
], "sigma2"];

(* Residual *)
Module[{nMin, lR, resE, resO},
  nMin = Min[Length /@ {sigma2Even["LValues"], sigma2Odd["LValues"], sigma2COE["LValues"]}];
  lR = sigma2COE["LValues"][[1 ;; nMin]];
  resE = sigma2Even["Sigma2TheorMean"][[1 ;; nMin]] - sigma2COE["Sigma2TheorMean"][[1 ;; nMin]];
  resO = sigma2Odd["Sigma2TheorMean"][[1 ;; nMin]] - sigma2COE["Sigma2TheorMean"][[1 ;; nMin]];
  SavePlot[Legended[Show[
    ListLinePlot[Transpose[{lR, resE}], PlotStyle -> {Blue, Thick}],
    ListLinePlot[Transpose[{lR, resO}], PlotStyle -> {Darker[Green], Thick}],
    Plot[0, {x, 0, Max[lR]}, PlotStyle -> {Red, Dashed}],
    $pOpts, FrameLabel -> {"L", "\[CapitalDelta]\!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\)(L)"},
    PlotLabel -> "Residual: QKT - COE"],
    Placed[LineLegend[{Directive[Blue, Thick], Directive[Darker[Green], Thick]},
      {"Even", "Odd"}, LabelStyle -> $legStyle], {Right, Bottom}]
  ], "residual"];
];

Export[FileNameJoin[{evenDir, "sigma2.wdx"}], sigma2Even, "WDX"];
Export[FileNameJoin[{oddDir, "sigma2.wdx"}], sigma2Odd, "WDX"];

Print["\n=== SFF ==="];
{t1, sffEven} = AbsoluteTiming[EnsembleSFF[qktEnsemble["Even"], tauMax, nPtsSFF]];
{t2, sffOdd}  = AbsoluteTiming[EnsembleSFF[qktEnsemble["Odd"], tauMax, nPtsSFF]];
Print["  Even: ", Round[t1, 0.1], "s, Odd: ", Round[t2, 0.1], "s"];

SavePlot[Legended[Show[
  ListLinePlot[Transpose[{sffCOE["Tau"], sffCOE["MeanK"]}], PlotStyle -> {Red, Thick}],
  ListLinePlot[Transpose[{sffEven["Tau"], sffEven["MeanK"]}], PlotStyle -> {Blue, Thick}],
  ListLinePlot[Transpose[{sffOdd["Tau"], sffOdd["MeanK"]}], PlotStyle -> {Darker[Green], Thick}],
  Plot[COESFFAnalytical[tau], {tau, 0.001, tauMax}, PlotStyle -> {Black, Thick, Dashed}],
  $pOpts, FrameLabel -> {"\[Tau]", "K(\[Tau])"},
  PlotLabel -> "SFF: J=" <> ToString[J], PlotRange -> {{0, tauMax}, {0, 1.5}}],
  Placed[LineLegend[{Directive[Blue, Thick], Directive[Darker[Green], Thick],
    Directive[Red, Thick], Directive[Black, Thick, Dashed]},
    {"Even", "Odd", "COE", "COE analyt."}, LabelStyle -> $legStyle], {Right, Bottom}]
], "sff"];


(* ================================================================ *)
(* MODULE 7: FLUCTUATION + LENGTH SPECTRUM                          *)
(* ================================================================ *)

Print["\n=== FLUCTUATION ANALYSIS ==="];
flucEven = GetFluctuations[qktEnsemble["Even"][[1]]];
flucOdd  = GetFluctuations[qktEnsemble["Odd"][[1]]];

SavePlot[Legended[Show[
  ListLinePlot[Transpose[{flucEven["x"], flucEven["delta"]}], PlotStyle -> {Blue, AbsoluteThickness[0.5]}],
  ListLinePlot[Transpose[{flucOdd["x"], flucOdd["delta"]}], PlotStyle -> {Darker[Green], AbsoluteThickness[0.5]}],
  $pOpts, FrameLabel -> {"x", "\!\(\*OverscriptBox[\(N\), \(~\)]\)(x)"},
  PlotLabel -> "Fluctuating Staircase",
  GridLines -> {None, {0}}, GridLinesStyle -> Directive[Red, Dashed]],
  Placed[LineLegend[{Directive[Blue], Directive[Darker[Green]]},
    {"Even", "Odd"}, LabelStyle -> $legStyle], {Right, Bottom}]
], "fluctuations"];

lsEven = LengthSpectrum[flucEven]; lsOdd = LengthSpectrum[flucOdd];
SavePlot[Legended[Show[
  ListLinePlot[Transpose[{lsEven["Period"][[1 ;; 50]], lsEven["Amplitude"][[1 ;; 50]]}],
    PlotStyle -> {Blue, AbsoluteThickness[1.5]}],
  ListLinePlot[Transpose[{lsOdd["Period"][[1 ;; 50]], lsOdd["Amplitude"][[1 ;; 50]]}],
    PlotStyle -> {Red, AbsoluteThickness[1.5]}],
  $pOpts, FrameLabel -> {"Period (kicks)", "Length Spectrum"}, PlotLabel -> "Length Spectrum",
  GridLines -> {None, {0}}, GridLinesStyle -> Directive[Gray, Dashed]],
  Placed[LineLegend[{Directive[Blue], Directive[Red]},
    {"Even", "Odd"}, LabelStyle -> $legStyle], {Right, Bottom}]
], "length_spectrum"];


(* ================================================================ *)
(* MODULE 8: EIGENVECTOR ANALYSIS                                   *)
(* ================================================================ *)

Print["\n=== EIGENVECTOR ANALYSIS ==="];
{tEig, eigData} = AbsoluteTiming[GetQKTParityEigensystem[J, alpha, kCenter]];
Print["  Eigensystem computed in ", Round[tEig, 0.1], "s"];

iprsEven = Map[IPR, eigData["Even"]["Vectors"]];
iprsOdd  = Map[IPR, eigData["Odd"]["Vectors"]];
Print["  Even <IPR>=", NumberForm[Mean[iprsEven], 4], " (COE: ", NumberForm[3.0 / dimE, 4], ")"];
Print["  Odd  <IPR>=", NumberForm[Mean[iprsOdd], 4], " (COE: ", NumberForm[3.0 / dimO, 4], ")"];

SavePlot[Legended[Show[
  Histogram[{iprsEven, iprsOdd}, Automatic, "ProbabilityDensity",
    ChartStyle -> {Directive[Blue, Opacity[0.5]], Directive[Darker[Green], Opacity[0.5]]}],
  Graphics[{Red, Dashed, AbsoluteThickness[2], Line[{{3.0 / dimE, 0}, {3.0 / dimE, 100}}]}],
  $pOpts, FrameLabel -> {"IPR", "Density"}, PlotLabel -> "IPR Distribution"],
  Placed[LineLegend[{Directive[Blue, Opacity[0.5]], Directive[Darker[Green], Opacity[0.5]], Directive[Red, Dashed]},
    {"Even", "Odd", "COE: 3/N"}, LabelStyle -> $legStyle], {Right, Top}]
], "ipr"];

Export[FileNameJoin[{evenDir, "ipr.wdx"}], iprsEven, "WDX"];
Export[FileNameJoin[{oddDir, "ipr.wdx"}], iprsOdd, "WDX"];


(* ================================================================ *)
(* MODULE 9: GoF ERRORS (parallelized per realization)              *)
(* ================================================================ *)

Print["\n=== GoF PER REALIZATION ==="];
DistributeDefinitions[ComputeGoFSingle, NNSpacings, IntWignerCOE, IntPrCOEFast,
  KSFromSorted, CvMFromSorted, MSDFromSorted];

{tE, gofEvenAll} = AbsoluteTiming[
  ParallelTable[ComputeGoFSingle[qktEnsemble["Even"][[i]]],
    {i, nReal}, Method -> "CoarsestGrained"]];
{tO, gofOddAll} = AbsoluteTiming[
  ParallelTable[ComputeGoFSingle[qktEnsemble["Odd"][[i]]],
    {i, nReal}, Method -> "CoarsestGrained"]];
Print["  Even: ", Round[tE, 0.1], "s, Odd: ", Round[tO, 0.1], "s"];

labels = {"KS_s", "CvM_s", "MSD_s", "KS_r", "CvM_r", "MSD_r"};
Do[
  Print["  ", labels[[idx]], ": Even mean=",
    NumberForm[Mean[gofEvenAll[[All, idx]]], 4],
    ", Odd mean=", NumberForm[Mean[gofOddAll[[All, idx]]], 4]];
, {idx, 6}];

(* Histograms *)
Do[
  Module[{eV = gofEvenAll[[All, idx]], oV = gofOddAll[[All, idx]]},
    SavePlot[Legended[Show[
      Histogram[{eV, oV}, Automatic, "ProbabilityDensity",
        ChartStyle -> {Directive[Blue, Opacity[0.5]], Directive[Darker[Green], Opacity[0.5]]}],
      $pOpts, FrameLabel -> {labels[[idx]], "Density"},
      PlotLabel -> labels[[idx]] <> " distribution (" <> ToString[nReal] <> " realizations)"],
      Placed[LineLegend[{Directive[Blue, Opacity[0.5]], Directive[Darker[Green], Opacity[0.5]]},
        {"Even", "Odd"}, LabelStyle -> $legStyle], {Right, Top}]
    ], "gof_" <> labels[[idx]]]];
, {idx, 6}];

Export[FileNameJoin[{evenDir, "gof.wdx"}], gofEvenAll, "WDX"];
Export[FileNameJoin[{oddDir, "gof.wdx"}], gofOddAll, "WDX"];


(* ================================================================ *)
(* MODULE 10: SUMMARY                                               *)
(* ================================================================ *)

Print["\n=== SUMMARY ==="];
Print["J=", J, ", \[Alpha]=", alpha, ", k=", kCenter, " \[PlusMinus] ", delta];
Print["Ensemble: ", nReal, " QKT + ", nCOE, " COE"];
Print["Even N=", dimE, ", Odd N=", dimO];
Print["Results: ", resultsDir];
Print["=== COMPLETE ==="];
