(* ================================================================ *)
(* FILE 1: INDIVIDUAL REALIZATION COLLAGE PLOTS                     *)
(* Requires: QKT_full.wl loaded, qktEnsemble defined               *)
(* ================================================================ *)


(* ================================================================ *)
(* DEFINITIONS                                                      *)
(* ================================================================ *)

ClearAll[IntPrCOETable, IntPrCOEFast, GetFluctuationsCentered,
  ExtractSpacingsAndRatios, MakeRowPlots];

(* --- Build fast COE ratio CDF lookup (replaces $IntPrCOEInterp) --- *)
IntPrCOETable = Table[
  {r, NIntegrate[PrCOE[x], {x, 0, r}]}, {r, 0, 1, 0.001}];
IntPrCOEFast = Interpolation[IntPrCOETable];
Print["IntPrCOEFast test: ", IntPrCOEFast[0.5]];  (* should be ~0.458 *)

(* --- Fluctuations centered on zero --- *)
GetFluctuationsCentered[specData_Association] :=
  Module[{x = specData["Spectrum"], n = specData["N"], raw},
    raw = Range[n] - x;
    <|"x" -> x, "delta" -> raw - Mean[raw], "N" -> n|>
  ];

(* --- Extract spacings + ratios from one spectrum --- *)
ExtractSpacingsAndRatios[specData_Association] :=
  Module[{spacings, sp, ratios},
    spacings = NNSpacings[specData];
    sp = Differences[specData["Spectrum"]];
    ratios = Min /@ Transpose[{Most[sp] / Rest[sp], Rest[sp] / Most[sp]}];
    <|"Spacings" -> spacings, "SortedS" -> Sort[spacings],
      "Ratios" -> ratios, "SortedR" -> Sort[ratios]|>
  ];

(* --- Compact plot options --- *)
$cOpts = Sequence[Frame -> True,
  FrameStyle -> Directive[Black, 16],
  LabelStyle -> Directive[Black, 16]];
$cSize = 380;

(* --- All plots for one sector --- *)
MakeRowPlots[sortedS_List, sortedR_List, flucData_Association, label_String] :=
  Module[{nS, nR, empS, empR, theoS, theoR, ksS, ksR,
          diffS, diffR, pPDFS, pPDFR, pCDFS, pCDFR,
          pMSDS, pMSDR, pCDFRplot, pFluc},

    nS = Length[sortedS]; nR = Length[sortedR];
    empS  = N[Range[nS] / nS];
    empR  = N[Range[nR] / nR];
    theoS = N[IntWignerCOE /@ sortedS];
    theoR = N[IntPrCOEFast /@ sortedR];

    ksS = Max[Abs[empS - theoS]];
    ksR = Max[Abs[empR - theoR]];
    diffS = (empS - theoS)^2;
    diffR = (empR - theoR)^2;

    (* P(s) *)
    pPDFS = Show[
      Histogram[sortedS, {0, 4, 0.05}, "ProbabilityDensity",
        ChartStyle -> Directive[LightGray, EdgeForm[Gray]],
        PlotRange -> {{0, 3.5}, {0, 1.1}}],
      Plot[WignerCOE[x], {x, 0, 3.5}, PlotStyle -> {Red, Dashed, Thick}],
      Plot[Exp[-x], {x, 0, 3.5}, PlotStyle -> {Black, Dashed}],
      $cOpts, FrameLabel -> {"s", "P(s)"},
      PlotLabel -> label <> " P(s)", ImageSize -> $cSize];

    (* P(r) *)
    pPDFR = Show[
      Histogram[sortedR, {0, 1, 0.02}, "ProbabilityDensity",
        ChartStyle -> Directive[LightGray, EdgeForm[Gray]],
        PlotRange -> {{0, 1}, {0, 2.2}}],
      Plot[PrCOE[x], {x, 0, 1}, PlotStyle -> {Red, Dashed, Thick}],
      Plot[PrPoisson[x], {x, 0, 1}, PlotStyle -> {Black, Dashed}],
      $cOpts, FrameLabel -> {"r", "P(r)"},
      PlotLabel -> label <> " P(r)", ImageSize -> $cSize];

    (* I(s) *)
    pCDFS = Show[
      ListLinePlot[Transpose[{sortedS, empS}],
        PlotStyle -> {Blue, Thick}, PlotRange -> {{0, 3.5}, {0, 1}}],
      Plot[IntWignerCOE[x], {x, 0, 3.5}, PlotStyle -> {Red, Dashed}],
      Plot[1 - Exp[-x], {x, 0, 3.5}, PlotStyle -> {Black, Dashed}],
      $cOpts, FrameLabel -> {"s", "I(s)"},
      PlotLabel -> label <> " I(s)  KS=" <> ToString[NumberForm[ksS, 4]],
      ImageSize -> $cSize];

    (* I(r) *)
    pCDFRplot = Show[
      ListLinePlot[Transpose[{sortedR, empR}],
        PlotStyle -> {Blue, Thick}, PlotRange -> {{0, 1}, {0, 1}}],
      ListLinePlot[Transpose[{sortedR, theoR}],
        PlotStyle -> {Red, Dashed}],
      $cOpts, FrameLabel -> {"r", "I(r)"},
      PlotLabel -> label <> " I(r)  KS=" <> ToString[NumberForm[ksR, 4]],
      ImageSize -> $cSize];

    (* CDF diff (s) *)
    pMSDS = ListLinePlot[Transpose[{sortedS, diffS}],
      PlotStyle -> {Blue, Thick}, Filling -> Axis,
      FillingStyle -> Directive[Opacity[0.3], Blue],
      $cOpts, FrameLabel -> {"s", "(I-I_COE)\[Superscript]2"},
      PlotLabel -> label <> " CDF diff (s)",
      PlotRange -> All, ImageSize -> $cSize];

    (* CDF diff (r) *)
    pMSDR = ListLinePlot[Transpose[{sortedR, diffR}],
      PlotStyle -> {Blue, Thick}, Filling -> Axis,
      FillingStyle -> Directive[Opacity[0.3], Blue],
      $cOpts, FrameLabel -> {"r", "(I-I_COE)\[Superscript]2"},
      PlotLabel -> label <> " CDF diff (r)",
      PlotRange -> All, ImageSize -> $cSize];

    (* Fluctuation *)
    pFluc = ListLinePlot[Transpose[{flucData["x"], flucData["delta"]}],
      PlotStyle -> {Blue, AbsoluteThickness[0.5]},
      $cOpts,
      FrameLabel -> {"x", "\!\(\*OverscriptBox[\(N\), \(~\)]\)(x)"},
      PlotLabel -> label <> " Fluctuation",
      GridLines -> {None, {0}},
      GridLinesStyle -> Directive[Red, Dashed],
      ImageSize -> $cSize];

    (* Return: P(s), P(r), I(s), I(r), diff(s), diff(r), fluc *)
    {pPDFS, pPDFR, pCDFS, pCDFRplot, pMSDS, pMSDR, pFluc}
  ];


(* ================================================================ *)
(* PARAMETERS                                                       *)
(* ================================================================ *)

collageDir = FileNameJoin[{resultsDir, "individual"}];
EnsureDir[collageDir];
nReal = Length[qktEnsemble["Even"]];


(* ================================================================ *)
(* EXECUTION: INDIVIDUAL REALIZATIONS                               *)
(* ================================================================ *)

Print["=== INDIVIDUAL COLLAGES ==="];

Do[
  Module[{stE, stO, flE, flO, colE, colO, collage},
    stE = ExtractSpacingsAndRatios[qktEnsemble["Even"][[i]]];
    stO = ExtractSpacingsAndRatios[qktEnsemble["Odd"][[i]]];
    flE = GetFluctuationsCentered[qktEnsemble["Even"][[i]]];
    flO = GetFluctuationsCentered[qktEnsemble["Odd"][[i]]];

    colE = MakeRowPlots[stE["SortedS"], stE["SortedR"], flE, "Even"];
    colO = MakeRowPlots[stO["SortedS"], stO["SortedR"], flO, "Odd"];

    collage = Grid[
      Prepend[
        Transpose[{colE, colO}],
        {Style["Even", Bold, 20], Style["Odd", Bold, 20]}
      ],
      Spacings -> {1, 0.5},
      Frame -> All, FrameStyle -> GrayLevel[0.85]
    ];

    Export[
      FileNameJoin[{collageDir,
        "realization_" <> StringPadLeft[ToString[i], 4, "0"] <> ".png"}],
      Rasterize[collage, ImageResolution -> 150, Background -> White]
    ];

    If[Mod[i, 50] == 0 || i == 1,
      Print["  Exported realization ", i, "/", nReal]];
  ];
, {i, nReal}];

Print["Individual collages complete: ", collageDir];


(* ================================================================ *)
(* EXECUTION: ENSEMBLE-POOLED COLLAGE                               *)
(* ================================================================ *)

Print["\n=== ENSEMBLE COLLAGE ==="];

Module[{spE, spO, rE, rO, flE, flO, colE, colO, ensCollage},
  spE = Sort[Flatten[Map[NNSpacings, qktEnsemble["Even"]]]];
  spO = Sort[Flatten[Map[NNSpacings, qktEnsemble["Odd"]]]];
  rE  = Sort[PoolSpacingRatios[qktEnsemble["Even"]]];
  rO  = Sort[PoolSpacingRatios[qktEnsemble["Odd"]]];

  Print["  Even: <s>=", NumberForm[Mean[spE], 5],
    ", <r>=", NumberForm[Mean[rE], 5]];
  Print["  Odd:  <s>=", NumberForm[Mean[spO], 5],
    ", <r>=", NumberForm[Mean[rO], 5]];

  flE = GetFluctuationsCentered[qktEnsemble["Even"][[1]]];
  flO = GetFluctuationsCentered[qktEnsemble["Odd"][[1]]];

  colE = MakeRowPlots[spE, rE, flE, "Even (ensemble)"];
  colO = MakeRowPlots[spO, rO, flO, "Odd (ensemble)"];

  ensCollage = Grid[
    Prepend[
      Transpose[{colE, colO}],
      {Style["Even (pooled " <> ToString[nReal] <> ")", Bold, 20],
       Style["Odd (pooled " <> ToString[nReal] <> ")", Bold, 20]}
    ],
    Spacings -> {1, 0.5},
    Frame -> All, FrameStyle -> GrayLevel[0.85]
  ];

  Export[
    FileNameJoin[{plotDir, "ensemble_collage.png"}],
    Rasterize[ensCollage, ImageResolution -> 150, Background -> White]
  ];
  Print["  Ensemble collage exported"];
];
