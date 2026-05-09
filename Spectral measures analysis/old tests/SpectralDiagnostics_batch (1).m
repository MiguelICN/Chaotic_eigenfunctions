(* ================================================================ *)
(* FULL SPECTRAL DIAGNOSTICS: MULTI-J, MULTI-k SCAN                *)
(* COE reference + QKT parity sectors with ensemble averaging       *)
(* ================================================================ *)


(* === SETUP === *)
baseDir = NotebookDirectory[];
Get[FileNameJoin[{baseDir, "QKT_full.wl"}]];
LaunchKernels[];


(* === PARAMETERS === *)
jValues   = {500, 1000, 1500, 2000};
alpha     = 0.84;
kValues   = Range[10., 30., 0.5];
nQKT      = 500;           (* QKT ensemble size per k value *)
delta     = 0.01;          (* k perturbation half-width *)
nCOE      = 500;           (* COE ensemble size per J *)
Lmax      = 200;
lDomain   = Range[0.5, Lmax, 0.5];


(* === HELPERS === *)

kLabel[k_] := Module[{s = ToString[N[k]]},
  If[StringEndsQ[s, "."], s <> "0", s]];

jLabel[j_] := "J" <> StringPadLeft[ToString[j], 4, "0"];

EnsureDir[dir_] := If[!DirectoryQ[dir],
  CreateDirectory[dir, CreateIntermediateDirectories -> True]];

PlotTitle[sector_, j_, k_] := If[k == 0,
  sector <> " reference (J=" <> ToString[j] <> ")",
  sector <> ": J=" <> ToString[j] <> ", k=" <> kLabel[k]];

$pOpts = {ImageSize -> 800, Frame -> True,
  FrameStyle -> Directive[30, Black],
  LabelStyle -> Directive[30, Black]};


(* === PLOT FUNCTIONS === *)
(* All QKT data now uses ensemble keys: Sigma2TheorMean, SETheorMean *)

MakeSigma2Plot[s2Even_, s2Odd_, s2COE_, j_, k_, dimE_, dimO_] :=
  Legended[
    Show[
      (* COE error band *)
      ListLinePlot[
        {Transpose[{s2COE["LValues"],
          s2COE["Sigma2TheorMean"] + s2COE["SETheorMean"]}],
         Transpose[{s2COE["LValues"],
          s2COE["Sigma2TheorMean"] - s2COE["SETheorMean"]}]},
        PlotStyle -> None, Filling -> {1 -> {2}},
        FillingStyle -> Directive[Opacity[0.15], Red]],
      ListLinePlot[
        Transpose[{s2COE["LValues"], s2COE["Sigma2TheorMean"]}],
        PlotStyle -> {Red, Thick}],
      (* Even error band *)
      ListLinePlot[
        {Transpose[{s2Even["LValues"],
          s2Even["Sigma2TheorMean"] + s2Even["SETheorMean"]}],
         Transpose[{s2Even["LValues"],
          s2Even["Sigma2TheorMean"] - s2Even["SETheorMean"]}]},
        PlotStyle -> None, Filling -> {1 -> {2}},
        FillingStyle -> Directive[Opacity[0.15], Blue]],
      ListLinePlot[
        Transpose[{s2Even["LValues"], s2Even["Sigma2TheorMean"]}],
        PlotStyle -> {Blue, Thick}],
      (* Odd error band *)
      ListLinePlot[
        {Transpose[{s2Odd["LValues"],
          s2Odd["Sigma2TheorMean"] + s2Odd["SETheorMean"]}],
         Transpose[{s2Odd["LValues"],
          s2Odd["Sigma2TheorMean"] - s2Odd["SETheorMean"]}]},
        PlotStyle -> None, Filling -> {1 -> {2}},
        FillingStyle -> Directive[Opacity[0.15], Darker[Green]]],
      ListLinePlot[
        Transpose[{s2Odd["LValues"], s2Odd["Sigma2TheorMean"]}],
        PlotStyle -> {Darker[Green], Thick}],
      (* Asymptotic + Poisson *)
      Plot[COEAsymptoticSigma2[L], {L, 0.5, Lmax},
        PlotStyle -> {Orange, Thick, Dashed}],
      Plot[L, {L, 0.5, Lmax},
        PlotStyle -> {Black, Dashed, Thin}],
      Evaluate@$pOpts,
      FrameLabel -> {"L",
        "\!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\)(L)"},
      PlotLabel -> PlotTitle["\!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\)", j, k],
      PlotRange -> {All, {0, Automatic}}
    ],
    Placed[LineLegend[
      {Directive[Blue, Thick], Directive[Darker[Green], Thick],
       Directive[Red, Thick], Directive[Orange, Thick, Dashed],
       Directive[Black, Dashed, Thin]},
      {"Even (N=" <> ToString[dimE] <> ")",
       "Odd (N=" <> ToString[dimO] <> ")",
       "COE(" <> ToString[nCOE] <> ")",
       "COE asymptotic", "Poisson"},
      LabelStyle -> Directive[Black, 18]
    ], {Right, Bottom}]
  ];

MakeResidualPlot[s2Even_, s2Odd_, s2COE_, j_, k_] :=
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
    Legended[
      Show[
        (* Even error band *)
        ListLinePlot[
          {Transpose[{lR, resE + seE}], Transpose[{lR, resE - seE}]},
          PlotStyle -> None, Filling -> {1 -> {2}},
          FillingStyle -> Directive[Opacity[0.2], Blue]],
        ListLinePlot[Transpose[{lR, resE}], PlotStyle -> {Blue, Thick}],
        (* Odd error band *)
        ListLinePlot[
          {Transpose[{lR, resO + seO}], Transpose[{lR, resO - seO}]},
          PlotStyle -> None, Filling -> {1 -> {2}},
          FillingStyle -> Directive[Opacity[0.2], Darker[Green]]],
        ListLinePlot[Transpose[{lR, resO}],
          PlotStyle -> {Darker[Green], Thick}],
        Plot[0, {x, 0, Max[lR]}, PlotStyle -> {Red, Dashed}],
        Evaluate@$pOpts,
        FrameLabel -> {"L",
          "\!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\)(QKT) - \!\(\*SuperscriptBox[\(\[Sigma]\), \(2\)]\)(COE)"},
        PlotLabel -> PlotTitle["Residual", j, k]
      ],
      Placed[LineLegend[
        {Directive[Blue, Thick], Directive[Darker[Green], Thick]},
        {"Even - COE", "Odd - COE"},
        LabelStyle -> Directive[Black, 18]
      ], {Right, Bottom}]
    ]
  ];

MakeDualEstimatorPlot[s2_, sector_, j_, k_] :=
  Module[{rd},
    rd = Abs[s2["Sigma2TheorMean"] - s2["Sigma2SampleMean"]] /
      (s2["Sigma2TheorMean"] + 10^-10);
    ListLinePlot[Transpose[{s2["LValues"], rd}],
      Filling -> Axis, FillingStyle -> Directive[Opacity[0.3], Blue],
      Evaluate@$pOpts, PlotRange -> {0, All},
      FrameLabel -> {"L", "Relative difference"},
      PlotLabel -> PlotTitle[sector <> " dual estimator", j, k]]
  ];

MakeNNSDPlot[spacings_, sector_, j_, k_] :=
  Legended[
    Show[
      Histogram[spacings, {0, 4, 0.05}, "ProbabilityDensity",
        ChartStyle -> Directive[LightBlue, EdgeForm[Gray]],
        PlotRange -> {{0, 3.5}, {0, 1.1}}],
      Plot[WignerCOE[s], {s, 0, 3.5}, PlotStyle -> {Red, Thick}],
      Plot[Exp[-s], {s, 0, 3.5}, PlotStyle -> {Black, Dashed}],
      Evaluate@$pOpts,
      FrameLabel -> {"s", "P(s)"},
      PlotLabel -> PlotTitle[sector, j, k]
    ],
    Placed[LineLegend[
      {Directive[Red, Thick], Directive[Black, Dashed]},
      {"COE surmise", "Poisson"},
      LabelStyle -> Directive[Black, 18]
    ], {Right, Bottom}]
  ];

MakeCumulativePlot[spacings_, sector_, j_, k_] :=
  Module[{sorted, nSp, empCDF},
    sorted = Sort[spacings];
    nSp = Length[sorted];
    empCDF = Range[nSp] / nSp;
    Legended[
      Show[
        ListLinePlot[Transpose[{sorted, empCDF}],
          PlotStyle -> {Blue, Thick},
          PlotRange -> {{0, 3.5}, {0, 1}}],
        Plot[IntWignerCOE[s], {s, 0, 3.5}, PlotStyle -> {Red, Thick}],
        Plot[1 - Exp[-s], {s, 0, 3.5}, PlotStyle -> {Black, Dashed}],
        Evaluate@$pOpts,
        FrameLabel -> {"s", "I(s)"},
        PlotLabel -> PlotTitle[sector, j, k]
      ],
      Placed[LineLegend[
        {Directive[Blue, Thick], Directive[Red, Thick],
         Directive[Black, Dashed]},
        {"Data", "COE surmise", "Poisson"},
        LabelStyle -> Directive[Black, 18]
      ], {Right, Bottom}]
    ]
  ];


(* ================================================================ *)
(* MAIN LOOP                                                        *)
(* ================================================================ *)

Do[
  Print["\n", StringJoin[ConstantArray["=", 60]]];
  Print["  J = ", j, "  (", jLabel[j], ")"];
  Print[StringJoin[ConstantArray["=", 60]]];

  (* --- Directories --- *)
  jDir    = FileNameJoin[{baseDir, "Results", jLabel[j]}];
  evenDir = FileNameJoin[{jDir, "Even"}];
  oddDir  = FileNameJoin[{jDir, "Odd"}];
  plotDir = FileNameJoin[{jDir, "plots"}];
  EnsureDir[evenDir]; EnsureDir[oddDir]; EnsureDir[plotDir];

  dimCOE = j + 1;

  (* ============================================================ *)
  (* A. COE REFERENCE (one ensemble per J)                         *)
  (* ============================================================ *)

  Print["  Generating ", nCOE, " COE(", dimCOE, ") matrices..."];
  {tCOE, coeSpectra} = AbsoluteTiming[GenerateCOEEnsemble[dimCOE, nCOE]];
  Print["  COE done in ", Round[tCOE, 0.1], "s"];

  Print["  Computing COE Sigma2..."];
  {tS2, sigma2COE} = AbsoluteTiming[EnsembleSigma2[coeSpectra, lDomain]];
  Print["  COE Sigma2 done in ", Round[tS2, 0.1], "s"];

  (* COE NNSD reference *)
  allSpacingsCOE = Sort[Flatten[Map[NNSpacings, coeSpectra]]];

  Export[FileNameJoin[{plotDir, "nnsd_coe.pdf"}],
    MakeNNSDPlot[allSpacingsCOE, "COE", j, 0]];
  Export[FileNameJoin[{plotDir, "cumulative_coe.pdf"}],
    MakeCumulativePlot[allSpacingsCOE, "COE", j, 0]];

  coeSpectra = .; allSpacingsCOE = .;

  (* ============================================================ *)
  (* B. QKT: ensemble at each k value                              *)
  (* ============================================================ *)

  Do[
    kStr = kLabel[k];
    Print["  k=", kStr, ": generating ", nQKT, " QKT matrices..."];

    (* Generate ensemble: nQKT matrices with k in [k-delta, k+delta] *)
    SeedRandom[Hash[{j, k, 42}]];
    kList = RandomReal[{k - delta, k + delta}, nQKT];

    {tDiag, qktResult} = AbsoluteTiming[
      GenerateQKTEnsemble[j, alpha, kList]
    ];

    evenSpectra = qktResult["Even"];
    oddSpectra  = qktResult["Odd"];
    dimE = evenSpectra[[1]]["N"];
    dimO = oddSpectra[[1]]["N"];

    Print["    diag: ", Round[tDiag, 0.1], "s  [Even N=", dimE,
      ", Odd N=", dimO, "]"];

    (* --- Save all spectra for later reanalysis --- *)
    Export[FileNameJoin[{evenDir, "spectra_k" <> kStr <> ".wdx"}],
      evenSpectra, "WDX"];
    Export[FileNameJoin[{oddDir, "spectra_k" <> kStr <> ".wdx"}],
      oddSpectra, "WDX"];

    (* --- Spacings: pool across ensemble --- *)
    spacingsEven = Sort[Flatten[Map[NNSpacings, evenSpectra]]];
    spacingsOdd  = Sort[Flatten[Map[NNSpacings, oddSpectra]]];

    Export[FileNameJoin[{evenDir, "spacings_k" <> kStr <> ".wdx"}],
      spacingsEven, "WDX"];
    Export[FileNameJoin[{oddDir, "spacings_k" <> kStr <> ".wdx"}],
      spacingsOdd, "WDX"];

    (* --- Sigma2: ensemble average --- *)
    {tS2E, s2Even} = AbsoluteTiming[EnsembleSigma2[evenSpectra, lDomain]];
    {tS2O, s2Odd}  = AbsoluteTiming[EnsembleSigma2[oddSpectra, lDomain]];
    Print["    sigma2: ", Round[tS2E + tS2O, 0.1], "s"];

    Export[FileNameJoin[{evenDir, "sigma2_k" <> kStr <> ".wdx"}],
      s2Even, "WDX"];
    Export[FileNameJoin[{oddDir, "sigma2_k" <> kStr <> ".wdx"}],
      s2Odd, "WDX"];

    (* --- Plots --- *)
    Export[FileNameJoin[{plotDir, "sigma2_k" <> kStr <> ".pdf"}],
      MakeSigma2Plot[s2Even, s2Odd, sigma2COE, j, k, dimE, dimO]];
    Export[FileNameJoin[{plotDir, "residual_k" <> kStr <> ".pdf"}],
      MakeResidualPlot[s2Even, s2Odd, sigma2COE, j, k]];
    Export[FileNameJoin[{plotDir, "dual_even_k" <> kStr <> ".pdf"}],
      MakeDualEstimatorPlot[s2Even, "Even", j, k]];
    Export[FileNameJoin[{plotDir, "dual_odd_k" <> kStr <> ".pdf"}],
      MakeDualEstimatorPlot[s2Odd, "Odd", j, k]];
    Export[FileNameJoin[{plotDir, "nnsd_even_k" <> kStr <> ".pdf"}],
      MakeNNSDPlot[spacingsEven, "Even", j, k]];
    Export[FileNameJoin[{plotDir, "nnsd_odd_k" <> kStr <> ".pdf"}],
      MakeNNSDPlot[spacingsOdd, "Odd", j, k]];
    Export[FileNameJoin[{plotDir, "cumulative_even_k" <> kStr <> ".pdf"}],
      MakeCumulativePlot[spacingsEven, "Even", j, k]];
    Export[FileNameJoin[{plotDir, "cumulative_odd_k" <> kStr <> ".pdf"}],
      MakeCumulativePlot[spacingsOdd, "Odd", j, k]];

    Print["    k=", kStr, " DONE"];

    (* Free per-k data *)
    qktResult = .; evenSpectra = .; oddSpectra = .;
    spacingsEven = .; spacingsOdd = .;
    s2Even = .; s2Odd = .;

  , {k, kValues}];

  (* Free COE Sigma2 for this J *)
  sigma2COE = .;

  Print["  ", jLabel[j], " COMPLETE"];

, {j, jValues}];

Print["\n=== ALL DONE ==="];
