(* ================================================================ *)
(* FILE 2: GOODNESS-OF-FIT ERROR COMPUTATION (PARALLELIZED)         *)
(* Requires: QKT_full.wl loaded, qktEnsemble defined               *)
(* Outputs: KS, CvM, MSD for I(s) and I(r), per realization + ens  *)
(* ================================================================ *)


(* ================================================================ *)
(* DEFINITIONS                                                      *)
(* ================================================================ *)

ClearAll[IntPrCOETable, IntPrCOEFast,
  KSFromSorted, CvMFromSorted, MSDFromSorted,
  ComputeGoFSingle, ComputeGoFEnsemble];

(* --- Build fast COE ratio CDF interpolator --- *)
IntPrCOETable = Table[
  {r, NIntegrate[PrCOE[x], {x, 0, r}]}, {r, 0, 1, 0.001}];
IntPrCOEFast = Interpolation[IntPrCOETable];
Print["IntPrCOEFast[0.5] = ", IntPrCOEFast[0.5]];

(* --- KS statistic from sorted data and CDF values --- *)
KSFromSorted[empCDF_List, theoCDFvals_List] :=
  Max[Abs[empCDF - theoCDFvals]];

(* --- Cramér–von Mises from sorted data --- *)
(* W^2 = (1/12n) + sum_i [(2i-1)/(2n) - F(x_i)]^2 / n *)
CvMFromSorted[n_Integer, theoCDFvals_List] :=
  Module[{midEmp},
    midEmp = N[(2 Range[n] - 1) / (2 n)];
    Total[(midEmp - theoCDFvals)^2] / n + 1 / (12 n)
  ];

(* --- MSD: mean of (F_emp - F_theo)^2 over data points --- *)
MSDFromSorted[empCDF_List, theoCDFvals_List] :=
  Mean[(empCDF - theoCDFvals)^2];

(* --- Compute all GoF metrics for one spectrum (pure numerical) --- *)
(* Returns: {KS_s, CvM_s, MSD_s, KS_r, CvM_r, MSD_r} *)
ComputeGoFSingle[specData_Association] :=
  Module[{spacings, sp, ratios, sortS, sortR, nS, nR,
          empS, empR, theoS, theoR,
          ksS, cvmS, msdS, ksR, cvmR, msdR},
    (* Spacings *)
    spacings = NNSpacings[specData];
    sortS = Sort[spacings];
    nS = Length[sortS];
    empS = N[Range[nS] / nS];
    theoS = N[IntWignerCOE /@ sortS];

    (* Ratios *)
    sp = Differences[specData["Spectrum"]];
    ratios = Min /@ Transpose[{Most[sp] / Rest[sp], Rest[sp] / Most[sp]}];
    sortR = Sort[ratios];
    nR = Length[sortR];
    empR = N[Range[nR] / nR];
    theoR = N[IntPrCOEFast /@ sortR];

    (* Metrics *)
    ksS  = KSFromSorted[empS, theoS];
    cvmS = CvMFromSorted[nS, theoS];
    msdS = MSDFromSorted[empS, theoS];

    ksR  = KSFromSorted[empR, theoR];
    cvmR = CvMFromSorted[nR, theoR];
    msdR = MSDFromSorted[empR, theoR];

    {ksS, cvmS, msdS, ksR, cvmR, msdR}
  ];

(* --- Compute GoF for pooled ensemble --- *)
ComputeGoFEnsemble[spectraList_List] :=
  Module[{allSpacings, allRatios, sortS, sortR, nS, nR,
          empS, empR, theoS, theoR,
          ksS, cvmS, msdS, ksR, cvmR, msdR, sp},
    allSpacings = Sort[Flatten[Map[NNSpacings, spectraList]]];
    sp = Flatten[Map[
      Function[{specData},
        Module[{d = Differences[specData["Spectrum"]]},
          Min /@ Transpose[{Most[d] / Rest[d], Rest[d] / Most[d]}]
        ]
      ], spectraList]];
    allRatios = Sort[sp];

    sortS = allSpacings; nS = Length[sortS];
    sortR = allRatios;   nR = Length[sortR];

    empS = N[Range[nS] / nS];
    empR = N[Range[nR] / nR];
    theoS = N[IntWignerCOE /@ sortS];
    theoR = N[IntPrCOEFast /@ sortR];

    ksS  = KSFromSorted[empS, theoS];
    cvmS = CvMFromSorted[nS, theoS];
    msdS = MSDFromSorted[empS, theoS];

    ksR  = KSFromSorted[empR, theoR];
    cvmR = CvMFromSorted[nR, theoR];
    msdR = MSDFromSorted[empR, theoR];

    <|"KS_s" -> ksS, "CvM_s" -> cvmS, "MSD_s" -> msdS,
      "KS_r" -> ksR, "CvM_r" -> cvmR, "MSD_r" -> msdR,
      "N_spacings" -> nS, "N_ratios" -> nR|>
  ];


(* ================================================================ *)
(* PARAMETERS                                                       *)
(* ================================================================ *)

nReal = Length[qktEnsemble["Even"]];
Print["Ensemble size: ", nReal];


(* ================================================================ *)
(* EXECUTION: PER-REALIZATION ERRORS (PARALLELIZED)                 *)
(* ================================================================ *)

Print["\n=== PER-REALIZATION GoF (parallel) ==="];

DistributeDefinitions[
  ComputeGoFSingle, NNSpacings, IntWignerCOE, IntPrCOEFast,
  KSFromSorted, CvMFromSorted, MSDFromSorted];

(* Even sector *)
Print["  Computing Even sector..."];
{tE, gofEvenAll} = AbsoluteTiming[
  ParallelTable[
    ComputeGoFSingle[qktEnsemble["Even"][[i]]],
    {i, nReal}, Method -> "CoarsestGrained"
  ]
];
Print["  Even done in ", Round[tE, 0.1], "s"];

(* Odd sector *)
Print["  Computing Odd sector..."];
{tO, gofOddAll} = AbsoluteTiming[
  ParallelTable[
    ComputeGoFSingle[qktEnsemble["Odd"][[i]]],
    {i, nReal}, Method -> "CoarsestGrained"
  ]
];
Print["  Odd done in ", Round[tO, 0.1], "s"];

(* Unpack into named columns *)
(* Each row: {KS_s, CvM_s, MSD_s, KS_r, CvM_r, MSD_r} *)
gofEvenTable = Association[
  "KS_s"  -> gofEvenAll[[All, 1]],
  "CvM_s" -> gofEvenAll[[All, 2]],
  "MSD_s" -> gofEvenAll[[All, 3]],
  "KS_r"  -> gofEvenAll[[All, 4]],
  "CvM_r" -> gofEvenAll[[All, 5]],
  "MSD_r" -> gofEvenAll[[All, 6]]
];

gofOddTable = Association[
  "KS_s"  -> gofOddAll[[All, 1]],
  "CvM_s" -> gofOddAll[[All, 2]],
  "MSD_s" -> gofOddAll[[All, 3]],
  "KS_r"  -> gofOddAll[[All, 4]],
  "CvM_r" -> gofOddAll[[All, 5]],
  "MSD_r" -> gofOddAll[[All, 6]]
];

(* Print summary statistics *)
Print["\n--- Even sector (per-realization) ---"];
Do[
  Print["  ", metric, ":  mean=",
    NumberForm[Mean[gofEvenTable[metric]], 4],
    "  std=", NumberForm[StandardDeviation[gofEvenTable[metric]], 4],
    "  max=", NumberForm[Max[gofEvenTable[metric]], 4]];
, {metric, {"KS_s", "CvM_s", "MSD_s", "KS_r", "CvM_r", "MSD_r"}}];

Print["\n--- Odd sector (per-realization) ---"];
Do[
  Print["  ", metric, ":  mean=",
    NumberForm[Mean[gofOddTable[metric]], 4],
    "  std=", NumberForm[StandardDeviation[gofOddTable[metric]], 4],
    "  max=", NumberForm[Max[gofOddTable[metric]], 4]];
, {metric, {"KS_s", "CvM_s", "MSD_s", "KS_r", "CvM_r", "MSD_r"}}];

(* Export per-realization data *)
Export[FileNameJoin[{evenDir, "gof_per_realization.wdx"}], gofEvenTable, "WDX"];
Export[FileNameJoin[{oddDir, "gof_per_realization.wdx"}], gofOddTable, "WDX"];
Print["\nPer-realization GoF exported"];


(* ================================================================ *)
(* EXECUTION: ENSEMBLE-POOLED ERRORS                                *)
(* ================================================================ *)

Print["\n=== ENSEMBLE-POOLED GoF ==="];

{tEE, gofEvenEns} = AbsoluteTiming[
  ComputeGoFEnsemble[qktEnsemble["Even"]]];
{tOE, gofOddEns} = AbsoluteTiming[
  ComputeGoFEnsemble[qktEnsemble["Odd"]]];

Print["\n--- Even ensemble ---"];
Do[Print["  ", k, " = ", NumberForm[gofEvenEns[k], 5]],
  {k, {"KS_s", "CvM_s", "MSD_s", "KS_r", "CvM_r", "MSD_r"}}];

Print["\n--- Odd ensemble ---"];
Do[Print["  ", k, " = ", NumberForm[gofOddEns[k], 5]],
  {k, {"KS_s", "CvM_s", "MSD_s", "KS_r", "CvM_r", "MSD_r"}}];

Export[FileNameJoin[{evenDir, "gof_ensemble.wdx"}], gofEvenEns, "WDX"];
Export[FileNameJoin[{oddDir, "gof_ensemble.wdx"}], gofOddEns, "WDX"];


(* ================================================================ *)
(* VISUALIZATION: HISTOGRAMS OF PER-REALIZATION ERRORS              *)
(* ================================================================ *)

Print["\n=== GoF HISTOGRAMS ==="];

$pOpts = Sequence[ImageSize -> 800, Frame -> True,
  FrameStyle -> Directive[30, Black],
  LabelStyle -> Directive[30, Black],
  Background -> White];

Do[
  Module[{eVals, oVals, plt},
    eVals = gofEvenTable[metric];
    oVals = gofOddTable[metric];
    plt = Legended[
      Show[
        Histogram[{eVals, oVals}, Automatic, "ProbabilityDensity",
          ChartStyle -> {Directive[Blue, Opacity[0.5]],
                         Directive[Darker[Green], Opacity[0.5]]}],
        $pOpts,
        FrameLabel -> {metric, "Density"},
        PlotLabel -> metric <> " distribution across " <>
          ToString[nReal] <> " realizations"
      ],
      Placed[LineLegend[
        {Directive[Blue, Opacity[0.5]], Directive[Darker[Green], Opacity[0.5]]},
        {"Even", "Odd"},
        LabelStyle -> Directive[Black, 18]
      ], {Right, Top}]
    ];
    Export[
      FileNameJoin[{plotDir, "gof_hist_" <> metric <> ".png"}],
      Rasterize[plt, ImageResolution -> 150, Background -> White]
    ];
  ];
, {metric, {"KS_s", "CvM_s", "MSD_s", "KS_r", "CvM_r", "MSD_r"}}];

Print["GoF histograms exported"];
Print["\n=== COMPLETE ==="];
