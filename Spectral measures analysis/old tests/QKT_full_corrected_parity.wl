(* ::Package:: *)

(* ================================================================ *)
(* QKT_full.wl                                                      *)
(* Quantum & Classical Kicked Top: Operators, Parity Decomposition,  *)
(* Spectral Statistics, and Phase-Space Tools                        *)
(*                                                                   *)
(* Convention: all spin operators use the DESCENDING basis            *)
(*   |j,j>, |j,j-1>, ..., |j,-j>   (index 1 = highest m)           *)
(*                                                                   *)
(* Requires: ForScience paclet (auto-installed if absent)            *)
(* ================================================================ *)

If[Length[PacletFind["ForScience"]] == 0,
    PacletInstall[FileNameJoin[{DirectoryName[$InputFileName], "ForScience-0.88.45.paclet"}]]
];

BeginPackage["QKT`"];
<<ForScience`;


(* ::Section:: *)
(*Usage definitions*)


(* ================================================================ *)
(* 1. SPIN OPERATORS                                                *)
(* ================================================================ *)

generateSpinOperators::usage =
  "generateSpinOperators[j] returns <|\"Sx\"->..., \"Sy\"->..., \"Sz\"->..., \"Sx2\"->...|> \
for total angular momentum j. All matrices are in the descending-m basis.";


(* ================================================================ *)
(* 2. FLOQUET OPERATORS                                             *)
(* ================================================================ *)

Floq::usage =
  "Floq[j, \[Alpha], k] returns the Floquet operator U = exp(-i k Sx^2/(2j)) \[CenterDot] exp(-i \[Alpha] Sz). \
Twist around Sx, free precession around Sz. Dimension (2j+1).";

Floqsym::usage =
  "Floqsym[j, \[Alpha], k] returns the symmetrized Floquet operator \
exp(-i \[Alpha] Sz/2) \[CenterDot] exp(-i k Sx^2/(2j)) \[CenterDot] exp(-i \[Alpha] Sz/2).";

Floqn::usage =
  "Floqn[j, \[Alpha], k, nVec] returns the Floquet operator with twist around Sx \
and rotation along an arbitrary axis nVec (3-vector, auto-normalized).";

Floqb::usage =
  "Floqb[j, \[Alpha], k, \[Gamma]] returns a Floquet operator that breaks \
one time-reversal symmetry (COE -> CUE transition).";

Floqbn::usage =
  "Floqbn[j, \[Alpha], k, nVec, \[Gamma]] returns a Floquet operator that \
can break both time-reversal symmetries (COE -> CUE).";


(* ================================================================ *)
(* 3. PARITY DECOMPOSITION                                         *)
(* ================================================================ *)

ParitySectorIndices::usage =
  "ParitySectorIndices[j] returns <|\"Even\"->idxList, \"Odd\"->idxList|> \
giving the row/column indices of the two decoupled parity sectors \
of the Floquet matrix. Works for both integer and half-integer j. \
\
The decoupling arises because Sx^2 couples m to m\[PlusMinus]2, so alternating \
array indices always decouple. For odd integer j, index 1 corresponds \
to odd m, so the Even-m sector starts at index 2. For even integer j \
and all half-integer j, index 1 starts the first sector.";

GetQKTParitySpectra::usage =
  "GetQKTParitySpectra[j, \[Alpha], k] returns <|\"Even\"->specData, \"Odd\"->specData|> \
where each specData is <|\"N\"->dim, \"Phases\"->sorted phases, \"Spectrum\"->unfolded|>. \
Uses ParitySectorIndices for correct parity separation.";

GetQKTParityEigensystem::usage =
  "GetQKTParityEigensystem[j, \[Alpha], k] returns the full eigensystem for both parity sectors: \
<|\"Even\"-><|\"Spec\"->specData, \"Vectors\"->vecs, \"FullVectors\"->reconVecs|>, \"Odd\"->...|>. \
FullVectors are reconstructed in the original (2j+1)-dimensional basis.";

GenerateQKTEnsemble::usage =
  "GenerateQKTEnsemble[j, \[Alpha], kList] generates QKT spectra for both parity sectors \
over all k values in kList (parallelized). Returns <|\"Even\"->list, \"Odd\"->list|>.";


(* ================================================================ *)
(* 4. COHERENT STATES                                               *)
(* ================================================================ *)

generateCoherentStateCompiler::usage =
  "generateCoherentStateCompiler[] returns a compiled function f[J, q, p] that \
generates the spin-j coherent state |q,p> in the stereographic disk representation. \
Uses the descending basis convention.";


(* ================================================================ *)
(* 5. SPECTRAL STATISTICS                                           *)
(* ================================================================ *)

UnfoldLinear::usage =
  "UnfoldLinear[phases] linearly unfolds a sorted list of quasienergy phases \
on [0,2\[Pi]) to mean spacing 1. Returns <|\"N\"->n, \"Phases\"->..., \"Spectrum\"->...|>.";

ValidateUnfolding::usage =
  "ValidateUnfolding[specData] computes staircase residuals and KS statistic \
to assess unfolding quality. Returns <|\"Phases\"->..., \"Residuals\"->..., \"KS\"->..., ...|>.";

NNSpacings::usage =
  "NNSpacings[specData] returns the nearest-neighbor spacings of the unfolded \
spectrum, including the wrap-around spacing on the circle.";

PoolSpacings::usage =
  "PoolSpacings[spectraList] pools and sorts all NN spacings from an ensemble of spectra.";

PoolSpacingRatios::usage =
  "PoolSpacingRatios[spectraList] pools all consecutive spacing ratios \
r_i = min(s_i,s_{i+1})/max(s_i,s_{i+1}) from an ensemble.";

MeanLevelSpacingRatio::usage =
  "MeanLevelSpacingRatio[eigenvalues] returns <r>, the mean consecutive spacing ratio.";

GetFluctuations::usage =
  "GetFluctuations[specData] returns the fluctuating part of the spectral staircase, \
centered on zero: delta_i = (i - x_i) - mean. Returns <|\"x\"->..., \"delta\"->..., \"N\"->n|>.";

LengthSpectrum::usage =
  "LengthSpectrum[fluc] returns the length spectrum (DFT of the fluctuating density). \
Periods in units of kick periods. Returns <|\"Period\"->..., \"Amplitude\"->...|>.";

FluctuationPowerSpectrum::usage =
  "FluctuationPowerSpectrum[fluc] returns |delta_k|^2, the power spectrum of \
staircase fluctuations. Returns <|\"Freq\"->..., \"Power\"->...|>.";

IPR::usage =
  "IPR[vec] returns the inverse participation ratio Sum[|v_i|^4] / (Sum[|v_i|^2])^2.";


(* ================================================================ *)
(* 6. ANALYTICAL RMT PREDICTIONS (COE, beta=1)                     *)
(* ================================================================ *)

WignerCOE::usage =
  "WignerCOE[s] returns the COE (beta=1) Wigner surmise: (Pi/2) s Exp[-Pi s^2/4].";

IntWignerCOE::usage =
  "IntWignerCOE[s] returns the integrated COE Wigner surmise: 1 - Exp[-Pi s^2/4].";

PrCOE::usage =
  "PrCOE[r] returns the COE spacing ratio distribution (Atas et al. 2013).";

PrCUE::usage =
  "PrCUE[r] returns the CUE (beta=2) spacing ratio distribution.";

PrPoisson::usage =
  "PrPoisson[r] returns the Poisson spacing ratio distribution: 2/(1+r)^2.";

COEAsymptoticSigma2::usage =
  "COEAsymptoticSigma2[L] returns the asymptotic COE number variance \
(2/Pi^2)(Ln[2 Pi L] + 1 + gamma - Pi^2/8). Valid for 1 << L << N.";

COESFFAnalytical::usage =
  "COESFFAnalytical[tau] returns the analytical COE spectral form factor (large-N limit).";


(* ================================================================ *)
(* 7. NUMBER VARIANCE                                               *)
(* ================================================================ *)

CompileLevelCounter::usage =
  "CompileLevelCounter[extSpectrum] returns a compiled binary-search function \
that counts levels <= z in the extended spectrum.";

ComputeSigma2Sliding::usage =
  "ComputeSigma2Sliding[counterFunc, nLevels, lValues] computes the number variance \
via systematic sliding windows (one center per integer). Returns {L, Sigma2_theor, Sigma2_sample} \
for the dual estimator diagnostic.";

EvaluateSpectrumSigma2::usage =
  "EvaluateSpectrumSigma2[specData, lValues] computes Sigma2 for one spectrum.";

EnsembleSigma2::usage =
  "EnsembleSigma2[spectraList, lValues] returns ensemble-averaged Sigma2 with \
dual estimator and standard errors.";

SingleSpectrumSigma2::usage =
  "SingleSpectrumSigma2[specData, lValues] returns Sigma2 from spectral averaging only.";


(* ================================================================ *)
(* 8. SPECTRAL FORM FACTOR                                         *)
(* ================================================================ *)

CompileSFF::usage =
  "CompileSFF is a compiled function: CompileSFF[spectrum, tauValues] -> K(tau).";

EnsembleSFF::usage =
  "EnsembleSFF[spectraList, tauMax, nPts] returns ensemble-averaged SFF with errors.";

SingleSFF::usage =
  "SingleSFF[specData, tauMax, nPts] computes K(tau) for a single spectrum.";


(* ================================================================ *)
(* 9. COE RANDOM MATRICES                                          *)
(* ================================================================ *)

RandomCOESpectrum::usage =
  "RandomCOESpectrum[n] generates one COE(n) matrix via U*U^T (Haar CUE from QR) \
and returns its unfolded spectrum.";

GenerateCOEEnsemble::usage =
  "GenerateCOEEnsemble[n, nReal] generates nReal COE(n) unfolded spectra (parallelized).";


(* ================================================================ *)
(* 10. CLASSICAL KICKED TOP                                         *)
(* ================================================================ *)

ClassicalToSphere::usage =
  "ClassicalToSphere[{q1,q2}] maps a point on the stereographic disk to (Sx,Sy,Sz) on S^2.";

ClassicalToDisk::usage =
  "ClassicalToDisk[{Sx,Sy,Sz}] maps a point on S^2 back to the stereographic disk.";

ClassicalStep::usage =
  "ClassicalStep[{Sx,Sy,Sz}, alpha, k, order] performs one stroboscopic step. \
order=0: rotation then kick; order=1: kick then rotation.";

ClassicalMapeo::usage =
  "ClassicalMapeo[{q1,q2}, alpha, k, n, order] generates a trajectory of n points \
on the disk starting from (q1,q2). Compiled, listable.";


(* ================================================================ *)
(* Implementations                                                  *)
(* ================================================================ *)

Begin["`Private`"];


(* ================================================================ *)
(* 1. SPIN OPERATORS (descending basis: index 1 = m=j)             *)
(* ================================================================ *)

ClearAll[sz, sp, sm, sx, sy, sx2];

sz[j_] := sz[j] = SparseArray[Band[{1, 1}] -> Range[j, -j, -1]];

sp[j_] := sp[j] = SparseArray[Band[{1, 2}] ->
    Table[Sqrt[j (j + 1) - m (m + 1)] // N, {m, j - 1, -j, -1}],
    {2 j + 1, 2 j + 1}];

sm[j_] := sm[j] = SparseArray[Band[{2, 1}] ->
    Table[Sqrt[j (j + 1) - m (m - 1)] // N, {m, j, -j + 1, -1}],
    {2 j + 1, 2 j + 1}];

sx[j_] := sx[j] = (1/2) (sp[j] + sm[j]);
sy[j_] := sy[j] = (1/(2 I)) (sp[j] - sm[j]);
sx2[j_] := sx2[j] = sx[j] . sx[j];

generateSpinOperators[j_] :=
  <|"Sx" -> sx[j], "Sy" -> sy[j], "Sz" -> sz[j], "Sx2" -> sx2[j]|>;


(* ================================================================ *)
(* 2. FLOQUET OPERATORS                                             *)
(* ================================================================ *)

twistPart[j_, k_] := MatrixExp[(-I k / (2 j)) sx2[j]];
freePart[j_, \[Alpha]_] :=
  SparseArray[Band[{1, 1}] -> Exp[-I \[Alpha] Range[j, -j, -1]]];

Floq[j_, \[Alpha]_, k_] := twistPart[j, k] . freePart[j, \[Alpha]];
Floqsym[j_, \[Alpha]_, k_] :=
  freePart[j, \[Alpha]/2] . twistPart[j, k] . freePart[j, \[Alpha]/2];

twistPartGeneral[j_, \[Alpha]_, nVec_] :=
  twistPartGeneral[j, \[Alpha], nVec] =
    Module[{u = Normalize[nVec], gen},
      gen = u . {sx[j], sy[j], sz[j]};
      MatrixExp[-I \[Alpha] gen]
    ];

Floqn[j_, \[Alpha]_, k_, nVec_] :=
  twistPart[j, k] . twistPartGeneral[j, \[Alpha], nVec];


(* ================================================================ *)
(* 3. PARITY DECOMPOSITION                                         *)
(*                                                                   *)
(* Sx^2 couples m to m +/- 2, so the Floquet matrix in the |j,m>    *)
(* basis decouples into two blocks: even-indexed and odd-indexed     *)
(* rows/columns. Which physical sector (even-m or odd-m) maps to     *)
(* which index set depends on whether j is odd-integer, even-integer *)
(* or half-integer.                                                  *)
(* ================================================================ *)

ParitySectorIndices[j_] :=
  Module[{d = 2 j + 1},
    If[IntegerQ[j] && OddQ[j],
      (* Odd integer j: index 1 has m=j (odd), so Even-m starts at index 2 *)
      <|"Even" -> Range[2, d, 2], "Odd" -> Range[1, d, 2]|>,
      (* Even integer j or half-integer j: index 1 starts the first sector *)
      <|"Even" -> Range[1, d, 2], "Odd" -> Range[2, d, 2]|>
    ]
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
    (* Reconstruct full-dimensional eigenvectors *)
    fullVecsEven = Table[
      Module[{vec = ConstantArray[0. + 0. I, d]},
        vec[[idx["Even"]]] = vecsEven[[i]]; vec],
      {i, Length[idx["Even"]]}];
    fullVecsOdd = Table[
      Module[{vec = ConstantArray[0. + 0. I, d]},
        vec[[idx["Odd"]]] = vecsOdd[[i]]; vec],
      {i, Length[idx["Odd"]]}];
    <|"Even" -> <|"Quasienergies" -> esEven[[1]],
                  "Vectors" -> vecsEven,
                  "FullVectors" -> fullVecsEven|>,
      "Odd"  -> <|"Quasienergies" -> esOdd[[1]],
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


(* ================================================================ *)
(* 4. COHERENT STATES                                               *)
(* ================================================================ *)

generateCoherentStateCompiler[] :=
  Compile[{{J, _Integer}, {q, _Real}, {p, _Real}},
    Module[{dim, \[Alpha], result, mVals, logs, maxLog, terms, norm},
      dim = 2 J + 1;
      result = ConstantArray[0.0 + 0. I, dim];
      \[Alpha] = (q + I p) / Sqrt[4 - (q^2 + p^2)];
      If[Abs[\[Alpha]] == 0,
        result[[dim]] = 1.0 + 0. I,  (* south pole at last index *)
        mVals = Range[J, -J, -1];    (* descending basis *)
        logs = 0.5 (LogGamma[2 J + 1] - LogGamma[J + 1 + mVals] -
               LogGamma[J + 1 - mVals]) + (J + mVals) Log[\[Alpha]];
        maxLog = Max[Re[logs]];
        terms = Exp[logs - maxLog];
        norm = Sqrt[Total[Abs[terms]^2]];
        result = terms / norm;
      ];
      result
    ],
    CompilationTarget -> "WVM", RuntimeOptions -> "Speed",
    Parallelization -> True, RuntimeAttributes -> {Listable}
  ];


(* ================================================================ *)
(* 5. SPECTRAL STATISTICS                                           *)
(* ================================================================ *)

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

MeanLevelSpacingRatio[eigenvalues_] :=
  Mean[Min /@ Transpose[{#, 1 / #}] &[
    Ratios[Differences[Sort[eigenvalues]]]]];

GetFluctuations[specData_Association] :=
  Module[{x = specData["Spectrum"], n = specData["N"], raw},
    raw = Range[n] - x;
    <|"x" -> x, "delta" -> raw - Mean[raw], "N" -> n|>
  ];

LengthSpectrum[fluc_Association] :=
  Module[{n, nHalf, ft, periods},
    n = fluc["N"]; nHalf = Floor[n / 2];
    ft = Fourier[fluc["delta"], FourierParameters -> {1, -1}];
    periods = Range[0, nHalf - 1];
    <|"Period" -> periods, "Amplitude" -> Re[ft[[1 ;; nHalf]]]|>
  ];

FluctuationPowerSpectrum[fluc_Association] :=
  Module[{ft, n, nHalf, freqs},
    n = fluc["N"]; nHalf = Floor[n / 2];
    ft = Abs[Fourier[fluc["delta"]]]^2;
    freqs = Range[0, n - 1] / n;
    <|"Freq" -> freqs[[1 ;; nHalf]], "Power" -> ft[[1 ;; nHalf]]|>
  ];

IPR[vec_] := Total[Abs[vec]^4] / Total[Abs[vec]^2]^2;


(* ================================================================ *)
(* 6. ANALYTICAL RMT PREDICTIONS (COE)                             *)
(* ================================================================ *)

WignerCOE[s_] := (Pi / 2) s Exp[-Pi s^2 / 4];
IntWignerCOE[s_] := 1 - Exp[-Pi s^2 / 4];

PrCOE[r_]     := (27 / 4) (r + r^2) / (1 + r + r^2)^(5/2);
PrCUE[r_]     := (81 Sqrt[3] / (4 Pi)) (r + r^2)^2 / (1 + r + r^2)^4;
PrPoisson[r_] := 2 / (1 + r)^2;

COEAsymptoticSigma2[L_] :=
  (2 / Pi^2) (Log[2 Pi L] + 1 + EulerGamma - Pi^2 / 8);

COESFFAnalytical[tau_] :=
  Piecewise[{
    {2 tau - tau Log[1 + 2 tau], 0 < tau < 1},
    {2 - tau Log[(2 tau + 1) / (2 tau - 1)], tau >= 1}
  }];


(* ================================================================ *)
(* 7. NUMBER VARIANCE                                               *)
(* ================================================================ *)

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

SingleSpectrumSigma2[specData_Association, lValues_List] :=
  Module[{results},
    results = EvaluateSpectrumSigma2[specData, lValues];
    <|"LValues" -> results[[All, 1]],
      "Sigma2Theor" -> results[[All, 2]],
      "Sigma2Sample" -> results[[All, 3]]|>
  ];


(* ================================================================ *)
(* 8. SPECTRAL FORM FACTOR                                         *)
(* ================================================================ *)

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

SingleSFF[specData_Association, tauMax_Real, nPts_Integer] :=
  Module[{taus},
    taus = N[Range[1, nPts] (tauMax / nPts)];
    <|"Tau" -> taus, "K" -> CompileSFF[specData["Spectrum"], taus]|>
  ];


(* ================================================================ *)
(* 9. COE RANDOM MATRICES                                          *)
(* ================================================================ *)

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


(* ================================================================ *)
(* 10. CLASSICAL KICKED TOP                                         *)
(* ================================================================ *)

ClassicalToSphere = Compile[{{q, _Real, 1}},
  Module[{b2 = q[[1]]^2 + q[[2]]^2},
    {q[[1]] Sqrt[1 - 0.25 b2], q[[2]] Sqrt[1 - 0.25 b2], 0.5 b2 - 1.0}
  ], CompilationTarget -> "WVM"];

ClassicalToDisk = Compile[{{s, _Real, 1}},
  Module[{b = Sqrt[2.0 / (1.0 - s[[3]])]},
    {s[[1]] b, s[[2]] b}
  ], CompilationTarget -> "WVM"];

ClassicalStep = Compile[{{s, _Real, 1}, {alpha, _Real}, {k, _Real}, {order, _Integer}},
  Module[{sx = s[[1]], sy = s[[2]], sz = s[[3]],
          ca = Cos[alpha], sa = Sin[alpha], sxc, sxs, r1, r2, r3},
    If[order == 0,
      r1 = ca sx - sa sy; r2 = sa sx + ca sy; r3 = sz;
      sx = r1; sy = r2; sz = r3;
      sxc = Cos[k sx]; sxs = Sin[k sx];
      r1 = sx; r2 = sxc sy - sxs sz; r3 = sxs sy + sxc sz;,
      sxc = Cos[k sx]; sxs = Sin[k sx];
      r1 = sx; r2 = sxc sy - sxs sz; r3 = sxs sy + sxc sz;
      sx = r1; sy = r2; sz = r3;
      r1 = ca sx - sa sy; r2 = sa sx + ca sy; r3 = sz;
    ];
    {r1, r2, r3}
  ], CompilationTarget -> "WVM"];

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
  CompilationTarget -> "WVM",
  CompilationOptions -> {"InlineCompiledFunctions" -> True},
  RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"
];


(* ================================================================ *)

End[];
EndPackage[];
