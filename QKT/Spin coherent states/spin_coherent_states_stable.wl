(* ============================================================================ *)
(* NUMERICALLY STABLE SPIN COHERENT STATES FOR THE QKT                        *)
(* Four alternative implementations + optimized original                       *)
(* ============================================================================ *)

(* ---- Spin operators (unchanged) ---- *)
ClearAll[sz, sp, sm, sx, sy, sx2];
sz[j_] := sz[j] = DiagonalMatrix[Table[-j + k - 1, {k, 1, 2 j + 1}]];
sp[j_] := sp[j] = DiagonalMatrix[Table[Sqrt[j (j + 1) - m (m + 1)], {m, -j, j - 1}], 1];
sm[j_] := sm[j] = DiagonalMatrix[Table[Sqrt[j (j + 1) - m (m - 1)], {m, -j + 1, j}], -1];
sx[j_] := sx[j] = SparseArray[(1/2) (sp[j] + sm[j])];
sy[j_] := sy[j] = SparseArray[(1/(2 I)) (sp[j] - sm[j])];
sx2[j_] := sx2[j] = sx[j] . sx[j];


(* ============================================================================ *)
(* METHOD 1: TRIGONOMETRIC REPRESENTATION (Arecchi et al. 1972)                *)
(*           with separated magnitude/phase and peak-only computation          *)
(*           This is the RECOMMENDED primary method.                           *)
(* ============================================================================ *)

generateCoherentStateTrig[] := Compile[
  {{J, _Integer}, {q, _Real}, {p, _Real}},
  Module[{dim, r2, cosHalf, sinHalf, phi, logMag, phase, 
          mVals, maxLogMag, magnitudes, result, norm,
          m0, mLow, mHigh, width, logBinom, sign},
    dim = 2 J + 1;
    r2 = q^2 + p^2;
    result = ConstantArray[0.0 + 0.0 I, dim];
    
    If[r2 < 1.0*^-30,
      (* South pole: |J,-J> *)
      result[[1]] = 1.0 + 0.0 I,
      
      If[r2 > 4.0 - 1.0*^-12,
        (* North pole: |J,+J> *)
        result[[dim]] = 1.0 + 0.0 I,
        
        (* General case *)
        cosHalf = Sqrt[4.0 - r2] / 2.0;
        sinHalf = Sqrt[r2] / 2.0;
        phi = ArcTan[q, p];
        
        (* Peak location and width *)
        (* m0 = -J*cosTheta = -J*(1 - r2/2) *)
        m0 = Round[-J (1.0 - r2/2.0)];
        If[m0 < -J, m0 = -J];
        If[m0 > J, m0 = J];
        width = Ceiling[8.0 Sqrt[J * sinHalf * cosHalf + 1.0]];
        mLow = Max[-J, m0 - width];
        mHigh = Min[J, m0 + width];
        
        (* Compute only significant components *)
        Do[
          logBinom = 0.5 (LogGamma[2.0 J + 1.0] - 
                          LogGamma[J + 1.0 + m] - 
                          LogGamma[J + 1.0 - m]);
          logMag = logBinom + 
                   (J - m) * Log[cosHalf] + 
                   (J + m) * Log[sinHalf];
          
          (* Phase: (-1)^{J+m} * e^{-i(J+m)*phi} *)
          (* (-1)^{J+m} = e^{i*pi*(J+m)} *)
          phase = (J + m) * (Pi - phi);
          
          result[[J + m + 1]] = Exp[logMag] (Cos[phase] + I Sin[phase]);
          , {m, mLow, mHigh}
        ];
        
        norm = Sqrt[Total[Abs[result]^2]];
        If[norm > 0.0, result = result / norm];
      ];
    ];
    result
  ],
  CompilationTarget -> "WVM",
  RuntimeAttributes -> {Listable},
  Parallelization -> True
];


(* ============================================================================ *)
(* METHOD 2: RECURRENCE FROM PEAK (Varshalovich et al. 1988)                   *)
(*           Avoids LogGamma entirely except for the seed value.               *)
(* ============================================================================ *)

generateCoherentStateRecurrence[] := Compile[
  {{J, _Integer}, {q, _Real}, {p, _Real}},
  Module[{dim, r2, tanHalf, phi, expIPhi, m0, 
          logSeed, result, ratio, cm, norm,
          cosHalf, sinHalf, mLow, mHigh, width},
    dim = 2 J + 1;
    r2 = q^2 + p^2;
    result = ConstantArray[0.0 + 0.0 I, dim];
    
    If[r2 < 1.0*^-30,
      result[[1]] = 1.0 + 0.0 I,
      
      If[r2 > 4.0 - 1.0*^-12,
        result[[dim]] = 1.0 + 0.0 I,
        
        cosHalf = Sqrt[4.0 - r2] / 2.0;
        sinHalf = Sqrt[r2] / 2.0;
        tanHalf = sinHalf / cosHalf;
        phi = ArcTan[q, p];
        expIPhi = Cos[phi] + I Sin[phi];
        
        (* Peak location *)
        m0 = Round[-J (1.0 - r2/2.0)];
        If[m0 < -J, m0 = -J];
        If[m0 > J, m0 = J];
        width = Ceiling[8.0 Sqrt[J * sinHalf * cosHalf + 1.0]];
        mLow = Max[-J, m0 - width];
        mHigh = Min[J, m0 + width];
        
        (* Seed at m0: compute c_{m0} in log-space *)
        logSeed = 0.5 (LogGamma[2.0 J + 1.0] - 
                        LogGamma[J + 1.0 + m0] - 
                        LogGamma[J + 1.0 - m0]) +
                  (J - m0) * Log[cosHalf] + 
                  (J + m0) * Log[sinHalf];
        
        (* Seed value with phase *)
        cm = Exp[logSeed] * (Cos[(J + m0)(Pi - phi)] + I Sin[(J + m0)(Pi - phi)]);
        result[[J + m0 + 1]] = cm;
        
        (* Recurrence UPWARD: m0 -> mHigh *)
        (* c_{m+1}/c_m = sqrt((J-m)/(J+m+1)) * tan(theta/2) * e^{i(pi-phi)} *)
        cm = result[[J + m0 + 1]];
        Do[
          ratio = Sqrt[(J - m + 0.0) / (J + m + 1.0)] * tanHalf;
          cm = cm * ratio * (-expIPhi);  (* -e^{i*phi} from the (-1) factor *)
          result[[J + m + 1]] = cm;
          , {m, m0 + 1, mHigh}
        ];
        
        (* Recurrence DOWNWARD: m0 -> mLow *)
        cm = result[[J + m0 + 1]];
        Do[
          ratio = Sqrt[(J + m + 0.0) / (J - m + 1.0)] / tanHalf;
          cm = cm * ratio / (-expIPhi);
          result[[J + m + 1]] = cm;
          , {m, m0 - 1, mLow, -1}
        ];
        
        norm = Sqrt[Total[Abs[result]^2]];
        If[norm > 0.0, result = result / norm];
      ];
    ];
    result
  ],
  CompilationTarget -> "WVM",
  RuntimeAttributes -> {Listable},
  Parallelization -> True
];


(* ============================================================================ *)
(* METHOD 3: KRYLOV SUBSPACE (Matrix Exponential)                              *)
(*           (Yan, Wang & Robnik 2024 / Sidje 1998)                            *)
(*           Uses MatrixExp acting on |J,-J> via sparse Jy or J_n.             *)
(*           NOT compiled — uses Mathematica's built-in sparse MatrixExp.       *)
(* ============================================================================ *)

(* Precompute the sparse rotation generator for a given J *)
buildRotationGenerator[j_] := Module[{sxMat, syMat},
  sxMat = SparseArray[(1/2) (sp[j] + sm[j])];
  syMat = SparseArray[(1/(2 I)) (sp[j] - sm[j])];
  {sxMat, syMat}
];

generateCoherentStateKrylov[j_, q_, p_] := Module[
  {r2, theta, phi, nVec, generator, southPole, result,
   sxMat, syMat, dim},
  dim = 2 j + 1;
  r2 = q^2 + p^2;
  
  If[r2 < 1.0*^-30,
    (* South pole *)
    SparseArray[{1 -> 1.0 + 0.0 I}, dim],
    
    If[r2 > 4.0 - 1.0*^-12,
      (* North pole *)
      SparseArray[{dim -> 1.0 + 0.0 I}, dim],
      
      (* General case *)
      theta = 2.0 ArcSin[Sqrt[r2] / 2.0];
      phi = ArcTan[q, p];
      
      {sxMat, syMat} = buildRotationGenerator[j];
      
      (* Generator: -i*theta*(sin(phi)*Jx - cos(phi)*Jy) *)
      generator = -I theta (Sin[phi] sxMat - Cos[phi] syMat);
      
      (* |J,-J> = first basis vector *)
      southPole = SparseArray[{1 -> 1.0 + 0.0 I}, dim];
      
      (* Apply matrix exponential via Krylov *)
      (* Mathematica's MatrixExp[M, v] uses Krylov internally for sparse M *)
      result = MatrixExp[N[generator], N[southPole]];
      result / Sqrt[Total[Abs[result]^2]]
    ]
  ]
];

(* Parallelized version for the full grid *)
generateCoherentStatesKrylov[j_, gridPoints_] := Module[
  {sxMat, syMat},
  {sxMat, syMat} = buildRotationGenerator[j];
  ParallelMap[
    Module[{q = #[[1]], p = #[[2]], r2, theta, phi, gen, sp0, res, dim},
      dim = 2 j + 1;
      r2 = q^2 + p^2;
      If[r2 < 1.0*^-30,
        SparseArray[{1 -> 1.0 + 0.0 I}, dim] // Normal,
        If[r2 > 4.0 - 1.0*^-12,
          SparseArray[{dim -> 1.0 + 0.0 I}, dim] // Normal,
          theta = 2.0 ArcSin[Sqrt[r2] / 2.0];
          phi = ArcTan[q, p];
          gen = -I theta (Sin[phi] sxMat - Cos[phi] syMat);
          sp0 = SparseArray[{1 -> 1.0 + 0.0 I}, dim];
          res = MatrixExp[N[gen], N[sp0]];
          Normal[res / Sqrt[Total[Abs[res]^2]]]
        ]
      ]
    ] &,
    gridPoints,
    Method -> "CoarsestGrained"
  ]
];


(* ============================================================================ *)
(* METHOD 4: OPTIMIZED ORIGINAL (fixes identified bugs)                        *)
(*           Key changes:                                                       *)
(*           - Removed RuntimeOptions -> "Speed"                                *)
(*           - Separated Re/Im computation                                      *)
(*           - Added peak-only computation                                      *)
(*           - Uses real-valued log for magnitude                               *)
(* ============================================================================ *)

generateCoherentStateOptimized[] := Compile[
  {{J, _Integer}, {q, _Real}, {p, _Real}},
  Module[{dim, r2, logR, logFactor, phi, m0, width, mLow, mHigh,
          logBinom, logMagAlpha, logMag, phase, maxLog,
          result, norm, mags, phases},
    dim = 2 J + 1;
    r2 = q^2 + p^2;
    result = ConstantArray[0.0 + 0.0 I, dim];
    
    If[r2 < 1.0*^-30,
      result[[1]] = 1.0 + 0.0 I,
      
      If[r2 > 4.0 - 1.0*^-12,
        result[[dim]] = 1.0 + 0.0 I,
        
        (* Decompose Log|alpha| = 0.5*Log(r2/(4-r2)) *)
        logMagAlpha = 0.5 Log[r2 / (4.0 - r2)];
        phi = ArcTan[q, p];
        
        (* Peak and width *)
        m0 = Round[-J (1.0 - r2 / 2.0)];
        If[m0 < -J, m0 = -J];
        If[m0 > J, m0 = J];
        width = Ceiling[8.0 Sqrt[J / 2.0 + 1.0]];
        mLow = Max[-J, m0 - width];
        mHigh = Min[J, m0 + width];
        
        (* First pass: compute log-magnitudes to find maximum *)
        mags = ConstantArray[-1.0*^300, mHigh - mLow + 1];
        phases = ConstantArray[0.0, mHigh - mLow + 1];
        Do[
          logBinom = 0.5 (LogGamma[2.0 J + 1.0] - 
                          LogGamma[J + 1.0 + m] - 
                          LogGamma[J + 1.0 - m]);
          mags[[m - mLow + 1]] = logBinom + (J + m) * logMagAlpha;
          phases[[m - mLow + 1]] = (J + m) * phi;
          , {m, mLow, mHigh}
        ];
        
        maxLog = Max[mags];
        
        (* Second pass: exponentiate and assemble *)
        Do[
          logMag = mags[[m - mLow + 1]] - maxLog;
          If[logMag > -36.0,  (* Only compute if above machine epsilon *)
            phase = phases[[m - mLow + 1]];
            result[[J + m + 1]] = Exp[logMag] (Cos[phase] + I Sin[phase]);
          ];
          , {m, mLow, mHigh}
        ];
        
        norm = Sqrt[Total[Abs[result]^2]];
        If[norm > 0.0, result = result / norm];
      ];
    ];
    result
  ],
  CompilationTarget -> "WVM",
  RuntimeAttributes -> {Listable},
  Parallelization -> True
  (* NOTE: No RuntimeOptions -> "Speed" *)
];


(* ============================================================================ *)
(* VALIDATION AND BENCHMARKING                                                  *)
(* ============================================================================ *)

validateCoherentStates[J_, nTests_: 100] := Module[
  {testPoints, gen1, gen2, gen4, states1, states2, states4,
   maxErr12, maxErr14, maxErr24, t1, t2, t4,
   norms1, norms2, norms4},
  
  (* Generate random test points on the disk *)
  SeedRandom[42];
  testPoints = Select[
    Table[{RandomReal[{-1.99, 1.99}], RandomReal[{-1.99, 1.99}]}, nTests],
    Norm[#] < 1.99 &
  ];
  Print["Test points: ", Length[testPoints]];
  
  (* Method 1: Trigonometric *)
  gen1 = generateCoherentStateTrig[];
  t1 = First@AbsoluteTiming[
    states1 = Map[gen1[J, #[[1]], #[[2]]] &, testPoints];
  ];
  
  (* Method 2: Recurrence *)
  gen2 = generateCoherentStateRecurrence[];
  t2 = First@AbsoluteTiming[
    states2 = Map[gen2[J, #[[1]], #[[2]]] &, testPoints];
  ];
  
  (* Method 4: Optimized original *)
  gen4 = generateCoherentStateOptimized[];
  t4 = First@AbsoluteTiming[
    states4 = Map[gen4[J, #[[1]], #[[2]]] &, testPoints];
  ];
  
  (* Check norms *)
  norms1 = Map[Sqrt[Total[Abs[#]^2]] &, states1];
  norms2 = Map[Sqrt[Total[Abs[#]^2]] &, states2];
  norms4 = Map[Sqrt[Total[Abs[#]^2]] &, states4];
  
  (* Cross-compare: max |⟨ψ_a|ψ_b⟩| - 1 *)
  (* For identical states, the overlap magnitude should be 1 *)
  maxErr12 = Max[MapThread[1 - Abs[Conjugate[#1] . #2] &, {states1, states2}]];
  maxErr14 = Max[MapThread[1 - Abs[Conjugate[#1] . #2] &, {states1, states4}]];
  maxErr24 = Max[MapThread[1 - Abs[Conjugate[#2] . #4] &, {states2, states4}]];
  
  Print["=== J = ", J, " ==="];
  Print["Method 1 (Trig):       ", t1, " s,  norm range: ", 
        {Min[norms1], Max[norms1]}];
  Print["Method 2 (Recurrence): ", t2, " s,  norm range: ", 
        {Min[norms2], Max[norms2]}];
  Print["Method 4 (Optimized):  ", t4, " s,  norm range: ", 
        {Min[norms4], Max[norms4]}];
  Print["Max overlap error (1 vs 2): ", maxErr12];
  Print["Max overlap error (1 vs 4): ", maxErr14];
  Print["Max overlap error (2 vs 4): ", maxErr24];
];

(* Run: validateCoherentStates[500] then validateCoherentStates[5000] *)


(* ============================================================================ *)
(* FULL GRID COMPUTATION (drop-in replacement for the user's workflow)         *)
(* ============================================================================ *)

computeCoherentStatesGrid[J_, radius_: 2, delta_: 0.025] := Module[
  {region, gen, states, t},
  
  region = Select[
    Flatten[Table[{x, y}, 
      {x, Range[-radius, radius, delta]}, 
      {y, Range[-radius, radius, delta]}], 1],
    Norm[#] < radius &
  ];
  Print["Grid points: ", Length[region]];
  
  gen = generateCoherentStateTrig[];
  Print["Generating coherent states (Trig method)..."];
  
  t = First@AbsoluteTiming[
    states = ParallelMap[
      gen[J, #[[1]], #[[2]]] &, 
      region, 
      Method -> "CoarsestGrained"
    ];
  ];
  
  Print["Time: ", t, " s"];
  Print["Sample norms: ", 
    Sqrt[Total[Abs[#]^2]] & /@ RandomSample[states, Min[5, Length[states]]]];
  
  {region, states}
];

(* Usage: {region, states} = computeCoherentStatesGrid[5000]; *)
