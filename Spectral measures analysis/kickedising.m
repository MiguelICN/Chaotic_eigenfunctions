(* ::Package:: *)

(* ::Title:: *)
(*Setup*)


(* ::Input:: *)
(*(* ================================================================ *)*)
(*(* MODULE 0: SETUP                                                  *)*)
(*(* ================================================================ *)*)
(**)
(*baseDir = NotebookDirectory[];*)
(*Get[FileNameJoin[{baseDir, "QKT_full.wl"}]];*)
(*$HistoryLength = 0;*)
(*EnsureDir[dir_] := If[!DirectoryQ[dir],*)
(*  CreateDirectory[dir, CreateIntermediateDirectories -> True]];*)
(*kStr[k_] := StringReplace[*)
(*  If[StringEndsQ[#, "."], # <> "0", #] &[ToString[N[k]]], "." -> "_"];*)


(* ::Input:: *)
(*LaunchKernels[14];*)


(* ::Title:: *)
(*Kicked ising*)


(* ::Input:: *)
(*Off[MakeExpression::boxfmt]*)
(**)
(*(* ============================================================ *)*)
(*(*  Mixed-Field Ising Model \[LongDash] Floquet texture test               *)*)
(*(* ============================================================ *)*)
(**)
(*L  = 8;            (* number of spins, dim = 2^L = 1024 *)*)
(*J  = 1.0;*)
(*hz = 0.8090;        (* Kim-Huse chaotic parameters *)*)
(*hx = 0.9045;*)
(*tau = 1.0;*)
(*dim = 2^L;*)


(* ::Input:: *)
(*(* Pauli matrices *)*)
(*Id2 = IdentityMatrix[2];*)
(*sx  = {{0, 1}, {1, 0}};*)
(*sz  = {{1, 0}, {0, -1}};*)
(**)
(*(* Single-site operator embedded in the L-site space *)*)
(*op[mat_, i_, len_] := *)
(*  KroneckerProduct @@ Table[If[k == i, mat, Id2], {k, 1, len}];*)
(**)
(*(* Diagonal part: Ising + longitudinal field, open BC *)*)
(*Hz = -J Sum[op[sz, i, L] . op[sz, i + 1, L], {i, 1, L - 1}]- hz Sum[op[sz, i, L], {i, 1, L}];*)
(**)
(*(* Off-diagonal part: transverse field *)*)
(*Hx = -hx Sum[op[sx, i, L], {i, 1, L}];*)


(* ::Input:: *)
(*diagHz   = Normal[Diagonal[Hz]];*)
(*expHz    = DiagonalMatrix[Exp[-I tau diagHz]];*)
(*sqrtHz   = DiagonalMatrix[Exp[-I tau diagHz/2]];*)
(*expHx    = MatrixExp[-I tau Hx];*)
(**)
(*UF    = expHz . expHx;*)
(*UFsym = sqrtHz . expHx . sqrtHz;*)


(* ::Input:: *)
(*{eval1, vec1} = Eigensystem[UF];*)
(*{eval2, vec2} = Eigensystem[UFsym];*)


(* ::Input:: *)
(*vec1[[1]][[1;;-1;;50]]//Im//Chop*)


(* ::Input:: *)
(*vec2[[1]][[1;;-1;;50]]//Im//Chop*)


(* ::Input:: *)
(*f1 = ConstantArray[1./Sqrt[dim], dim];*)
(**)
(*P1 = Abs[vec1 . f1]^2;*)
(*P2 = Abs[vec2 . f1]^2;*)
(**)
(*sat1 = dim * Total[P1^2];*)
(*sat2 = dim * Total[P2^2];*)


(* ::Input:: *)
(*Print["Asymmetric U_F  :  \[CapitalSigma]\[Infinity] = ", sat1, "   (expect ~2, CUE-like)"];*)
(*Print["Symmetric U_F   :  \[CapitalSigma]\[Infinity] = ", sat2, "   (expect ~3, COE)"];*)
(**)
(*(* Optional sanity check: maximum imaginary component of UFsym eigenvectors *)*)
(*Print["max |Im(vec_sym)| = ", Max[Abs[Im[vec2]]],*)
(*      "   (should be ~ machine epsilon up to a global phase)"];*)


(* ::Chapter:: *)
(*Change of basis*)


(* ::Input:: *)
(*gaugeSweep[theta_] := Module[{D, U, V, P},*)
(*  D = DiagonalMatrix[Exp[-I theta diagHz]];*)
(*  U = D . UFsym . ConjugateTranspose[D];*)
(*  V = Eigenvectors[U];*)
(*  P = Abs[V . f1]^2;*)
(*  dim * Total[P^2]*)
(*];*)


(* ::Input:: *)
(*sweep = Table[{theta, gaugeSweep[theta]}, {theta, 0, tau/2, tau/40}];*)


(* ::Input:: *)
(*ListLinePlot[sweep,*)
(*  PlotMarkers -> Automatic,*)
(*  Frame -> True,*)
(*  FrameLabel -> {"\[Theta]", "\[CapitalSigma]\[Infinity]"},*)
(*  GridLines -> {None, {2, 3}},*)
(*  PlotRange -> All,*)
(*  PlotLabel -> "U(\[Theta]) = exp(-i\[Theta] H_z) \[CenterDot] U_sym \[CenterDot] exp(i\[Theta] H_z)",*)
(*  ImageSize -> 600*)
(*]*)


(* ::Input:: *)
(*sweep*)


(* ::Chapter:: *)
(*Entanglement*)


(* ::Input:: *)
(*EntanglementEntropySVD[psi_, dA_, dB_] :=*)
(*  Module[{mat, sv2},*)
(*    mat = ArrayReshape[psi, {dA, dB}];*)
(*    sv2 = SingularValueList[mat]^2;*)
(*    (* Hard threshold at machine epsilon level *)*)
(*    sv2 = Select[sv2, # > 1.*^-14 &];*)
(*    (* Renormalize to absorb discarded numerical noise *)*)
(*    sv2 = sv2 / Total[sv2];*)
(*    -sv2 . Log[sv2]*)
(*  ]*)


(* ::Input:: *)
(*entanglement=ParallelTable[EntanglementEntropySVD[i, L/2,L/2 ],{i,vec1}];*)


(* ::Input:: *)
(*entanglement2=ParallelTable[EntanglementEntropySVD[i, L/2,L/2 ],{i,vec2}];*)


(* ::Input:: *)
(*ListPlot[Transpose[{Arg[eval1],entanglement}],PlotRange->All]*)


(* ::Input:: *)
(*ListPlot[Transpose[{Arg[eval2],entanglement2}],PlotRange->All]*)


(* ::Input:: *)
(*PageEntropy[La_,Lb_]:=N[PolyGamma[0,2^La*2^Lb+1]-PolyGamma[0,2^Lb+1]-(2^La-1)/(2*2^Lb)]*)


(* ::Input:: *)
(*La=4;*)
(*Lb=4;*)


(* ::Input:: *)
(*PageEntropy[La,Lb]//N*)


(* ::Input:: *)
(*Print["dA, dB = ",{2^La,2^Lb}];*)
(*Print["Max possible entropy = ",N[Log[Min[2^La,2^Lb]]]];*)
(*Print["Page value = ",PageEntropy[La,Lb]];*)


(* ::Input:: *)
(*Print["Mean Schmidt rank: ",Mean[Length[Select[SingularValueList[ArrayReshape[#,{2^La,2^Lb}]]^2,#>10^-14&]]&/@vec1]];*)


(* ::Input:: *)
(*phases=Mod[Arg[eval1],2 Pi];*)
(*spacings=Differences[Sort[phases]];*)
(*ratios=Min[#1/#2,#2/#1]&@@@Partition[spacings,2,1];*)
(*Mean[ratios]*)


(* ::Input:: *)
(*{J,hz,hx,tau} (*print these to verify*)*)


(* ::Input:: *)
(*sv=SingularValueList[ArrayReshape[vec1[[1]],{64,64}]];*)
(*ListLogPlot[sv^2,PlotRange->All]*)
