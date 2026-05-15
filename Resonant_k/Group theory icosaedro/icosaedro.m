(* ::Package:: *)

(* ::Title:: *)
(*Generalized QKT*)


(* ::Chapter::Closed:: *)
(*Setup*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]];*)


(* ::Input:: *)
(*Get["C:\\Users\\Miguel\\Github\\Chaotic_eigenfunctions\\Resonant_k\\QKT_full.wl"];*)


(* ::Input:: *)
(*LaunchKernels[14];*)


(* ::Chapter::Closed:: *)
(*Grid*)


(* ::Input:: *)
(*radius=2;*)
(*delta=0.02;*)
(*xCoords=Range[-radius,radius,delta];*)
(*yCoords=Range[-radius,radius,delta];*)
(*gridPointsRaw=Flatten[Table[{x,y},{x,xCoords},{y,yCoords}],1];*)
(*(*Keep only those inside disk*)*)
(*region=Select[gridPointsRaw,Norm[#]<radius&];*)
(*region//Dimensions*)


(* ::Input:: *)
(*GetQKTParityEigensystem[jParam_, alphaParam_, kParam_] :=*)
(*  Module[{mat, d, idx, matEven, matOdd, esEven, esOdd,*)
(*          phEven, phOdd, vecsEven, vecsOdd, idxE, idxO,*)
(*          fullVecsEven, fullVecsOdd},*)
(*    mat = Floqsym[jParam, alphaParam, kParam];*)
(*    d = 2 jParam + 1;*)
(*    idx = ParitySectorIndices[jParam];*)
(*    matEven = mat[[idx["Even"], idx["Even"]]];*)
(*    matOdd  = mat[[idx["Odd"], idx["Odd"]]];*)
(*    esEven = Eigensystem[matEven, Method -> "Direct"];*)
(*    esOdd  = Eigensystem[matOdd, Method -> "Direct"];*)
(*    phEven = Mod[Arg[esEven[[1]]], 2 Pi];*)
(*    phOdd  = Mod[Arg[esOdd[[1]]], 2 Pi];*)
(*    idxE = Ordering[phEven]; idxO = Ordering[phOdd];*)
(*    vecsEven = esEven[[2, idxE]];*)
(*    vecsOdd  = esOdd[[2, idxO]];*)
(*    fullVecsEven = Table[*)
(*      Module[{vec = ConstantArray[0. + 0. I, d]},*)
(*        vec[[idx["Even"]]] = vecsEven[[i]]; vec],*)
(*      {i, Length[idx["Even"]]}];*)
(*    fullVecsOdd = Table[*)
(*      Module[{vec = ConstantArray[0. + 0. I, d]},*)
(*        vec[[idx["Odd"]]] = vecsOdd[[i]]; vec],*)
(*      {i, Length[idx["Odd"]]}];*)
(*    <|"Even" -> <|"Quasienergies" -> esEven[[1]],*)
(*                  "Phases" -> phEven,*)
(*                  "Vectors" -> vecsEven,*)
(*                  "FullVectors" -> fullVecsEven|>,*)
(*      "Odd"  -> <|"Quasienergies" -> esOdd[[1]], *)
(*                  "Phases" -> phOdd,*)
(*                  "Vectors" -> vecsOdd,*)
(*                  "FullVectors" -> fullVecsOdd|>|>*)
(*  ];*)


(* ::Chapter:: *)
(*Method*)


(* ::Input:: *)
(*(*Floquet parameters*)*)
(*J=100;*)
(*\[Alpha]=ArcTan[2.];*)
(*k=J*Pi(2/5.);*)


(* ::Input:: *)
(*ops=generateSpinOperators[J];*)


(* ::Input:: *)
(*Print["Generating coherent states over ",Length[region]," points..."];*)
(*coherentStateGenerator=generateCoherentStateCompiler[];*)
(*statesCompiled=ParallelMap[coherentStateGenerator[J,#[[1]],#[[2]]]&,region,Method->"CoarsestGrained"  (*Optional:better batching*)];*)


(* ::Input:: *)
(*Print["Computing Floquet matrix..."];*)
(*data=GetQKTParityEigensystem[J,\[Alpha],k];*)
(*neweigenvec=Join[data["Even","FullVectors"],data["Odd","FullVectors"]];*)
(*Print["Projecting states into Floquet basis..."];*)
(*overlaps=Conjugate[neweigenvec] . Transpose[statesCompiled];*)
(*Print["Computing IPRs..."];*)
(*iprCompiled=Total[Abs[overlaps]^4,{1}];*)


(* ::Input:: *)
(*newdata=Table[-Log[i]/Log[2 J+1],{i,iprCompiled}];*)
(*plot=ListDensityPlot[Transpose[{region[[All,1]],region[[All,2]],newdata}],PlotRange->All,PlotTheme->"Detailed",ColorFunction->"SunsetColors",PlotLabel->Style["D\:2082, J = "<>ToString[J]<>", \[Alpha] = "<>ToString[N[\[Alpha]]]<>", k = "<>ToString[NumberForm[k,{4,2}]],25,Black],ImageSize->600,FrameLabel->{Style["Q",40,Black],Rotate[Style["P",40,Black],-90 Degree]},LabelStyle->Directive[30,Black],AspectRatio->1,Frame->True,InterpolationOrder->3]*)
(*Export["d2_J_"<>ToString[J]<>"alpha_"<>ToString[\[Alpha]]<>"_k_"<>ToString[k]<>"_.png",plot]*)


(* ::Chapter:: *)
(*Husimis*)


(* ::Input:: *)
(*dir="HusimiPlots";*)
(*If[!DirectoryQ[dir],CreateDirectory[dir]];*)


(* ::Input:: *)
(*husimiValues=Abs[overlaps]^2;*)
(*husimiValuesSorted=husimiValues;*)


(* ::Input:: *)
(*(*should be 2J+1*)*)
(*plots=Table[*)
(*Module[{husimiList,plot},*)
(*husimiList=MapThread[Append,{region,husimiValuesSorted[[q]]}];*)
(*plot=ListDensityPlot[husimiList,PlotRange->All,PlotTheme->"Detailed",ColorFunction->"BlueGreenYellow",PlotLabel->Style["Husimi, J = "<>ToString[J]<>", \[Alpha] = "<>ToString[N[\[Alpha]]]<>", k = "<>ToString[NumberForm[k,{4,2}]]<>", Eigenvec_q = "<>ToString[q],25,Black],ImageSize->1000,FrameLabel->{Style["Q",40,Black],Rotate[Style["P",40,Black],-90 Degree]},LabelStyle->Directive[30,Black],AspectRatio->1,Frame->True,InterpolationOrder->1];*)
(*Export[FileNameJoin[{dir,"Husimi_J_"<>ToString[J]<>"_q"<>ToString[q]<>".png"}],plot];*)
(*plot],{q,2J+1}];*)
