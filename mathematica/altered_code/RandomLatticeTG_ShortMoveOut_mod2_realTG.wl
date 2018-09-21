(* ::Package:: *)

#!/usr/local/bin/MathematicaScript -script

Print ["Hello World!"]




(*One particle case
Here I define all the parameters of the system.*)
numBar=2; (* Total number of the delta barriers*)
numBar = ToExpression[$ScriptCommandLine[[2]]]
ys = Symbol["y"<>ToString[#]]&/@Range[1,numBar] (*List of the variables (placeholders) for the positions of the barriers*)
ds=Symbol["d"<>ToString[#]]&/@Range[1,numBar](*List of the variables (placeholders) for the heights of the barriers*)
specificRules = {L->10, \[HBar]->1, m->1}; (*Rules to fix the length of the box L, and the \[HBar]/m ratio*)
thisNum = ToExpression[$ScriptCommandLine[[4]]]
specificRules = {L->thisNum, \[HBar]->1, m->1};
$Assumptions=Element[ds, Reals]&&Element[ys, Reals]&&Element[L, Reals] &&Element[m, Reals]&& Element[\[HBar], Reals]&&Element[k, Reals](*General assumptions about parameters. They are all real.*)
(*Here I define cyclic permutation of the barriers with \[CapitalDelta]\[Element][-L, L] shift. Note that for negative \[CapitalDelta] the effective shift is equivalent to  L+\[CapitalDelta] shift. This is precisely what function Mod[\[CapitalDelta], L] is doing.*)
(* Cyclic permutation of the barrier heights *)
cyclicBarrierShiftHeights[\[CapitalDelta]_, ys_, ds_, L_]:=Module[{shift, signs, dels, resultDs, Nb},(Nb = Length[ds];dels=Flatten[{L/2-#&/@ Reverse[ys], L}];signs=Sign[Mod[\[CapitalDelta],L] - dels]; shift=Floor[Total[signs + 1]/2];resultDs = Flatten[{ds[[Nb-shift+1;;Nb]],ds[[1;;Nb-shift]]}])]
(* Cyclic permutation of the barrier positions *)
cyclicBarrierShiftPos[\[CapitalDelta]_, ys_, L_]:=Module[{Nb, resultYs},(Nb = Length[ds];resultYs = Sort[Mod[ys+\[CapitalDelta] + L/2, L]-L/2])]
(*The Bethe equations
Here we construct a Bethe Equation *)
\[ScriptCapitalT][n_, ys_, ds_]:= Module[{combs},(combs=Subsets[Range[1,Length[ys]], {n}];2^(2n) Total[Times@@(ds[[#]] m/\[HBar]^2)If[Length[#]>=1, Sin[k(ys[[#]][[1]]+L/2)]Sin[k(L/2-ys[[#]][[n]])], Sin[k L]]Product[Sin[k(ys[[#]][[j]]-ys[[#]][[j-1]])],{j,2,n}]&/@combs])]
\[ScriptCapitalT][1, ys, ds];
makeBE[ys_, ds_]:=Sum[\[ScriptCapitalT][n, ys, ds]k^(Length[ds]-n), {n, 0, Length[ds]}]==0
AbsoluteTiming[TGBE=makeBE[ys, ds];]
(*The Coefficients
Here I reconstruct the coefficients and the normaliation constant*)
(*n is the number of the well in which the particle is in, for 1-particle case. The coefficient we are calculating is Subscript["\[ScriptCapitalA]", {1,0,0,0,0}][{k1}], with 1 in the well number n.  ds - heights of the barriers, ys - positions of the barriers*)
makeCoef[n_, ds_, ys_]:=Module[{combs, posCombs, elems, curProd, curYP},(combs=Subsets[Range[1,n-1], n-1]; elems = Total[(Times@@(ds[[#]]*(-2 m I)/(k \[HBar]^2))If[Length[#]>0,(1-Exp[-I k (L + 2 ys[[#]][[1]])]), 1])Product[(1-Exp[-I 2 k (+ys[[#]][[j]]-ys[[#]][[j-1]])]), {j, 2, Length[#]}]&/@combs])];
makeCoef[2, ds, ys];
(*Real part of a coefficient*)
realPartMakeCoef[n_, ds_, ys_]:=Module[{combs, posCombs, elems, curProd, curYP},(combs=Subsets[Range[1,n-1], n-1]; elems = Total[(Times@@(ds[[#]]) ((-4 m )/(k \[HBar]^2))^Length[#] If[Length[#]>0,(-Sin[k(L/2+ys[[#]][[1]])]Cos[k(ys[[#]][[Length[#]]]+L/2)]), 1])Product[(Sin[ k (-ys[[#]][[j]]+ys[[#]][[j-1]])]), {j, 2, Length[#]}]&/@combs])];
realPartMakeCoef[2, ds, ys];
(*Imaginary part of a coefficient *)
imPartMakeCoef[n_, ds_, ys_]:=Module[{combs, posCombs, elems, curProd, curYP},(combs=Subsets[Range[1,n-1], {1,n-1}]; elems = Total[(Times@@(ds[[#]]) ((-4 m )/(k \[HBar]^2))^Length[#] If[Length[#]>0,(Sin[k(L/2+ys[[#]][[1]])]Sin[k(ys[[#]][[Length[#]]]+L/2)]), 1])Product[(Sin[ k (-ys[[#]][[j]]+ys[[#]][[j-1]])]), {j, 2, Length[#]}]&/@combs])];
imPartMakeCoef[2, {d1,d2,d3}, {y1,y2,y3}];
(* Calculate the normalization constant*)
calcNorm[ds_, ys_]:= Module[{ReA, ImA, ReA1, ImA1},(ReA1 = realPartMakeCoef[numBar+1, ds, ys]; ImA1=imPartMakeCoef[numBar+1, ds, ys];2 (ys[[1]]+L/2 + -(Sin[k(L+2ys[[1]])])/(2 k)) + 2 Total[(ReA =realPartMakeCoef[#, ds, ys];ImA=imPartMakeCoef[#, ds, ys]; ReA^2 (ys[[#]]-ys[[#-1]] + -(Sin[k(L+2ys[[#]])]-Sin[k (L+2ys[[#-1]])])/(2 k)) + ImA^2 (ys[[#]]-ys[[#-1]] - -(Sin[k(L+2ys[[#]])]-Sin[k (L+2ys[[#-1]])])/(2 k)) - ReA ImA (Cos[k(L+2ys[[#]])]-Cos[k (L+2ys[[#-1]])])/k)&/@Range[2, numBar]]+2(ReA1^2 (L/2-ys[[numBar]] + -(Sin[k 2 L]-Sin[k (L+2ys[[numBar]])])/(2 k)) + ImA1^2 (L/2-ys[[numBar]] - -(Sin[k 2 L]-Sin[k (L+2ys[[numBar]])])/(2 k)) - ReA1 ImA1 (Cos[k 2 L]-Cos[k (L+2ys[[numBar]])])/k))]
Print["So far so good 1"]


(*Actual solution*)
a = \[Gamma] * L/(numBar-1);
y0  = -L/2+(L-(numBar-1)a)/2;randYsAllDelta=(y0+a*(#) + \[Delta] (y0+L/2))/.specificRules&/@Range[0, numBar-1]//FullSimplify;
(*randYsAllDelta=((-L/2 + L/(2(numBar)))+L/numBar(#) + \[Delta]L/(2(numBar))(**Mod[#,2]*))/.specificRules&/@Range[0, numBar-1];*)
thumorse={0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0};
randDs =(100+0*Mod[#+0,2])&/@Range[1, numBar](* RandomReal[{0, 3}, Length[ds]]*)(*(1)&/@Range[1, numBar]*) ;
(*randDs[[{-7, -6}]] = 0;*)
randRulesAllDelta= Flatten[{ys[[#]]->randYsAllDelta[[#]],ds[[#]]->randDs[[#]]}&/@Range[1, Length[ds]]]
allYsAllDelta=ys[[#]]->randYsAllDelta[[#]]&/@Range[1, Length[ds]];
Plot[randYsAllDelta/.Flatten[{specificRules, \[Delta]->0}], {\[Gamma], 0, 3}]
(*cyclicBarrierShiftHeights[ -2, randYsAllDelta, ds, L/.specificRules];
Plot[cyclicBarrierShiftPos[\[CapitalDelta], randYsAllDelta, L/.specificRules][[1]], {\[CapitalDelta], -5,5}];*)
randYsAllDelta[[2]]
bEF = Experimental`CreateNumericalFunction[{k,g},{TGBE[[1]]/.Flatten[{randRulesAllDelta, specificRules}]/.{\[Delta]->0, \[Gamma]->g}}, {1}, Compiled->True];
bEF[{1,0.5}]
bEF["CompiledFunction"][1,0.5]

betheEquationsF[k_,g_]=TGBE[[1]]/.Flatten[{randRulesAllDelta, specificRules}]/.{\[Delta]->0, \[Gamma]->g};
bECompiled = Compile[{k, g}, betheEquationsF[k,g]]
rationalGammas = Range[0,3, 1/100];
momenta = Range[1,10,0.01];
(*butterflyValues = Flatten[Module[{g, kk},(g = rationalGammas[[#]];(kk=momenta[[#]];{g,kk,If[Abs[bEF["CompiledFunction"][kk, g][[1]]]<1^(-2), 1, 0]})&/@Range[1, Length[momenta]])&/@Range[1, Length[rationalGammas]]],1];
ListDensityPlot[butterflyValues, PlotLegends\[Rule]Automatic, InterpolationOrder\[Rule]0]*)

betheEquations=TGBE[[1]]/.Flatten[{randRulesAllDelta, specificRules}]/.{\[Delta]->\[CapitalDelta],\[Gamma]->(numBar-1)/numBar};
SetOptions[ContourPlot,BaseStyle->{FontFamily->"Times",FontSize->14}];
plotObj1=ContourPlot[{betheEquations==0},{\[CapitalDelta], -1,1}, {k, 1, 10}, PlotRange->Full, PlotLegends->Automatic, PlotPoints->100, FrameLabel->{"\[CapitalDelta]", "k"}, AspectRatio->2 , Epilog->{Thin,{Green,Point[{-0.7,8.68}],Green,Point[{-0.7,4.38}], Black, Point[{-0.4,8.72}], Black,Point[{-0.4,4.95}]}}]
(*orderedPairs=First[Cases[plotObj,x_GraphicsComplex \[RuleDelayed]  First@x,Infinity]];
orderedSqrP={#[[1]], #[[2]]^2} &/@ orderedPairs;
lineData=Cases[plotObj,x_Line \[RuleDelayed]  x,Infinity];
linesIndex = Cases[lineData,Line[x_] \[RuleDelayed]  x,Infinity];
allLines = orderedSqrP[[#]]&/@ linesIndex;
Length[orderedPairs]
Length[lineData]*)
Delta=0;
GGamma = 1/3(*(numBar-1)/numBar*);
GGamma = ToExpression[$ScriptCommandLine[[3]]];
randYs=randYsAllDelta/.Flatten[{specificRules}](*Sort[RandomReal[{-L/2, L/2}/.specificRules, Length[ys]]]*)(*{3}*);
randRules =randRulesAllDelta/.Flatten[{specificRules}]

(*plotObj=ContourPlot[Evaluate[{TGBE/.k\[Rule]k1, TGBE/.k\[Rule]k2}//.Flatten[{specificRules, \[Delta]\[Rule]Delta, \[Gamma]\[Rule]GGamma, randRules}]], {k1, -0.1, 9}, {k2, -0.1, 5}, PlotPoints\[Rule]300,PlotLegends\[Rule]Automatic, FrameLabel\[Rule]{"k1", "k2"}]*)
(*sampleroots=Cases[(*Normal@*)plotObj,(*Line[pts_,___]\[RuleDelayed]pts*)GraphicsComplex[pts__]\[RuleDelayed]pts,Infinity];
ListPlot[#[[1]]&/@sampleroots[[1]]]*)
initPt ={3.9};
rootPt = FindRoot[TGBE//.Flatten[{specificRules, \[Delta]->Delta,\[Gamma]->GGamma, randRules}],{k, #}]&/@initPt 
rootPt=NSolve[(TGBE//.Flatten[{specificRules, \[Delta]->Delta,\[Gamma]->GGamma, randRules}])&&k > 0 && k <4,k, Reals, WorkingPrecision->20]
(TGBE[[1]]//.Flatten[{specificRules, \[Delta]->Delta,\[Gamma]->GGamma, randRules, rootPt[[#]]}])&/@Range[1,Length[rootPt]]
(*The wave function*)
Acoeff[block_, P_]=Subscript["\[ScriptCapitalA]", block][P];
allBlocks[N_, barN_] := DeleteDuplicates[Flatten[Permutations[#]&/@(PadRight[#, barN+1]&/@IntegerPartitions[N, barN+1]),1]]
\[Chi][wellN_, ds_, ys_] := Module[{ Ep}, (Ep ={k,-k}; \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(ne = 1\), \(\(Dimensions[Ep]\)[\([1]\)]\)]\(
\*SuperscriptBox[\((\(-Exp[\(-I\)\ k\ L]\))\), \(If[Refine[Sign[Ep[\([ne]\)]], k > 0] == 1, 0, 1]\)] \((makeCoef[\ wellN, ds, \ ys] /. k -> Ep[\([ne]\)])\) Exp[I\ Ep[\([ne]\)]\ x]\)\))]
\[CapitalPsi][ds_,ys_, barN_] := Module[{bN, allB, N, dels, curWell, barHs, barPos,LL}, (N=1;allB=allBlocks[N, barN];
barPos=ys;
(*barPos = Flatten[{y0,ys,Symbol["y"<>ToString[barN+1]]}];LL=2\[Pi]/kCeiling[L k/(2 \[Pi])];Total[(curWell=#;\[Chi][curWell, ds, ys](Product[HeavisideTheta[-barPos[[i+1]]+Sign[x]Mod[Abs[x], LL]], {i, 0, curWell-1}]Product[HeavisideTheta[barPos[[i+1]]-Sign[x]Mod[Abs[x], LL]], {i, curWell,barN+1}]))&/@Range[1,barN+1]]*)
HeavisideTheta[x+L/2]HeavisideTheta[-x+L/2]Total[(curWell=#;\[Chi][curWell, ds, ys]Product[HeavisideTheta[-barPos[[i]]+x], {i, 1, curWell-1}]Product[HeavisideTheta[barPos[[i]]-x], {i, curWell,barN}])&/@Range[1,barN+1]])]
onePartWF[x_] =\[CapitalPsi][ds,ys, numBar]//.Flatten[{specificRules, \[Delta]->Delta,\[Gamma]->GGamma, randRules, #}]&/@rootPt;
rootList = (k/.#)&/@rootPt
ntoplot=Range[1,Length[rootPt]](*{numBar,2*numBar}*);

Flatten[{specificRules, \[Delta]->Delta, \[Gamma]->GGamma,randRules, rootPt}]
nrmSymb = calcNorm[ds, ys];
normSqr = (nrmSymb//.Flatten[{\[Delta]->Delta,\[Gamma]->GGamma,{randRules},{specificRules}, {#}}])&/@rootPt;
normSqr[[ntoplot]]
(*NIntegrate[Abs[onePartWF[x][[#]]]^2/1, Evaluate[ ({x, -L/2, L/2}/.specificRules)]]&/@ntoplot*)
TGDensity=(Total@(Abs[onePartWF[x][[#]]]^2/normSqr[[#]]&/@ntoplot)[[1;;#]])&/@ntoplot;
SinglePartDensity = Abs[onePartWF[x][[#]]]^2/normSqr[[#]]&/@ntoplot;
gridlns = Flatten[{randYs[[1;;-1]]/.{\[Delta]->Delta, \[Gamma]->GGamma}, (-L/.specificRules)/2,(L/.specificRules)/2}];
glStyled=Flatten[{{gridlns[[#]], If[Mod[#,2]==1,{(*Dashed*)}, {(*Directive[Thick]*)}]}&/@ Range[1,Length[gridlns]-2],{{gridlns[[-2]], Directive[Black, Thick]}, {gridlns[[-1]], Directive[Black, Thick]}}},1]
legnd=Evaluate[ToString[#[[1]]]&/@rootPt]
SetOptions[Plot,BaseStyle->{FontFamily->"Times",FontSize->16, FontColor->Black}];
TGdensPlot = Plot[Evaluate[TGDensity[[2]](*[[{6,7}]]*)],Evaluate[ ({x, -L/1.8, L/1.8}/.specificRules)],GridLines -> {glStyled, {}}, PlotRange->Full(*, PlotStyle\[Rule]{Blue}*), Evaluated->True, PlotLegends->{"n=7", "n=14"}(*legnd[[ntoplot]]*), PlotPoints->150, AxesLabel->{x, ""}, PlotStyle->{Thickness[0.005], Thickness[0.002]}, Filling->Bottom, Ticks->Automatic, Exclusions->None]
tgplotname = StringJoin["Pictures/OIST/sze/TGdensPlot_TG_nbar_", ToString[numBar],"_gamma_", ToString[N[GGamma, 2]], ".png"];
Export[tgplotname, TGdensPlot];


(*TGLimit WF, 2p, simple symmetrization of WF*)
intRange = 5; (*should range from -L/2 to + L/2*)
(*WF and norm*)
DoublePartWF = Abs[onePartWF[x1][[1]]*onePartWF[x2][[2]] - onePartWF[x2][[1]]*onePartWF[x1][[2]]];
DPWFNorm = NIntegrate[Abs[DoublePartWF]^2, {x1, -intRange, intRange}, {x2, -intRange, intRange}]
(*plot WF and norm*)
DPWFplot = DensityPlot[2 Abs[(DoublePartWF)]^2/.HeavisideTheta[0]->1, {x1, -intRange, intRange}, {x2, -intRange, intRange}, FrameLabel->{"x1", "x2"}, PlotLegends->Automatic,PlotRange->Full, PlotPoints->50]
plotname = StringJoin["Pictures/OIST/sze/DPWFplot_TG_nbar_", ToString[numBar],"_gamma_", ToString[N[GGamma, 2]], ".png"]
Export[plotname, DPWFplot];
Print["So far so good 2"]
Print[rootPt]
Print["----------"]
Print[DoublePartWF]
Print["----------"]
Print[DPWFNorm]
Print["----------"]
Print[gridlns]
(*Print["----------"]
Print[rootPt]*)


(* ::Title:: *)
(*Tonks-Girardeau gas: the Bethe Ansatz solution in a random finite lattice, analytical form*)


(* ::Input:: *)
(*Import["/Users/irina 1/qsu_repo/bethe_ansatz/RandomLattice.png"]*)


(* ::Chapter:: *)
(*One particle case*)


(* ::Text:: *)
(*Here I define all the parameters of the system.*)


(* ::Input:: *)
(*numBar=2; (* Total number of the delta barriers*)*)
(*numBar = ToExpression[$ScriptCommandLine[[2]]]*)
(*ys = Symbol["y"<>ToString[#]]&/@Range[1,numBar] (*List of the variables (placeholders) for the positions of the barriers*)*)
(*ds=Symbol["d"<>ToString[#]]&/@Range[1,numBar](*List of the variables (placeholders) for the heights of the barriers*)*)
(*specificRules = {L->10, \[HBar]->1, m->1}; (*Rules to fix the length of the box L, and the \[HBar]/m ratio*)*)
(*$Assumptions=Element[ds, Reals]&&Element[ys, Reals]&&Element[L, Reals] &&Element[m, Reals]&& Element[\[HBar], Reals]&&Element[k, Reals](*General assumptions about parameters. They are all real.*)*)


(* ::Text:: *)
(*Here I define cyclic permutation of the barriers with \[CapitalDelta]\[Element][-L, L] shift. Note that for negative \[CapitalDelta] the effective shift is equivalent to  L+\[CapitalDelta] shift. This is precisely what function Mod[\[CapitalDelta], L] is doing.*)


(* ::Input:: *)
(*(* Cyclic permutation of the barrier heights *)*)
(*cyclicBarrierShiftHeights[\[CapitalDelta]_, ys_, ds_, L_]:=Module[{shift, signs, dels, resultDs, Nb},(Nb = Length[ds];dels=Flatten[{L/2-#&/@ Reverse[ys], L}];signs=Sign[Mod[\[CapitalDelta],L] - dels]; shift=Floor[Total[signs + 1]/2];resultDs = Flatten[{ds[[Nb-shift+1;;Nb]],ds[[1;;Nb-shift]]}])]*)
(*(* Cyclic permutation of the barrier positions *)*)
(*cyclicBarrierShiftPos[\[CapitalDelta]_, ys_, L_]:=Module[{Nb, resultYs},(Nb = Length[ds];resultYs = Sort[Mod[ys+\[CapitalDelta] + L/2, L]-L/2])]*)


(* ::Section:: *)
(*The Bethe equations*)


(* ::Text:: *)
(*Here we construct a Bethe Equation *)


(* ::Input:: *)
(*\[ScriptCapitalT][n_, ys_, ds_]:= Module[{combs},(combs=Subsets[Range[1,Length[ys]], {n}];2^(2n) Total[Times@@(ds[[#]] m/\[HBar]^2)If[Length[#]>=1, Sin[k(ys[[#]][[1]]+L/2)]Sin[k(L/2-ys[[#]][[n]])], Sin[k L]]Product[Sin[k(ys[[#]][[j]]-ys[[#]][[j-1]])],{j,2,n}]&/@combs])]*)
(*\[ScriptCapitalT][1, ys, ds];*)


(* ::Input:: *)
(*makeBE[ys_, ds_]:=Sum[\[ScriptCapitalT][n, ys, ds]k^(Length[ds]-n), {n, 0, Length[ds]}]==0*)
(*AbsoluteTiming[TGBE=makeBE[ys, ds];]*)


(* ::Section:: *)
(*The Coefficients*)


(* ::Text:: *)
(*Here I reconstruct the coefficients and the normaliation constant*)


(* ::Input:: *)
(*(*n is the number of the well in which the particle is in, for 1-particle case. The coefficient we are calculating is Subscript["\[ScriptCapitalA]", {1,0,0,0,0}][{k1}], with 1 in the well number n.  ds - heights of the barriers, ys - positions of the barriers*)*)
(*makeCoef[n_, ds_, ys_]:=Module[{combs, posCombs, elems, curProd, curYP},(combs=Subsets[Range[1,n-1], n-1]; elems = Total[(Times@@(ds[[#]]*(-2 m I)/(k \[HBar]^2))If[Length[#]>0,(1-Exp[-I k (L + 2 ys[[#]][[1]])]), 1])Product[(1-Exp[-I 2 k (+ys[[#]][[j]]-ys[[#]][[j-1]])]), {j, 2, Length[#]}]&/@combs])];*)
(*makeCoef[2, ds, ys];*)


(* ::Input:: *)
(*(*Real part of a coefficient*)*)


(* ::Input:: *)
(*realPartMakeCoef[n_, ds_, ys_]:=Module[{combs, posCombs, elems, curProd, curYP},(combs=Subsets[Range[1,n-1], n-1]; elems = Total[(Times@@(ds[[#]]) ((-4 m )/(k \[HBar]^2))^Length[#] If[Length[#]>0,(-Sin[k(L/2+ys[[#]][[1]])]Cos[k(ys[[#]][[Length[#]]]+L/2)]), 1])Product[(Sin[ k (-ys[[#]][[j]]+ys[[#]][[j-1]])]), {j, 2, Length[#]}]&/@combs])];*)
(*realPartMakeCoef[2, ds, ys];*)


(* ::Input:: *)
(*(*Imaginary part of a coefficient *)*)


(* ::Input:: *)
(*imPartMakeCoef[n_, ds_, ys_]:=Module[{combs, posCombs, elems, curProd, curYP},(combs=Subsets[Range[1,n-1], {1,n-1}]; elems = Total[(Times@@(ds[[#]]) ((-4 m )/(k \[HBar]^2))^Length[#] If[Length[#]>0,(Sin[k(L/2+ys[[#]][[1]])]Sin[k(ys[[#]][[Length[#]]]+L/2)]), 1])Product[(Sin[ k (-ys[[#]][[j]]+ys[[#]][[j-1]])]), {j, 2, Length[#]}]&/@combs])];*)
(*imPartMakeCoef[2, {d1,d2,d3}, {y1,y2,y3}];*)


(* ::Input:: *)
(*(* Calculate the normalization constant*)*)


(* ::Input:: *)
(*calcNorm[ds_, ys_]:= Module[{ReA, ImA, ReA1, ImA1},(ReA1 = realPartMakeCoef[numBar+1, ds, ys]; ImA1=imPartMakeCoef[numBar+1, ds, ys];2 (ys[[1]]+L/2 + -(Sin[k(L+2ys[[1]])])/(2 k)) + 2 Total[(ReA =realPartMakeCoef[#, ds, ys];ImA=imPartMakeCoef[#, ds, ys]; ReA^2 (ys[[#]]-ys[[#-1]] + -(Sin[k(L+2ys[[#]])]-Sin[k (L+2ys[[#-1]])])/(2 k)) + ImA^2 (ys[[#]]-ys[[#-1]] - -(Sin[k(L+2ys[[#]])]-Sin[k (L+2ys[[#-1]])])/(2 k)) - ReA ImA (Cos[k(L+2ys[[#]])]-Cos[k (L+2ys[[#-1]])])/k)&/@Range[2, numBar]]+2(ReA1^2 (L/2-ys[[numBar]] + -(Sin[k 2 L]-Sin[k (L+2ys[[numBar]])])/(2 k)) + ImA1^2 (L/2-ys[[numBar]] - -(Sin[k 2 L]-Sin[k (L+2ys[[numBar]])])/(2 k)) - ReA1 ImA1 (Cos[k 2 L]-Cos[k (L+2ys[[numBar]])])/k))]*)


(* ::Chapter:: *)
(*Actual solution*)


(* ::Input:: *)
(*a = \[Gamma] * L/(numBar-1);*)
(*y0  = -L/2+(L-(numBar-1)a)/2;randYsAllDelta=(y0+a*(#) + \[Delta] (y0+L/2))/.specificRules&/@Range[0, numBar-1]//FullSimplify;*)
(*(*randYsAllDelta=((-L/2 + L/(2(numBar)))+L/numBar(#) + \[Delta]L/(2(numBar))(**Mod[#,2]*))/.specificRules&/@Range[0, numBar-1];*)*)
(*thumorse={0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0};*)
(*randDs =(100+0*Mod[#+0,2])&/@Range[1, numBar](* RandomReal[{0, 3}, Length[ds]]*)(*(1)&/@Range[1, numBar]*) ;*)
(*(*randDs[[{-7, -6}]] = 0;*)*)
(*randRulesAllDelta= Flatten[{ys[[#]]->randYsAllDelta[[#]],ds[[#]]->randDs[[#]]}&/@Range[1, Length[ds]]]*)
(*allYsAllDelta=ys[[#]]->randYsAllDelta[[#]]&/@Range[1, Length[ds]];*)


(* ::Input:: *)
(*Plot[randYsAllDelta/.Flatten[{specificRules, \[Delta]->0}], {\[Gamma], 0, 3}]*)


(* ::Input:: *)
(*(*cyclicBarrierShiftHeights[ -2, randYsAllDelta, ds, L/.specificRules];*)
(*Plot[cyclicBarrierShiftPos[\[CapitalDelta], randYsAllDelta, L/.specificRules][[1]], {\[CapitalDelta], -5,5}];*)*)


(* ::Input:: *)
(*randYsAllDelta[[2]]*)


(* ::Input:: *)
(*Solve[randYsAllDelta[[2]]== 10/6, \[Gamma]]*)


(* ::Input:: *)
(*bEF = Experimental`CreateNumericalFunction[{k,g},{TGBE[[1]]/.Flatten[{randRulesAllDelta, specificRules}]/.{\[Delta]->0, \[Gamma]->g}}, {1}, Compiled->True];*)


(* ::Input:: *)
(*bEF[{1,0.5}]*)
(*bEF["CompiledFunction"][1,0.5]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*betheEquationsF[k_,g_]=TGBE[[1]]/.Flatten[{randRulesAllDelta, specificRules}]/.{\[Delta]->0, \[Gamma]->g};*)
(*bECompiled = Compile[{k, g}, betheEquationsF[k,g]]*)
(*rationalGammas = Range[0,3, 1/100];*)
(*momenta = Range[1,10,0.01];*)
(*(*butterflyValues = Flatten[Module[{g, kk},(g = rationalGammas[[#]];(kk=momenta[[#]];{g,kk,If[Abs[bEF["CompiledFunction"][kk, g][[1]]]<1^(-2), 1, 0]})&/@Range[1, Length[momenta]])&/@Range[1, Length[rationalGammas]]],1];*)
(*ListDensityPlot[butterflyValues, PlotLegends\[Rule]Automatic, InterpolationOrder\[Rule]0]*)*)
(**)


(* ::Input:: *)
(*betheEquations=TGBE[[1]]/.Flatten[{randRulesAllDelta, specificRules}]/.{\[Delta]->\[CapitalDelta],\[Gamma]->(numBar-1)/numBar};*)
(*SetOptions[ContourPlot,BaseStyle->{FontFamily->"Times",FontSize->14}];*)
(*plotObj1=ContourPlot[{betheEquations==0},{\[CapitalDelta], -1,1}, {k, 1, 10}, PlotRange->Full, PlotLegends->Automatic, PlotPoints->100, FrameLabel->{"\[CapitalDelta]", "k"}, AspectRatio->2 , Epilog->{Thin,{Green,Point[{-0.7,8.68}],Green,Point[{-0.7,4.38}], Black, Point[{-0.4,8.72}], Black,Point[{-0.4,4.95}]}}]*)
(*(*orderedPairs=First[Cases[plotObj,x_GraphicsComplex \[RuleDelayed]  First@x,Infinity]];*)
(*orderedSqrP={#[[1]], #[[2]]^2} &/@ orderedPairs;*)
(*lineData=Cases[plotObj,x_Line \[RuleDelayed]  x,Infinity];*)
(*linesIndex = Cases[lineData,Line[x_] \[RuleDelayed]  x,Infinity];*)
(*allLines = orderedSqrP[[#]]&/@ linesIndex;*)
(*Length[orderedPairs]*)
(*Length[lineData]*)*)


(* ::Input:: *)
(*Delta=0;*)
(*GGamma = 1/3(*(numBar-1)/numBar*);*)
(*GGamma = 0;*)
(*randYs=randYsAllDelta/.Flatten[{specificRules}](*Sort[RandomReal[{-L/2, L/2}/.specificRules, Length[ys]]]*)(*{3}*);*)
(*randRules =randRulesAllDelta/.Flatten[{specificRules}]*)
(**)


(* ::Input:: *)
(*(*plotObj=ContourPlot[Evaluate[{TGBE/.k\[Rule]k1, TGBE/.k\[Rule]k2}//.Flatten[{specificRules, \[Delta]\[Rule]Delta, \[Gamma]\[Rule]GGamma, randRules}]], {k1, -0.1, 9}, {k2, -0.1, 5}, PlotPoints\[Rule]300,PlotLegends\[Rule]Automatic, FrameLabel\[Rule]{"k1", "k2"}]*)*)


(* ::Input:: *)
(*(*sampleroots=Cases[(*Normal@*)plotObj,(*Line[pts_,___]\[RuleDelayed]pts*)GraphicsComplex[pts__]\[RuleDelayed]pts,Infinity];*)
(*ListPlot[#[[1]]&/@sampleroots[[1]]]*)*)


(* ::Input:: *)
(*initPt ={3.9};*)
(*rootPt = FindRoot[TGBE//.Flatten[{specificRules, \[Delta]->Delta,\[Gamma]->GGamma, randRules}],{k, #}]&/@initPt *)


(* ::Input:: *)
(*rootPt=NSolve[(TGBE//.Flatten[{specificRules, \[Delta]->Delta,\[Gamma]->GGamma, randRules}])&&k > 0 && k <4,k, Reals, WorkingPrecision->20]*)


(* ::Input:: *)
(*(TGBE[[1]]//.Flatten[{specificRules, \[Delta]->Delta,\[Gamma]->GGamma, randRules, rootPt[[#]]}])&/@Range[1,Length[rootPt]]*)


(* ::Section:: *)
(*The wave function*)


(* ::Input:: *)
(*Acoeff[block_, P_]=Subscript["\[ScriptCapitalA]", block][P];*)
(*allBlocks[N_, barN_] := DeleteDuplicates[Flatten[Permutations[#]&/@(PadRight[#, barN+1]&/@IntegerPartitions[N, barN+1]),1]]*)
(*\[Chi][wellN_, ds_, ys_] := Module[{ Ep}, (Ep ={k,-k}; \!\( *)
(*\*UnderoverscriptBox[\(\[Sum]\), \(ne = 1\), \(\(Dimensions[Ep]\)[\([1]\)]\)]\( *)
(*\*SuperscriptBox[\((\(-Exp[\(-I\)\ k\ L]\))\), \(If[Refine[Sign[Ep[\([ne]\)]], k > 0] == 1, 0, 1]\)] \((makeCoef[\ wellN, ds, \ ys] /. k -> Ep[\([ne]\)])\) Exp[I\ Ep[\([ne]\)]\ x]\)\))]*)
(*\[CapitalPsi][ds_,ys_, barN_] := Module[{bN, allB, N, dels, curWell, barHs, barPos,LL}, (N=1;allB=allBlocks[N, barN];*)
(*barPos=ys;*)
(*(*barPos = Flatten[{y0,ys,Symbol["y"<>ToString[barN+1]]}];LL=2\[Pi]/kCeiling[L k/(2 \[Pi])];Total[(curWell=#;\[Chi][curWell, ds, ys](Product[HeavisideTheta[-barPos[[i+1]]+Sign[x]Mod[Abs[x], LL]], {i, 0, curWell-1}]Product[HeavisideTheta[barPos[[i+1]]-Sign[x]Mod[Abs[x], LL]], {i, curWell,barN+1}]))&/@Range[1,barN+1]]*)*)
(*HeavisideTheta[x+L/2]HeavisideTheta[-x+L/2]Total[(curWell=#;\[Chi][curWell, ds, ys]Product[HeavisideTheta[-barPos[[i]]+x], {i, 1, curWell-1}]Product[HeavisideTheta[barPos[[i]]-x], {i, curWell,barN}])&/@Range[1,barN+1]])]*)


(* ::Input:: *)
(*onePartWF[x_] =\[CapitalPsi][ds,ys, numBar]//.Flatten[{specificRules, \[Delta]->Delta,\[Gamma]->GGamma, randRules, #}]&/@rootPt;*)


(* ::Input:: *)
(*rootList = (k/.#)&/@rootPt*)


(* ::Input:: *)
(*ntoplot=Range[1,Length[rootPt]](*{numBar,2*numBar}*);*)
(**)
(*Flatten[{specificRules, \[Delta]->Delta, \[Gamma]->GGamma,randRules, rootPt}]*)
(*nrmSymb = calcNorm[ds, ys];*)
(*normSqr = (nrmSymb//.Flatten[{\[Delta]->Delta,\[Gamma]->GGamma,{randRules},{specificRules}, {#}}])&/@rootPt;*)
(*normSqr[[ntoplot]]*)
(*(*NIntegrate[Abs[onePartWF[x][[#]]]^2/1, Evaluate[ ({x, -L/2, L/2}/.specificRules)]]&/@ntoplot*)*)


(* ::Input:: *)
(*TGDensity=(Total@(Abs[onePartWF[x][[#]]]^2/normSqr[[#]]&/@ntoplot)[[1;;#]])&/@ntoplot;*)


(* ::Input:: *)
(*SinglePartDensity = Abs[onePartWF[x][[#]]]^2/normSqr[[#]]&/@ntoplot;*)


(* ::Input:: *)
(*gridlns = Flatten[{randYs[[1;;-1]]/.{\[Delta]->Delta, \[Gamma]->GGamma}, (-L/.specificRules)/2,(L/.specificRules)/2}];*)
(*glStyled=Flatten[{{gridlns[[#]], If[Mod[#,2]==1,{(*Dashed*)}, {(*Directive[Thick]*)}]}&/@ Range[1,Length[gridlns]-2],{{gridlns[[-2]], Directive[Black, Thick]}, {gridlns[[-1]], Directive[Black, Thick]}}},1]*)
(*legnd=Evaluate[ToString[#[[1]]]&/@rootPt]*)
(*SetOptions[Plot,BaseStyle->{FontFamily->"Times",FontSize->16, FontColor->Black}];*)
(*Plot[Evaluate[TGDensity[[2]](*[[{6,7}]]*)],Evaluate[ ({x, -L/1.8, L/1.8}/.specificRules)],GridLines -> {glStyled, {}}, PlotRange->Full(*, PlotStyle\[Rule]{Blue}*), Evaluated->True, PlotLegends->{"n=7", "n=14"}(*legnd[[ntoplot]]*), PlotPoints->150, AxesLabel->{x, ""}, PlotStyle->{Thickness[0.005], Thickness[0.002]}, Filling->Bottom, Ticks->Automatic, Exclusions->None]*)


(* ::Input:: *)
(**)


(* ::Subsection:: *)
(*TGLimit WF, 2p, simple symmetrization of WF*)


(* ::Input:: *)
(*intRange = 5; (*should range from -L/2 to + L/2*)*)


(* ::Input:: *)
(*(*WF and norm*)*)
(*DoublePartWF = onePartWF[x1][[1]]*onePartWF[x2][[2]] + onePartWF[x2][[1]]*onePartWF[x1][[2]];*)


(* ::Input:: *)
(*DPWFNorm = NIntegrate[Abs[DoublePartWF]^2, {x1, -intRange, intRange}, {x2, -intRange, intRange}]*)


(* ::Input:: *)
(*(*plot WF and norm*)*)
(*DensityPlot[2 Abs[(DoublePartWF)]^2/.HeavisideTheta[0]->1, {x1, -intRange, intRange}, {x2, -intRange, intRange}, FrameLabel->{"x1", "x2"}, PlotLegends->Automatic,PlotRange->Full, PlotPoints->50]*)


(* ::Chapter:: *)
(*Other stuff *)


(* ::Subsection:: *)
(*3 Particle WF*)


(* ::Input:: *)
(*(*TripplePartWF = onePartWF[x1][[1]]*onePartWF[x2][[2]] + onePartWF[x2][[1]]*onePartWF[x1][[2]];*)*)


(* ::Input:: *)
(*expr = Det[{{A, b, c}, {F, G, H}, {J, d, k}}]*)
(*expr = ToString[expr]*)
(*expr = StringReplace[ expr, "-" -> "+"]*)
(*expr = ToExpression[expr]*)


(* ::Input:: *)
(*testMatrix = {{A, b, c}, {F, G, H}, {J, d, k}}*)


(* ::Input:: *)
(*testMatrix[[1,1]]*)


(* ::Input:: *)
(*myS = "this is string "*)
(*StringReplace[myS, "string" -> "elephant"]*)


(* ::Input:: *)
(*DensityPlot3D[2 Abs[(TriplePartWF)]^2/.HeavisideTheta[0]->1, {x1, -intRange, intRange}, {x2, -intRange, intRange},{x3, -intRange, intRange }, FrameLabel->{"x1", "x2","x3"}, PlotLegends->Automatic,PlotRange->Full, PlotPoints->50]*)


(* ::Subsection:: *)
(*Probability calculations*)


(* ::Input:: *)
(*gridlns*)


(* ::Input:: *)
(*barrierPairs = {{-5, -5/3}, {-5/3, 5/3}, {5/3, 5}}*)
(*ntoprob = Range[1, Length[barrierPairs]]*)


(* ::Input:: *)
(*(*2p WF*)*)


(* ::Input:: *)
(*(*TGDensity*)*)
(*nprob = (numBar+1)^2;*)
(*ngap = numBar +1;*)


(* ::Input:: *)
(*pij[i_, j_]:= NIntegrate[Abs[DoublePartWF]^2, Flatten[{x1,barrierPairs[[i]]}], Flatten[{x2,barrierPairs[[j]]}],  WorkingPrecision-> 10]/DPWFNorm*)


(* ::Input:: *)
(*probSingle= NIntegrate[Abs[DoublePartWF]^2, Flatten[{x1,barrierPairs[[1]]}], Flatten[{x2,barrierPairs[[1]]}]]/DPWFNorm*)


(* ::Input:: *)
(* *)


(* ::Subchapter:: *)
(*Determining Szilard work output*)


(* ::Input:: *)
(*fQSEwork[p_, pout_] =p*Log[p/ pout]*)


(* ::Input:: *)
(*pij[1,1]*)


(* ::InheritFromParent:: *)
(*0.32340812731684887`*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*QSEwork = 0 ;*)
(*q = 0;*)
(*For[i =1, i<=   ngap,i++,*)
(*For[j= 1, j<= ngap, j++,*)
(*(*syntax: If[condition, case true, case false]*)*)
(*If[i==j, out = 1, out = 1/2];*)
(*Print[i, j];*)
(*Print[out];*)
(*Print[pij[i,j]];*)
(*QSEwork= QSEwork + fQSEwork[pij[i,j], out];*)
(*q = q +pij[i,j];*)
(*]*)
(*]*)
(**)
(*Print["QSEwork"]*)
(*Print[QSEwork]*)
(*Print["-----"]*)
(*Print["Rel QSEwork"]*)
(*QSEwork/N[Log[1/2]]*)
(*Print["-----"]*)
(*Print["total prob"]*)
(*Print[q]*)


(* ::Input:: *)
(**)
