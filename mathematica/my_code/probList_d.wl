(* ::Package:: *)

#!/usr/local/bin/MathematicaScript -script

Print["creating prob array"]
fname = $ScriptCommandLine[[2]]
text = Import[fname];
text2 = StringSplit[text, "----------"]; 
t = StringReplace[fname, ".txt" -> ""]
textname =  StringSplit[Last[StringSplit[t, "/"]], "_"];

Print[t]
Print[textname]
Print[text2]
Print[text2[[2]]]


(*define parameters from text file*)
DoublePartWF = ToExpression[text2[[2]]]
barrierPositions = ToExpression[text2[[4]]]
DPWFNorm = ToExpression[text2[[3]]]
Print[barrierPositions]
gamma = ToExpression[textname[[7]]];
barHeight = ToExpression[textname[[9]]];

Print[barrierPositions]


(*obtain integration intervals*)
CreateIntervals[list_]:=
(
Clear[l];
l = {};
 For[i = 1, i<Length[list],i ++,
AppendTo[l, {list[[i]], list[[i+1]]}]

]
l
)
Print[barrierPositions]

CreateIntervals[Sort[barrierPositions]]
barrierPairs = l
Print["barrierPairs"]
Print[barrierPairs]
For[i=1, i<= Length[barrierPairs], i++,
If[barrierPairs[[i,1]] == barrierPairs[[i,2]], barrierPairs =  Delete[barrierPairs, i]]
]
Print[barrierPairs]


(*define and calculate probability*)
pij[i_, j_]:= NIntegrate[Abs[DoublePartWF]^2, Flatten[{x1,barrierPairs[[i]]}], Flatten[{x2,barrierPairs[[j]]}],  WorkingPrecision-> 10]/DPWFNorm
Print["defined"]
probSingle= NIntegrate[Abs[DoublePartWF]^2, Flatten[{x1,barrierPairs[[1]]}], Flatten[{x2,barrierPairs[[1]]}]]/DPWFNorm;
Print["probSingle"];
Print[probSingle];


getProbArray[ngap_]:=
(
Clear[l];
l = {};
For[i =1, i<=   ngap,i++,
Clear[subl];
subl = {};
For[j= 1, j<= ngap, j++,
(*syntax: If[condition, case true, case false]*)
Print[i, j];
prob = pij[i,j];
Print[prob];
AppendTo[subl, prob];
]
AppendTo[l, subl];
]

l;
)

ngap = Length[barrierPairs];
Print["gamma"]
Print[gamma]
Print["ngap"]
Print[ngap]
Print["barHeight"]
Print[barHeight]
getProbArray[ngap];
probList = l;


Print["probList"];
Print[probList];

reliefplot= ReliefPlot[probList];
relplotname = StringJoin["Pictures/OIST/sze/probListPlot_nbar_", ToString[ngap-1],"_gamma_", ToString[N[gamma, 2]],"_height_",ToString[barHeight], ".png"];
Export[relplotname, reliefplot];
(*visual representation of prob list*)







