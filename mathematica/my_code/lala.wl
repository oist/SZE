(* ::Package:: *)

#!/usr/local/bin/MathematicaScript -script
(*Playscript to try fancy functions and stuffs*)
Print["Hello World!"]
Plot[x^2, {x, 0, 2}]
(*Print[ToExpression[$ScriptCommandLine[[2]]]^2]*)

(*
popup[gr_] := CreateDocument[gr, "CellInsertionPointCell" -> Cell[], 
                             ShowCellBracket -> False, WindowElements -> {}, 
                             WindowFrame -> "Generic", WindowSize -> All,
                             WindowTitle -> None, WindowToolbars -> {}]
      
Plot[Sin[x], {x, 0, 2 \[Pi]}, DisplayFunction -> popup]*)

(*test2 = Plot[Sin[x], {x, 0, 10}];
Export["test2.jpg", test2];*)



(*Import parameters*)
text = Import[$ScriptCommandLine[[2]]];
text2 = StringSplit[text, "----------"];


(*define parameters from text file*)
DoublePartWF = ToExpression[text2[[2]]]
barrierPositions = ToExpression[text2[[4]]]
DPWFNorm = ToExpression[text2[[3]]]


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

CreateIntervals[Sort[barrierPositions]]
barrierPairs = l
Print["barrierPairs"]
Print[barrierPairs]

(*define and calculate probability*)
pij[i_, j_]:= NIntegrate[Abs[DoublePartWF]^2, Flatten[{x1,barrierPairs[[i]]}], Flatten[{x2,barrierPairs[[j]]}],  WorkingPrecision-> 10]/DPWFNorm
Print["defined"]
probSingle= NIntegrate[Abs[DoublePartWF]^2, Flatten[{x1,barrierPairs[[1]]}], Flatten[{x2,barrierPairs[[1]]}]]/DPWFNorm
Print["probSingle"]
Print[probSingle]


getProbArray[ngap_]:=
(
Clear[l];
l = {};
For[i =1, i<=   ngap,i++,
For[j= 1, j<= ngap, j++,
(*syntax: If[condition, case true, case false]*)
Print[i, j];
prob = pij[i,j];
Print[prob];
AppendTo[l, prob];

]
]

l
)

ngap = Length[barrierPairs]
getProbArray[ngap]
probList = l

Print["probList"]
Print[probList]



