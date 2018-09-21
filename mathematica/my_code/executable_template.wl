(* ::Package:: *)

#!/usr/local/bin/MathematicaScript -script
(*Executable mathematica script*)
(*can be run from terminal with e.g. (2 \[Rule] $ScriptCommandLine[[2]])
wolframscript -file executable_template.wl 2
*)


(*simple print*)
Print["Hello World!"]
(*take argument from terminal command line*)
Print[ToExpression[$ScriptCommandLine[[2]]]^2]

(*save plot to file*)
testPlot = Plot[Sin[x], {x, 0, 10}];
Export["test.jpg", testPlot];


(*to import files:
myFile = Import[filename];
*)
