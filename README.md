# SZE
Code written for internship on the Quantum Szilard Engine
The curious physicist's guide to Lisa's Szilard code

20/09/2018
If questions and confusions persist, email me at lisa.mickel@gmail.com


##
PYTHON

makeArrays: make arrays with all possible matches for slurm 

readSango: read sango output (work) and plot [code not so nice]

sze_2p and sze_2p mod: runs mathematica scripts with Irina's code to generate WF and generate probability arrays by integrating over the WF. Generated WF and probArrays are saved in txt files and needn't be regenerated each time. Different plotting options available; beauty of the code was reduced due to tryouts / looking at the change in k values that were manually imported to the script
needs: 
mathematica/altered_code/RandomLatticeTG_ShortMoveOut_mod2[...].wl
mathematica/my_code/probList.wl


sze_entropy[] : Use work/ probArrays to see how observables behave during the cycle; I would not trust it...


(**)
MATHEMATICA

my_code:

executable_template: template to create executable mathematica code; please note that code needs to be made executable from terminal with the 'chmod' command

lala: examples for things that can be used in other scripts

probList: created probability array (2d) from wavefunction saved by mathematica/altered_code/RandomLatticeTG_ShortMoveOut_mod2[...].wl
probList_d: barrier height can also be given as input argument

quickCalc: not important, I used it to plot graphs and try mathematica things

k_WF: don't remember, looks like a shortened version of Irina's code...

altered_code:

versions of Irina's code, either shortened or made executable.
mathematica/altered_code/RandomLatticeTG_ShortMoveOut_mod2[...].wl
--> mod2: WF not in TG limit
--> mod2_realTG : WF in TG limit
--> mod2_realTG_d: height of barriers can be given as input argument


%%
MATLAB

mossy_code: ask Mossy for details, contains the magic code which was altered to sparce_magic_nd_exec[] as found in lisa_code


lisa_code:
note: alot of code, not all works or is relevant, please treat with care!

important code:
thermalSZE_InteractingBosons_estates_slurm
script to be run from sango, gets energy levels for interacting bosons  in infinite well(variable interaction strength and barrier number, which influences the box size) from magic code (numerical diagonalization) and saves them in a) output file b) whole workspace that can be loaded locally to obtain work output of SZE

thermalSZE_Bosons
W(T, nb) for non-interacting bosons, analytical code

thermalSZE_int_loadws
W(T, nb, g) --> interactions possible, energy states for scenarios with particles in same well are loaded from matlab workspace

checkSymmetry2D
checks symmetry of WF array as produced by magic code to eliminate fermionic states

getElevelsDP
returns array with combinations of 2particles over energy levels in the indistinguishable case (i.e. (1,1), (1,2), (1,3), (2,2), (2,3) etc.)

getDistElevelsDP
returns array with combinations of 2particles over energy levels in the distinguishable case (i.e. (1,1), (1,2), (2,1), (1,3), (3,1), (2,2), (2,3), (3,2) etc.)

getDiffA
returns partition function and energy states obtained analytically from elevels passed to function

getZ
calculated the partition function from given energies (i.e. 2p energy state of a system) and beta 
