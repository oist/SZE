%work area!
clear all
clf

%% Settings
L = 10;
nx = 2^9;
%% Kraken code

%% wrem
nbar = 1;
Nx= (nbar+1)*nx + nbar;
%% make WF
[pInf1, pInf2, eInf1, eInf2] = sparce_magic_nd_exec(0, Nx, L);
[pSplit1, pSplit2, eS1, eS2] = sparce_magic_nd_exec(1, nx, L);

%% getProb Arrays
figure(4)
pAInf = getPosArray(pInf1, pInf2, 2);
figure(5)
pASplit = getPosArray(pSplit1, pSplit2, 2);

wInf = calcWork(pAInf)
wSplit = calcWork(pASplit)
EInf = eInf1 + eInf2
ESplit= eS1+ eS2

%% gap case

[p1, p2, e1, e2] = sparce_magic_nd_exec(2, 2^9, L);
figure(6)
pA = getPosArray(p1,p2, 2);
w = calcWork(pA)
E = e1+ e2
%% cycle


