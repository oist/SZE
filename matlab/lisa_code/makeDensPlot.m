nbAr = [2];
for nb = nbAr
    densPlot_noTun(nb);
end

%% simple

for nb = nbAr
    [psi1, psi2, e1, e2] = sparce_magic_nd_exec(nb, nx, L, 100, 150);
    densPlot_gen(psi1, psi2, nb, 30);
    pA= getPosArray(psi1, psi2, nb);
end

%% kraken
clear all;
clf;
nbAr = [2];
states = getElevelsDP(10) %around 50 states;
L = 10;
nx = 2^7;
nb  = 2;
Nx= (nb+1)*nx + nb;
nstates = length(states);
psi1Ar = zeros(nstates, Nx);
psi2Ar = zeros(nstates, Nx);
E2Ar = zeros(nstates, 1);
E1Ar = zeros(nstates, 1);
%pAr = zeros(nstates, Nx, Nx);

%{
for nb = nbAr
    [psi1, psi2, e1, e2] = sparce_magic_nd_exec(nb, nx, L, 1, 2);
    densPlot_gen(psi1, psi2, nb, 30);
    pA= getPosArray(psi1, psi2, nb);
end
%}



pAr = zeros(Nx, Nx);

for i=1:length(states)
   disp('here')
   [psi1Ar(i,:),psi2Ar(i,:),E1Ar(i),E2Ar(i)]  = sparce_magic_nd_exec(nb, nx, L, states(i,1),states(i,2)); 
   disp('here')
   b = probArray(psi1Ar(i,:), psi2Ar(i,:));
   pAr = pAr + probArray(psi1Ar(i,:), psi2Ar(i,:));
end

imagesc(pAr);