% thermal SZE - interacting bosons
%% Settings
npart = 2; %thou shalt not change the particle number
hb = 1;
m  = 1;
L = 10;
Escale = getEnInf(1, L, hb, m);
maxlevel = 10;
nstates = 10;
%%
Estates = getElevelsDP(maxlevel, 1);
 
%% nstates
nstates = min(length(Estates), nstates);
%% vary T, nb, g
%gAr, temperature, nbAr defined as input in slurm

%% pn
nx = 2^9; %cannot be increased for npart = 2; if changed to npart =3: decrease to 6 or 7
%% different wells
%% analytic version
for k = 1:length(nbAr)
        nbar = nbAr(k);
        l = L/(nbar+1);
            %get Elevels
            eSingleA = zeros(1,maxlevel);
            for i = 1:maxlevel
                eSingleA(i) = getEnInf(i, l, hb, m);
            end
            eDa = zeros(1, length(Estates));
            for i = 1:length(Estates)
                eDa(i) = eSingleA(Estates(i,1)) + eSingleA(Estates(i,2));
            end
            eDa = sort(eDa);
            eDa = eDa(1:nstates);
end


%% same wells
for k = 1:length(nbAr)
    nbar = nbAr(k);
    g = gAr*(nbar +1);
    l = L/(nbar+1);
        %get Elevels
        [psiSame, eS] = sparce_magic_nd_exec_all(0, nx, l, 2, nstates, g);
end


save(char(strcat('thermalSZE_InteractingBosons_slurm_workspace_gAr_', string(gAr), '_nbAr_', string(nbAr(1)), '.mat')))

fprintf('eDa_\n');
fprintf('%u\n', eDa);
fprintf('_\n');
fprintf('eS\n');
fprintf('%u\n', eS);
fprintf('_\n');
fprintf('g_%u_\n', gAr);
fprintf('nbar_%u_\n', nbAr);
fprintf('nx_%u_\n', nx);
fprintf('nstates_%u_\n',nstates);