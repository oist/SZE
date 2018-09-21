% p0 with new formula
L = 10;
maxlevel = 1000; 
hb = 1;
m =1;
nb=2;
nstat = 500000;
%%
EstatesI = getElevelsDP(maxlevel, 1);
EstatesD = getDistElevelsDP(maxlevel);
%%
TAr = [0.01, 0.1, 1, 10, 100, 1000,1e4];%, 1e5, 1e6,1e7, 1e8, 1e9];%, 1e10, 1e11, 1e12, 1e13,1e14,  1e15, 1e16, 1e17, 1e18, 1e19, 1e20];%1e11, 1e15, 1e20];

npart =2;

%non interacting bosons
Zsame = zeros(1, length(TAr));
Zdiff = zeros(1, length(TAr));
for b=1:length(TAr)
    beta = 1/TAr(b);
    [Zsame(b), Esame(b)] = getDiffA(beta, L/2, maxlevel, EstatesI, nstat, hb, m);
    [Zdiff(b), Ediff(b)] = getDiffA(beta, L/2, maxlevel, EstatesD, nstat, hb, m);
end
Ztot = 2*Zsame + Zdiff;
prem = Zdiff./Ztot;

figure(1)
plot(TAr, Zsame./Ztot);
title('p0');
set(gca, 'XScale', 'log');

%% manual Work

%1 barrier case
l = L/2;
eS = getEnergiesA(l, maxlevel, EstatesI, nstat, hb, m);
eD = getEnergiesA(l, maxlevel, EstatesD, nstat, hb, m);
figure(2);
plot(eS); hold on; plot(eD);
title('energies');

Zs=zeros(1, length(TAr));
Zd=zeros(1, length(TAr));
for b=1:length(TAr)
    beta = 1/TAr(b);
    Zs(b) = getZ(eS, beta);
    Zd(b) = getZ(eD, beta);
end
Ztotal = 2*Zs + Zd;
ps = Zs./Ztotal;
pd = Zd./Ztotal;

figure(3)
plot(TAr, Zs, 'g'); hold on; plot(TAr, Zd); hold on; plot(TAr, Ztotal, 'c');
title('Z');

figure(4)
plot(TAr, ps, 'g'); hold on; plot(TAr, pd); hold on; plot(TAr, 2*ps+pd, 'c');
title('ps and pd');
%% work
 WAr = -(2*ps.*log(ps) + pd.*log(pd./prem));
 figure(5)
 plot(TAr, WAr/log(2));
 title('WT');
 set(gca, 'XScale', 'log');

 
 %% automated try
nb = 2;
l = L/(nb+1);
eS = getEnergiesA(l, maxlevel, EstatesI, nstat, hb, m);
eD = getEnergiesA(l, maxlevel, EstatesD, nstat, hb, m);
figure(2);
plot(eS); hold on; plot(eD);
title('energies');

%%

Zs=zeros(1, length(TAr));
Zd=zeros(1, length(TAr));
for b=1:length(TAr)
    beta = 1/TAr(b);
    Zs(b) = getZ(eS, beta);
    Zd(b) = getZ(eD, beta);
end
Ztotal = 2*Zs + Zd;
ps = Zs./Ztotal;
pd = Zd./Ztotal;

figure(6)
plot(TAr, Zs, 'g'); hold on; plot(TAr, Zd); hold on; plot(TAr, Ztotal, 'c');
title('Z');

figure(7)
plot(TAr, ps, 'g'); hold on; plot(TAr, pd); hold on; plot(TAr, 2*ps+pd, 'c');
title('ps and pd');
%% work
 WAr = -(2*ps.*log(ps) + pd.*log(pd));
 figure(8)
 plot(TAr, WAr); hold on;
 plot(TAr, -(2*ps.*log(ps) + pd.*log(pd./prem))); hold on;
 plot(TAr, prem)
 title('WT');
 set(gca, 'XScale', 'log');
