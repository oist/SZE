% thermal SZE - interacting bosons

clear all;

%% Settings
npart = 2; %thou shalt not change the particle number
hb = 1;
m  = 1;
L = 10;
Escale = getEnInf(1, L, hb, m);
maxlevel = 1000;

%%
EstatesI = getElevelsDP(maxlevel, 1);
EstatesD = getDistElevelsDP(maxlevel);
 
%% nstates
nstates = 40000;
nstates = (length(EstatesI));
%nstates = min(length(Estates), nstates);
%% check sufficiency of nstates!
TAr0 = [0.01, 0.1, 1, 10, 100, 1000,1e4, 1e5];%, 1e6,1e7, 1e8, 1e9];%, 1e10, 1e11, 1e12, 1e13,1e14,  1e15, 1e16, 1e17, 1e18, 1e19, 1e20];%1e11, 1e15, 1e20];

npart =2;

%non interacting bosons
Zsame = zeros(1, length(TAr0));
Zdiff = zeros(1, length(TAr0));
for b=1:length(TAr0)
    beta = 1/TAr0(b);
    [Zsame(b), Esame(b)] = getDiffA(beta, L/2, maxlevel, EstatesI, nstates, hb, m);
    [Zdiff(b), Ediff(b)] = getDiffA(beta, L/2, maxlevel, EstatesD, nstates, hb, m);
end
Ztot = 2*Zsame + Zdiff;
prem = Zdiff./Ztot;

figure(1)
plot(TAr0, Zsame./Ztot);
title('p0');
set(gca, 'XScale', 'log');

%% vary T, nb, g
TAr = [0.01,0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50 ,100, 200];
nbAr = [1,2,3, 4, 5, 6, 7, 8, 9, 10];
nstatAr = nstates*ones(1, length(TAr)); % # of states used for Z can vary with T


%% pnrem 
Zsame = zeros(1, length(TAr));
Zdiff = zeros(1, length(TAr));
for b=1:length(TAr)
    beta = 1/TAr(b);
    [Zsame(b), Esame] = getDiffA(beta, L/2, maxlevel, EstatesI, nstates, hb, m);
    [Zdiff(b), Ediff] = getDiffA(beta, L/2, maxlevel, EstatesD, nstates, hb, m);
end
Ztot = 2*Zsame + Zdiff;
prem = Zdiff./Ztot;
%% pn

eS=zeros(length(nbAr), nstates);
eD=zeros(length(nbAr), nstates);
Zs=zeros(length(nbAr), length(TAr));
Zd=zeros(length(nbAr), length(TAr));
Ztotal=zeros(length(nbAr), length(TAr));

for k = 1:length(nbAr)
    nb = nbAr(k);
    l = L/(nb+1);
    eS(k,:) = getEnergiesA(l, maxlevel, EstatesI, nstates, hb, m);
    eD(k,:) = getEnergiesA(l, maxlevel, EstatesD, nstates, hb, m);
    for b=1:length(TAr)
        beta = 1/TAr(b);
        Zs(k,b) = getZ(eS(k,:), beta);
        Zd(k,b) = getZ(eD(k,:), beta);
    end
    
end
%%
ns = (nbAr+1);
nd = (nbAr+1).*nbAr/2;
Ztotal = ns.'.*Zs + nd.'.*Zd;
ps = Zs./Ztotal;
pd = Zd./Ztotal;

betaW = - (ns.'.*ps.*log(ps) + nd.'.*pd.*log(pd./prem));

%% you will now enter the reign of plotting
%%
figure(2)
colormap jet
surf(betaW); hold on;
imagesc(betaW);
xlabel('T 10^x');
ylabel('nbar');
title('W(T, nbar)/T');

%%
figure(5)
%for k=1:6%length(nbAr)
%    plot(TAr, betaW(k,:)/log(k+1)); hold on;
%end
nbplot = [1,2,4,8];
for k = 1:length(nbplot)
    n = nbplot(k);
    %analytical W T0
    plot(TAr(1:8), getT0LimW(n)/log(n+1)*ones(length(TAr(1:8))), 'k--'); hold on;
    %analytical W TInf
    plot([TAr(end-3:end), 300], getTInfLimW(n)/log(n+1)*ones(length(TAr(end-3:end))+1), 'k--'); hold on;
    
    plot(TAr, betaW(n,:)/log(n+1)); hold on;
end
xlabel('T 10^x');
ylabel('W/ T /log(nb+1)');
title('W(T) varying nb: 1,2, 4, 8');
set(gca, 'XScale', 'log');
%legend('nb 1','nb 2', 'nb 4', 'nb 8');

%%
figure(3)
for k=2:2%length(nbAr)
    %plot(TAr, betaW(k,:)/sum(betaW(k,:))); hold on;
    plot(TAr, ps(k,:)); hold on;
    plot(TAr, pd(k,:)); hold on;
end
plot(TAr, prem);
legend('p_{same}','p_{diff}', 'p_r');
ylabel('probability');
xlabel('T [10^x]');
title('p(T), nb = 2')
set(gca, 'XScale', 'log');
%% 
figure(4)
for k=1:5%length(nbAr
    nb = nbAr(k);
    %analytical ps and pd T0
    plot(TAr(1:9), 2/(nb+1)/(nb+2)*ones(length(TAr(1:9))), 'k--'); hold on;
    %analytical ps TInf)
    plot([TAr(end-3:end), 300], 1/(nb+1)/(nb+1)*ones(length(TAr(end-3:end))+1), 'k--'); hold on;
    %analytical pd TInf)
    plot([TAr(end-3:end),300], 2/(nb+1)/(nb+1)*ones(length(TAr(end-3:end))+1), 'k--'); hold on;
    
    plot(TAr, ps(k,:)); hold on;
    plot(TAr, pd(k,:)); hold on;
end
%legend('nb 1','nb 2', 'nb 3','nb 4');
ylabel('probability');
xlabel('T [10^x]');
title('p_{same}(T) and p_{diff}(T), varying nb')
set(gca, 'XScale', 'log');
%plot(TAr, prem);

%%
figure(6)
tInf = zeros(1,length(nbAr));
t0 = zeros(1,length(nbAr));
for k = 1:length(nbAr)
    nb = nbAr(k);
    tInf(k) = getTInfLimW(nb)/log(nb+1);
    t0(k) = getT0LimW(nb)/log(nb+1);
    
end
%for b=1:length(TAr)
plot(nbAr, betaW(:,1)./log(nbAr+1).', '-o'); hold on;
plot(nbAr, betaW(:,7)./log(nbAr+1).', '-o'); hold on;
plot(nbAr, betaW(:,end)./log(nbAr+1).', '-o'); hold on;
%plot(nbAr,tInf, '-o'); hold on;
%plot(nbAr, t0, '-o'); hold on;
%end
legend('T= 0.01','T= 1','T= 100')%, 'limTInf', 'limT0');
ylabel('W/ T /log(nb+1)');
xlabel('nb');
title('W(nb) varying T')
%plot(nbAr, betaW(:,end))
