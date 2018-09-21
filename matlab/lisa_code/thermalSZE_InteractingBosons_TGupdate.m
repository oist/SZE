% thermal SZE - interacting bosons

clear all;

%% Settings
npart = 2; %thou shalt not change the particle number
hb = 1;
m  = 1;
L = 10;
Escale = getEnInf(1, L, hb, m);
maxlevel = 500;
%pnrem specific
lminP = 0.01;
npoints = 10;
%%
Estates = getElevelsDP(maxlevel, 1); % bosonic states, same occupanices
EstatesSame = getElevelsDP(maxlevel, 0); %TG states, no same occupancies

%% nstates
nstates = length(EstatesSame);
%% vary T, nb, g
TAr = [0.01, 0.1, 1, 10, 100, 1000,1e4, 1e5, 1e6,1e7, 1e8, 1e9]%, 1e10, 1e11, 1e12, 1e13,1e14,  1e15, 1e16, 1e17, 1e18, 1e19, 1e20];%1e11, 1e15, 1e20];

%TAr = [1e16, 1e17, 1e18, 1e19, 1e20];%1e11, 1e15, 1e20];
nbAr = [1,2,3,4,5, 6, 7, 8, 9, 10]%, 11, 12, 13, 14, 15,16,17,18,19,20, 21, 22, 23,24,25, 26, 27, 28, 29, 30];
%% pnrem 
premMaxSame = zeros(1, length(TAr));
premMaxDiff = zeros(1, length(TAr));

lmin = lminP*L;
lhalf = linspace(lmin, L/2, npoints);
lhalf2 = lhalf + L/2 - lhalf(1);
lvalues = [lhalf, lhalf2(2:end)];

for b =1:length(TAr)
    beta = 1/TAr(b);
    psame = zeros(1, length(lvalues));
    pdiff = zeros(1, length(lvalues));
    Zs =  zeros(1, length(lvalues));
    Es =  zeros(1, length(lvalues));
    Zda =  zeros(1, length(lvalues));
    Eda =  zeros(1, length(lvalues));
    for p=1:length(lvalues)
        l = lvalues(p);
        for n =1:nstates
            thisE = 1/2*(getEnInf(Estates(n,1), L-l,hb, m) + ... 
            getEnInf(Estates(n,2), l, hb, m) + getEnInf(Estates(n,2), L-l, hb, m) + ...
            getEnInf(Estates(n,1), l,hb, m));
            Zda(p) = Zda(p) + exp(-beta*thisE);
            Eda(p) = Eda(p) + thisE*exp(-beta*thisE);
        end
        Eda = Eda./Zda;
        
        %[Zda, Eda] = getDiffA(beta, l, maxlevel, Estates, nstatAr(b), hb, m);
        %[Zd, Ed] = getDiffK(beta, nx, l, maxlevel, Estates, nstatAr(b));
        [Zs(p), Es(p)] = getDiffA(beta, l, maxlevel, EstatesSame, nstates, hb, m);
        
    end

    Zs2 = fliplr(Zs);
    
    figure(1)
    Ztotal = Zs + Zda + Zs2;  
    pnremArS = Zs./Ztotal;
    pnremArD = Zda./Ztotal;
    subplot(2,1,1)
    plot(lvalues, pnremArS); hold on;
    plot(lvalues, pnremArD); hold on;
    xlabel('l');
    ylabel('pnrem');
    subplot(2,1,2)
    plot(lvalues, Es); hold on;
    plot(lvalues, Eda); hold on;
    xlabel('l');
    ylabel('E');
    premMaxSame(b) = max(pnremArS);
    premMaxDiff(b) = max(pnremArD);
end

%% pn
%% different wells
%% analytic version
Zdiffa = zeros(length(TAr), length(nbAr));
Ediffa = zeros(length(TAr), length(nbAr));
for k = 1:length(nbAr)
    nbar = nbAr(k);
    l = L/(nbar+1);
    for b = 1:length(TAr)
        beta= 1/TAr(b);
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
        %get Z
        [Zdiffa(b,k),Ediffa(b,k)] = getZ(eDa, beta);
    end
end

%% same wells

Zsame = zeros(length(TAr), length(nbAr));
Esame = zeros(length(TAr), length(nbAr));
for k = 1:length(nbAr)
    nbar = nbAr(k);
    l = L/(nbar+1);
    for b = 1:length(TAr)
        beta= 1/TAr(b);    
        %get Z
        [Zsame(b,k),Esame(b,k)] = getDiffA(beta, l, maxlevel, EstatesSame, nstates, hb, m);
    end
end


%%
%figure(1)
%plot(eDa); hold on; plot(eD); hold on; plot(eS);


%% ingredients are ready
%boring step by step due to lack of smartness
WAr = zeros(length(TAr), length(nbAr));
%using analytic solution for different gap case


for b=1:length(TAr)
    beta = 1/TAr(b);
    for k = 1:length(nbAr)
        nb = nbAr(k);
        nS = nb+1;
        nD = (nb+1)*nb/2;
        Ztot = nS*Zsame(b,k) + nD*Zdiffa(b,k);
        pnS= Zsame(b,k)/Ztot;
        pnD = Zdiffa(b,k)/Ztot;
        pnremS = premMaxSame(b);
        pnremD = premMaxDiff(b);
        WAr(b,k) = -1/beta*(nS*pnS*log(pnS/pnremS) + nD*pnD*log(pnD/pnremD));
    end
end
return
%{gGAr
%matrix calculation!!!CHECK
nS = ones(length(gAr), length(TAr), length(nbAr)).*(nbAr+1);
nD = ones(length(gAr), length(TAr), length(nbAr)).*(nbAr+1).*nbar./2;
Ztotal = nS.*Zsame + nD.*Zdiff;
pnS= Zsame./Ztotal;
pnD = Zdiff./Ztotal;
%W = nS.*pnS.*log(pnS./pnremS) + nD.*pnD.*log(pnD./pnremD);
W = -(nS.*pnS.*log(pnS) + nD.*pnD.*log(pnD));
%}
%% normalize W: by classical, by T, by escale
Wclass = ones(length(gAr), length(TAr), length(nbAr)).*log(nbAr+1);
t3d = ones(length(gAr), length(TAr), length(nbAr)).*(TAr.');
Wnorm = WAr./Escale./Wclass./t3d %check for dimensions
%%
for j=1:length(gAr)
    figure(j+1)
    surf(Wnorm(j,:,:)); hold on;
    imagesc(Wnorm(j,:,:))
    colormap jet
    colorbar
    title(['g = ', string(gAr(j))])
    xlabel('T')
    ylabel('nb')
end

