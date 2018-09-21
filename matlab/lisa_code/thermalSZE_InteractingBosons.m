% thermal SZE - interacting bosons

clear all;

%% Settings
npart = 2; %thou shalt not change the particle number
hb = 1;
m  = 1;
L = 10;
Escale = getEnInf(1, L, hb, m);
maxlevel = 100;
nstates = 300;
%pnrem specific
lminP = 0.01;
npoints = 2;
%%
Estates = getElevelsDP(maxlevel, 1);
 
%% nstates
nstates = min(length(Estates), nstates);
%% vary T, nb, g
TAr = [1, 10];
nbAr = [1,2,3];
gAr = [-1,0, 1];
nstatAr = nstates*ones(1, length(TAr)); % # of states used for Z can vary with T
%% pnrem 
premMaxSame = zeros(length(gAr), length(TAr));
premMaxDiff = zeros(length(gAr), length(TAr));

lmin = lminP*L;
lhalf = linspace(lmin, L/2, npoints);
lhalf2 = lhalf + L/2 - lhalf(1);
lvalues = [lhalf, lhalf2(2:end)];

nx = 2^9;

for b =1:length(TAr)
    beta = 1/TAr(b);
    pnremArS = zeros(length(gAr), length(lvalues));
    pnremArD = zeros(length(gAr), length(lvalues));
    psame = zeros(length(gAr), length(lvalues));
    pdiff = zeros(length(gAr), length(lvalues));
    Zs =  zeros(length(gAr), length(lvalues));
    Es =  zeros(length(gAr), length(lvalues));
    Zs2 =  zeros(length(gAr), length(lvalues));
    Ztotal =  zeros(length(gAr), length(lvalues));
    Zda =  zeros(1, length(lvalues));
    Eda =  zeros(1, length(lvalues));
    for p=1:length(lvalues)
        l = lvalues(p);
        for n =1:length(nstatAr(b))
            thisE = 1/2*(getEnInf(Estates(n,1), L-l,hb, m) + ... 
            getEnInf(Estates(n,2), l, hb, m) + getEnInf(Estates(n,2), L-l, hb, m) + ...
            getEnInf(Estates(n,1), l,hb, m));
            Zda(p) = Zda(p) + exp(-beta*thisE);
            Eda(p) = Eda(p) + thisE*exp(-beta*thisE);
        end
        %[Zda, Eda] = getDiffA(beta, l, maxlevel, Estates, nstatAr(b), hb, m);
        %[Zd, Ed] = getDiffK(beta, nx, l, maxlevel, Estates, nstatAr(b));
        

        for j=1:length(gAr)
            [Zs(j,p), Es(j,p)] = getSame(beta, nx, l, nstatAr(b), gAr(j));
            %[Zs(j,p), Es(j,p)] = getDiffA(beta, l, maxlevel, Estates, nstatAr(b), hb, m);
        end        
    end

    Zs2 = fliplr(Zs);
    
    figure(1)
    for j = 1:length(gAr)
        Ztotal(j,:) = Zs(j,:) + Zda + Zs2(j,:);  
        pnremArS(j,:) = Zs(j,:)./Ztotal(j,:);
        pnremArD(j,:) = Zda./Ztotal(j,:);
        
        plot(lvalues, pnremArS(j,:)); hold on; 
        plot(lvalues, pnremArD(j,:)); hold on;
        xlabel('L');
        ylabel('pnrem(l)');
        premMaxSame(j, b) = max(pnremArS(j,:));
        premMaxDiff(j, b) = max(pnremArD(j,:));
    end
end

%% pn
nx = 2^10; %cannot be increased for npart = 2; if changed to npart =3: decrease to 6 or 7
%% different wells
%% kraken version 
%{
Zdiff = zeros(length(gAr), length(TAr), length(nbAr));
Ediff = zeros(length(gAr), length(TAr), length(nbAr));
for j = 1:length(gAr)
    g = gAr(j);
    for k = 1:length(nbAr)
        nbar = nbAr(k);
        l = L/(nbar+1);
        for b = 1:length(TAr)
            beta= 1/TAr(b);
            %get Elevels
            %eToCompute = max(max(Estates([1:nstatAr(b)],:)));
            [psiDiff, eSingle] = sparce_magic_nd_exec_all(0, nx, l, 1,maxlevel,0);
            eD = zeros(1, length(Estates));
            for i = 1:length(Estates)
                eD(i) = eSingle(Estates(i,1)) + eSingle(Estates(i,2));
            end
            eD = sort(eD);
            eD = eD(1:nstatAr(b));
            %get Z
            [Zdiff(j,k,b),Ediff(j,k,b)] = getZ(eD, beta);
        end
    end
end
%}
%% analytic version
Zdiffa = zeros(length(gAr), length(TAr), length(nbAr));
Ediffa = zeros(length(gAr), length(TAr), length(nbAr));
for j = 1:length(gAr)
    g = gAr(j);
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
            eDa = eDa(1:nstatAr(b));
            %get Z
            [Zdiffa(j,k,b),Ediffa(j,k,b)] = getZ(eDa, beta);
        end
    end
end

%%
%{
subplot(2,1,1)
plot(sort(eD)), hold on;
plot(sort(eDa)); 
subplot(2,1,2)
plot((sort(eDa) -sort(eD))./sort(eDa))
%}
%% same wells

Zsame = zeros(length(gAr), length(TAr), length(nbAr));
Esame = zeros(length(gAr), length(TAr), length(nbAr));
for j = 1:length(gAr)
    g = gAr(j);
    for k = 1:length(nbAr)
        nbar = nbAr(k);
        l = L/(nbar+1);
        for b = 1:length(TAr)
            beta= 1/TAr(b);
            %get Elevels
            [psiSame, eS] = sparce_magic_nd_exec_all(0, nx, l, 2, nstatAr(b), g);
            %get Z
            [Zsame(j,k,b),Esame(j,k,b)] = getZ(eS, beta);
        end
    end
end

%%
%figure(1)
%plot(eDa); hold on; plot(eD); hold on; plot(eS);


%% ingredients are ready
%boring step by step due to lack of smartness
WAr = zeros(length(gAr), length(TAr), length(nbAr));
%using analytic solution for different gap case

for j=1:length(gAr)
    for b=1:length(TAr)
        beta = 1/TAr(b);
        for k = 1:length(nbAr)
            nb = nbAr(k);
            nS = nb+1;
            nD = (nb+1)*nb/2;
            Ztot = nS*Zsame(j,b,k) + nD*Zdiffa(j,b,k);
            pnS= Zsame(j,b,k)/Ztot;
            pnD = Zdiffa(j,b,k)/Ztot;
            pnremS = premMaxSame(j,b);
            pnremD = premMaxDiff(j,b);
            WAr(j,b,k) = -1/beta*(nS*pnS*log(pnS/pnremS) + nD*pnD*log(pnD/pnremD));
        end
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

