%% Debugging energy expectation value 
clear all;
L = 10;
hb = 1;
m = 1;
nbar = 1;
nx = 2^8;
beta = 0.1;
%% 1) 1 particle case
%<E>(L) plot

npoints = 50;
count = 1;
zArr = zeros(1,npoints);
ezArr = zeros(1,npoints);
L = npoints;
beta = 1;
eLevels = linspace(1, npoints, npoints);
eLevels = eLevels.';
for i =1:round(L/npoints):L
    [zArr(count),  ezArr(count)] = getPartFkt(eLevels,  50, 1, i, hb,m,beta);
    count= count+1;
end
xzvals = linspace(1,L,length(zArr));
plot(xzvals, zArr, 'b'); %hold on; plot(xzvals, ezArr,'r');

%% second version
npoints = 20;
nstates = 50;
beta = 1;
bArr = [0.01];
zArrman = zeros(1, npoints);
for b = 1: length(bArr)
    beta = bArr(b)
    for i = 1:npoints
        thisZ = 0;
        eStates = zeros(1, nstates);
        pInZ = zeros(1, nstates);
        for j= 1:nstates
            eStates(j) =getEnInf(j, i, hb, m);
            pInZ(j)= exp(-beta*eStates(j));
            thisZ = thisZ +pInZ(j);
        end
        xv = linspace(1, nstates, nstates);
        figure(1)
        subplot(2,1,1);
        xlabel('nth Energy state');
        ylabel('Energy of state')
        plot(xv, eStates); hold on;
        subplot(2,1,2);
        xlabel('nth Energy state');
        ylabel('exp(-beta*En)')
        plot(xv, pInZ); hold on;    
        zArrman(i) = thisZ;
    end 
    figure(2)
    xv = linspace(1,npoints, npoints);
    xlabel('L')
    ylabel(['Z'])
    plot(xv, zArrman); hold on;
end


%% 2) fermi dirac statistics
%F(E) = 1/(exp((E-mu)*beta)+1) = nr
%sum(nr) = N
%dF/dN = mu
%mu = -1/beta d(lnZ)/dN

nstates = 100;
beta = 1;

npoints = 20;
El = zeros(1, npoints);
Zl = zeros(1, npoints);
bArr = [1, 1e-1, 1e-2,1e-3,1e-5, 1e-7, 1e-8, 1e-10];
for b=1:length(bArr)
    beta = bArr(b);
    for l=1:npoints
        mu = getEnInf(1,l, hb,m);
        FD = zeros(1, nstates);
        Estates = zeros(1, nstates);
        figure(1);
        xlabel('En');
        ylabel('FD(En)');
        for i = 1:nstates
            Estates(i) = getEnInf(i, l, hb, m);
            FD(i) = 1/(exp((Estates(i)-mu)*beta)+1);
        end
        xv = linspace(1,nstates, nstates);
        disp(['sum of nr, beta:', string(beta), 'length', string(l)]);
        disp(sum(FD));
        FD = FD/sum(FD); %ehem....

        plot(xv, FD); hold on;
        El(l) = sum(FD.*Estates);
        Zl(l) = sum(log(1+exp(-beta.*(El(l))-mu)));
    end
figure(2);
xv = linspace(1,npoints, npoints);
subplot(2,1,1);
plot(xv, El); hold on;
xlabel('L');
ylabel('E');
subplot(2,1,2);
plot(xv, Zl); hold on;
xlabel('L');
ylabel('Z');
end

%E = sum(nr*Er)
%% 3) compare energy levels
Elevels  = getElevelsDP(50, 0); %not degenerate
nstates = 50;
infEnergies = zeros(1,nstates);
krakenEnergies = zeros(1,nstates);
infExp = zeros(1,nstates);
krakenExp = zeros(1, nstates);
nx = 2^10;
nbar = 0;
[krakenState, krakenE] = sparce_magic_nd_exec_all(nbar, nx, L);

for i = 1:nstates
    infEnergies(i) = getEnInf(Elevels(i,1), L, hb, m) + getEnInf(Elevels(i,2), L, hb, m);
    infExp(i) = exp(-beta*infEnergies(i));
    krakenEnergies(i) = krakenE(Elevels(i,1)) + krakenE(Elevels(i, 2));
    krakenExp(i) = exp(-beta*krakenEnergies(i));
end

%% 4) En(n) and exp(-beta*En)(n)
% generate occupation factor for different temperatures

betaArr = [1, 0.5, 0.1, 0.01];
ntemp = length(betaArr);
infExp = zeros(ntemp,nstates);
krakenExp = zeros(ntemp, nstates);

for i=1:ntemp
    for j =1:nstates
        infExp(i,j) = exp(-betaArr(i)*infEnergies(j));
        krakenExp(i,j) = exp(-betaArr(i)*krakenEnergies(j));
    end
end
%% plot
xv = linspace(1,50,50);
subplot(1,2,1)
plot(xv, infEnergies, '-or'); hold on;
plot(xv, krakenEnergies, '-og'); hold on;
xlabel('n'); 
ylabel('En');
subplot(1,2,2)
for i = 1:ntemp
plot(xv, infExp(i,:), '-r'); hold on;
plot(xv, krakenExp(i,:), '-g'); hold on;

end
xlabel('n'); 
ylabel(['exp(-beta*En)', string(betaArr)]);
%imagesc(probArray(p1,p2))


%% 5) grand canonical ensemble basics

%% 6) manual canonical ensemble for 1p in infinite well

clear all;
nstates = 10000;
L= 10;
hb = 1;
m = 1;
beta = 1;

%bArr = [100, 10, 1]
%bArr = [1e-1, 1e-2,1e-3];
bArr = [1e-5, 1e-7, 1e-8, 1e-10];
npoints = 20;
Z = zeros(1,npoints);
dZ = zeros(1,npoints);
E = zeros(1,npoints);
S = zeros(1,npoints);
for b = 1:length(bArr)
    beta = bArr(b);
    for l = 1:npoints
        lambda = 1/l*hb*(2*pi*beta/m)^(3/2);
        if lambda>exp(5/2)
            disp(beta)
            disp(lambda)
        end

        for i=1:nstates
            Z(l) = Z(l) + exp(-beta*hb^2*pi^2*i^2/(2*m*l^2));
            dZ(l) = dZ(l) +  beta*hb^2*pi^2*i^2/(2*m*l^3)*exp(-beta*hb^2*pi^2*i^2/(2*m*l^2));
        end

        for i = 1:nstates
            thisE =  hb^2*pi^2*i^2/(2*m*l^2);
            E(l) = E(l) + thisE*exp(-beta*thisE);
        end

        E(l) = E(l)/Z(l);
        S(l) = log(Z(l)) + beta*E(l);

    end
    figure(3);
    xv = linspace(1, npoints, npoints);
    subplot(2,2,1);plot( xv, Z);ylabel('Z(l)'); hold on;
    subplot(2,2,2);plot(xv, dZ);ylabel('dZ/dl(l)'); hold on;
    subplot(2,2,3);plot(xv, E); ylabel('E(l)'); hold on;
    subplot(2,2,4);plot(xv, S); ylabel('S(l)'); hold on;

end

%% reproduce from paper: N=3, g=0, W(T), double well


%% a) general formula, only calculate pn

%W*beta = sum(pn*ln(pn/pnrem))
%pn = sum(exp(- Ei*beta))/Z
%Z = sum(N sum(exp(-Ei*beta)))

%assumption 1: Eni = n*Ei dame! -> Eni = Ei for g= 0
clear all
nstates = 1e5;
L= 10;
npart = 2;
hb =1;
m=1;
beta = 1;
betaAr = [1, 1e-1, 1e-2, 1e-3, 1e-4];
TAr = [0.01, 0.1, 1, 10, 100, 1000];
%TAr = [10];
%TAr=[1e-6]
%% quick calc for Z, E T=10 in cycle

TAr = [10];
for b=1:length(TAr)
    beta = 1/TAr(b);
    Z = 0;
    E = 0;
  for n =1:length(Estates)
       thisE = getEnInf(Estates(n,1), L,hb, m) + getEnInf(Estates(n,2), L, hb, m);
        Z = Z + exp(-beta*thisE);
        E = E + thisE*exp(-beta*thisE);
  end
  E = E/Z
end

%% compariso, passed!
[Zp, Ep] = getPartFkt(Estates,  length(Estates), 2, L, hb,m,beta);

%% W
TAr = [10];
WAr = zeros(1, length(betaAr));
EAr = zeros(1, length(betaAr));
for b =1:length(TAr)
    beta = 1/TAr(b);
    Z = zeros(1, npart+1);
    E= zeros(1, npart+1);
    pnAr = zeros(1, npart+1);
    %{
    for i=1:npart+1
        for n = 1:nstates
        En = getEnInf(n, L/2, hb, m);
        Z(i) = Z(i)+ exp(-En*beta);
        end
    Ztotal = sum(Z);
    end
    %}
    l  = L/2;
    for i = 1:npart+1 %i particles on the left side
            %n p on left side = i-1
            for n =1:length(Estates)
                thisE = 1/2*((i-1)*getEnInf(Estates(n,1), L-l,hb, m) + ... 
                (npart+1-i)*getEnInf(Estates(n,2), l, hb, m) + (i-1)*getEnInf(Estates(n,2), L-l, hb, m) + ...
                (npart+1-i)*getEnInf(Estates(n,1), l,hb, m));
                Z(i) = Z(i) + exp(-beta*thisE);
                E(i) = E(i) + thisE*exp(-beta*thisE);
            end
    end
    E = E./Z
    Ztotal = sum(Z);
    pnAr = Z/Ztotal;

    WAr(b)= -1/beta* sum(pnAr.*log(pnAr./pnMaxAr(:,b).'));
end
%% plot 
plot(TAr, WAr./TAr/log(2))
set(gca, 'XScale', 'log')
xlabel('T')
ylabel('W/T/ln(2)')


 %% get pnrem-> vary Lrem
 npoints = 9;
 pnMaxAr = zeros(npart+1, length(TAr));
 %Estates = getElevelsTP(nstates, 1);
 for b =1:length(TAr)
    beta = 1/TAr(b);
    pnremAr = zeros(npart+1, npoints);
     for l=1:npoints
        Z= zeros(1, npart+1);
        for i = 1:npart+1 %i particles on the left side
            %n p on left side = i-1
            for n =1:nstates
                thisE = (i-1)*getEnInf(n, L-l,hb, m) + (npart+1-i)*getEnInf(n, l, hb, m);
                Z(i) = Z(i) + exp(-beta*thisE);
            end
        end
        Ztotal = sum(Z);
        pnremAr(:,l) = Z/Ztotal;
     end
     xv= linspace(1,npoints,npoints);
     for i = 1:npart+1
        plot(xv, pnremAr(i,:)); hold on; 
        xlabel('L');
        ylabel('pnrem(l)');
        %legend('0', '1', '2', '3');
        pnMaxAr(i, b) = max(pnremAr(i,:))
     end
     
 end
% next step: find maximum pnrem for each temperature and create array to
% use for w calculation
 %% 2 particle, sum over energy of microstate, not single state only, get pnrem-> vary Lrem
 npoints = 9;
 pnMaxAr = zeros(npart+1, length(TAr));
 maxlevel = 1000;
 Estates = getElevelsDP(maxlevel, 1);
 Etemp = zeros(length(TAr), npoints);
 
 for b =1:length(TAr)
    beta = 1/TAr(b);
    pnremAr = zeros(npart+1, npoints);
    EAr = zeros(npart+1, npoints);
     for l=1:npoints
        Z= zeros(1, npart+1);
        E= zeros(1, npart+1);
        for i = 1:npart+1 %i particles on the left side
            %n p on left side = i-1
            for n =1:length(Estates)
                thisE = 1/2*((i-1)*getEnInf(Estates(n,1), L-l,hb, m) + ... 
                (npart+1-i)*getEnInf(Estates(n,2), l, hb, m) + (i-1)*getEnInf(Estates(n,2), L-l, hb, m) + ...
                (npart+1-i)*getEnInf(Estates(n,1), l,hb, m));
                Z(i) = Z(i) + exp(-beta*thisE);
                E(i) = E(i) + thisE*exp(-beta*thisE);
            end
        end
        Ztotal = sum(Z);
        pnremAr(:,l) = Z/Ztotal;
        EAr(:,l) = E./Z
        EAr(isnan(EAr)) = 0;
     end
     
     xv= linspace(1,npoints,npoints);
     figure(2)
     for i = 1:npart+1
        plot(xv, pnremAr(i,:)); hold on; 
        xlabel('L');
        ylabel('pnrem(l)');
        %legend('0', '1', '2', '3');
        pnMaxAr(i, b) = max(pnremAr(i,:))
     end
     
     %{
     %%
    
     figure(3)
     for i = 1:npart+1
        plot(xv, EAr(i,:)); hold on; 
        xlabel('L');
        ylabel('E(l)');
        legend('0', '1', '2')
     end
     pause
    %}
 end
% next step: find maximum pnrem for each temperature and create array to
% use for w calculation

