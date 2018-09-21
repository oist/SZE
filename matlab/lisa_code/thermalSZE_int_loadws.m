% read workspace or read .out file and obtain elevels

workspaceDir= 'w';
gvals= [-1, 0, 1];
nbvals = [1,2,3];
TAr = [0.01,0.02, 0.05, 0.1, 0.2, 0.5, 1, 2];
%% FROM OTHER CODE

%% pnrem 
Zsame = zeros(1, length(TAr), length(gvals));
Zdiff = zeros(1, length(TAr),length(gvals));
for j=1:length(gvals)
    %load workspace
    thisWS = load(char(strcat('thermalSZE_InteractingBosons_slurm_workspace_gAr_', string(gvals(j)),'_nbAr_1.mat')));
    nstat = length(thisWS.eS);
    for b=1:length(TAr)
        beta = 1/TAr(b);
        Zsame(1,b,j) = getZ(thisWS.eS, beta);
        [Zdiff(1,b,j), Ediff] = getDiffA(beta, L/2, maxlevel, EstatesD, nstat, hb, m);
    end
end
Ztot = 2*Zsame + Zdiff;
premI = Zdiff./Ztot;

%% pn

eD=zeros(length(nbvals), nstates);
Zs=zeros(length(nbvals), length(TAr), length(gvals));
Zd=zeros(length(nbvals), length(TAr));

for j=1:length(gvals)
    for k = 1:length(nbvals)
        nb = nbvals(k);
        thisWS = load(char(strcat('thermalSZE_InteractingBosons_slurm_workspace_gAr_', string(gvals(j)),'_nbAr_',string(nb),'.mat')));
        nstat = length(thisWS.eS);
        l = L/(nb+1);
        %eD(k,:) = getEnergiesA(l, maxlevel, EstatesD, nstates, hb, m);
        for b=1:length(TAr)
            beta = 1/TAr(b);
            Zs(k,b,j) = getZ(thisWS.eS, beta);
            %Zd(k,b) = getZ(eD(k,:), beta);
        end
    end
end
for k = 1:length(nbvals)
    nb = nbvals(k);
    l = L/(nb+1);
    eD(k,:) = getEnergiesA(l, maxlevel, EstatesD, nstates, hb, m);
    for b=1:length(TAr)
        beta = 1/TAr(b);
        Zd(k,b) = getZ(eD(k,:), beta);
    end
end


%%
ns = (nbvals+1);
nd = (nbvals+1).*nbvals/2;
Ztotal = ns.'.*Zs + nd.'.*Zd;
ps = Zs./Ztotal;
pd = Zd./Ztotal;

betaWI = - (ns.'.*ps.*log(ps) + nd.'.*pd.*log(pd./premI));

%% plotting area

%%
figure(8)
for k=1:length(nbvals)
    plot(TAr, betaWI(k,:,1)/log(nbvals(k)+1)); hold on;
end
title('W(T), varying nb, g = -0.1');
set(gca, 'XScale', 'log');

figure(9)
for k=1:length(nbvals)
    plot(TAr, betaWI(k,:,2)/log(nbvals(k)+1)); hold on;
end
title('W(T), varying nb, g = 0');
set(gca, 'XScale', 'log');

figure(10)
for k=1:length(nbvals)
    plot(TAr, betaWI(k,:,3)/log(nbvals(k)+1)); hold on;
end
title('W(T), varying nb, g = 0.1');
set(gca, 'XScale', 'log');

%%

figure(11)
for j=1:length(gvals)
    plot(TAr, betaWI(1,:,j)/log(1+1)); hold on;
    gvals(j)
end
legend('g=-1', 'g=0', 'g=1');
xlabel('T [10^x]');
ylabel('W/T/log(nb+1)');
title('W(T), varying g, nb = 1');
set(gca, 'XScale', 'log');

figure(12)
for j=1:length(gvals)
    plot(TAr, betaWI(2,:,j)/log(2+1)); hold on;
end
xlabel('T [10^x]');
ylabel('W/T/log(nb+1)');
legend('g=-1', 'g=0', 'g=1');
title('W(T), varying g, nb = 2');
set(gca, 'XScale', 'log');

figure(13)
for j=1:length(gvals)
    plot(TAr, betaWI(3,:,j)/log(3+1)); hold on;
end
xlabel('T [10^x]');
ylabel('W/T/log(nb+1)');
legend('g=-1', 'g=0', 'g=1');
title('W(T), varying g, nb = 3');
set(gca, 'XScale', 'log');

