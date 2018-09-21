% two particles, TG limit interaction in nbarrier infinite well potential
% measurement determines position of both particles
% approach: partition function for canonical ensemble
clear all;

%% output settings
plotWT = 0;
constantT = 0.1;
plotWnb = 0;
constantnb = 3;

plotWTnb = 1;

%% code settings
npart = 2; %thou shalt not change the particle number
L = 10;
hb = 1;
m = 1;
maxlevel = 5000;
Escale = getEnInf(1,L, hb, m); % ground state of system, use to scale

TAr = [0.01, 0.1, 1, 10, 100, 1000,1e4, 1e5, 1e6,1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13,1e14,  1e15, 1e16, 1e17, 1e18, 1e19, 1e20];%1e11, 1e15, 1e20];

%TAr = [1e16, 1e17, 1e18, 1e19, 1e20];%1e11, 1e15, 1e20];
nbAr = [1,2,3,4,5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,16,17,18,19,20, 21, 22, 23,24,25, 26, 27, 28, 29, 30];

%TAr = [0.01, 0.1];% 1, 10, 100, 1000,1e4, 1e5, 1e6,1e7, 1e8, 1e9, 1e10];%1e11, 1e15, 1e20];
%nbAr = [1,2,3];%,4,5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,16,17,18,19,20];
%% energy states
Estates = getElevelsDP(maxlevel, 1); % n occupation numbersfor particles in different wells
EstatesSame = getElevelsDP(maxlevel, 0); % n occupation numbers for particles in same well

%% get pnrem 
% situation is after expansion process, only one barrier left, removal
% position needs to be determined, therefore same as 1barrier case

 npoints = 10; % values to be taken for l
 lmin = 0.001*L;
 lhalf = linspace(lmin, L/2, npoints);
 lhalf2 = lhalf + L/2 - lhalf(1);
 lvalues = [lhalf, lhalf2(2:end)];
 pnMaxAr = zeros(npart+1, length(TAr));
 
for b =1:length(temperature)
    beta = 1/temperature(b);
    pnremArS = zeros(1, length(lvalues));
    pnremArD = zeros(1, length(lvalues));
    psame = zeros(1, length(lvalues));
    pdiff = zeros(1, length(lvalues));
    Zs =  zeros(1, length(lvalues));
    Es =  zeros(1, length(lvalues));
    Ztotal =  zeros(1, length(lvalues));
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
        

        
        %[Zs(j,p), Es(j,p)] = getSame(beta, nx, l, nstatAr(b), gAr(j));
        [Zs(p), Es(p)] = getDiffA(beta, l, maxlevel, EstatesSame,nstates, hb, m);
            
    end

    Zs2 = fliplr(Zs);
    
    
    Ztotal = Zs + Zda + Zs2;  
    pnremArS = Zs./Ztotal;
    pnremArD = Zda./Ztotal;
    premMaxSame(j, b) = max(pnremArS(j,:));
    premMaxDiff(j, b) = max(pnremArD(j,:));
    
end


%% old
 npoints = 10; % values to be taken for l
 lmin = 0.001*L;
 lhalf = linspace(lmin, L/2, npoints);
 lhalf2 = lhalf + L/2 - lhalf(1);
 lvalues = [lhalf, lhalf2(2:end)];
 pnMaxAr = zeros(npart+1, length(TAr));
 
 figure(7)
 for b =1:length(TAr)
    beta = 1/TAr(b);
    pnremAr = zeros(npart+1, length(lvalues));
    EAr = zeros(npart+1, length(lvalues));
    for p=1:length(lvalues)
        l = lvalues(p);
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
        pnremAr(:,p) = Z/Ztotal;
        EAr(:,p) = E./Z;
        EAr(isnan(EAr)) = 0;
     end
     
     for i = 1:npart+1
        subplot(2,1,1)
        plot(lvalues, pnremAr(i,:)); hold on; 
        xlabel('L');
        ylabel('pnrem(l)');
        %legend('0', '1', '2', '3');
        subplot(2,1,2)
        plot(lvalues, EAr(i,:)); hold on;
        xlabel('L')
        ylabel('E')
        
        pnMaxAr(i, b) = max(pnremAr(i,:));
     end
 end
save('../../thermalWorkspace/pnrem')
%% plot W(T)
if plotWT
    nbK = [1,2,3,4];
    for k=1:length(nbK)
    %nb = constantnb;
    nb = nbK(k);
    nscenarios = (nb+1)*nb/2 +nb+1;
    WAr = zeros(1, length(TAr));
    EAr = zeros(1, length(TAr));
    for b =1:length(TAr)
        beta = 1/TAr(b);
        pnAr = ones(1, nscenarios);
        l  = L/(nb+1); 
        Z = 0;
        E = 0;
        for n =1:length(Estates)
            thisE = getEnInf(Estates(n,1), l,hb, m) + getEnInf(Estates(n,2), l, hb, m);
            Z = Z + exp(-beta*thisE);
            E = E + thisE*exp(-beta*thisE);
        end
        E = E/Z;
        Ztotal = nscenarios*Z;
        pnAr = pnAr*Z/Ztotal;
        remAr = zeros(1, nscenarios);
        for i = 1:nb+1
            remAr(i) = pnMaxAr(1,b);
        end
        remAr(remAr ==0) = pnMaxAr(2,b);
        WAr(b)= -1/beta* sum(pnAr.*log(pnAr./remAr));
    end
    
    %% plot W for different temperatures
    figure(3)
    plot(TAr, WAr./TAr/log(nb+1)/Escale); hold on;
    plot(TAr, log(nscenarios/3)*ones(1,length(TAr))/log(nb+1)/Escale);
    set(gca, 'XScale', 'log')
    xlabel('T')
    ylabel('W/T/ln(nb+1)/Escale') 
    legend('W(T)', 'W(T-> Inf)');
    title(['W(T) nbar 1 2 3 ' + string(nb)])
    save('../../thermalWorkspace/WT')
    end
end

%% 1 T varying nbar
if plotWnb
    TK = TAr(1:4);
    for b=1:length(TK)
    beta = 1/TK(b);
    WAr = zeros(1, length(nbAr));
    Wclass = zeros(1, length(nbAr));
    EAr = zeros(1, length(nbAr));
    for k =1:length(nbAr)
        nb = nbAr(k);
        nscenarios = (nb+1)*nb/2 +nb+1;
        pnAr = ones(1, nscenarios);
        l  = L/(nb+1); 
        Z = 0;
        E = 0;
        for n =1:length(Estates)
            thisE = getEnInf(Estates(n,1), l,hb, m) + getEnInf(Estates(n,2), l, hb, m);
            Z = Z + exp(-beta*thisE);
            E = E + thisE*exp(-beta*thisE);
        end
        E = E/Z;
        Ztotal = nscenarios*Z;
        pnAr = pnAr*Z/Ztotal;
        remAr = zeros(1, nscenarios);
        for i = 1:nb+1
            remAr(i) = pnMaxAr(1,b);
        end
        remAr(remAr ==0) = pnMaxAr(2,b);
        WAr(k)= -1/beta* sum(pnAr.*log(pnAr./remAr));
        Wclass(k) = 1/beta*log(nb+1);
    
    end

    %% plot W for varying nb
    
    figure(4)
    plot(nbAr, WAr*beta/Escale); hold on;
    plot(nbAr, Wclass*beta/Escale, 'r'); hold on;
    %figure(2)
    %plot(nbAr, WAr./Wclass,'g'); hold on;
    nscenAr = (nbAr+1).*nbAr./2 +nbAr+1;
    plot(nbAr, -1/beta*beta/Escale*(1./(nbAr./2 +1).*log(1./nscenAr./pnMaxAr(1,b)) + nbAr./(nbAr+2).*log(1./nscenAr./pnMaxAr(2,b))));
    %set(gca, 'XScale', 'log')
    xlabel('nb');
    ylabel('W/T/Escale');
    legend('W(nb)', 'Wclass', 'W(nb) formula');
    title(['W(nb)/T/Escale, T 0.1 1 10' + string(constantT)]);
    save('../../thermalWorkspace/Wnb')
    end
end
%%  W(T -y, nbar - x)
if plotWTnb
    WAr = zeros(length(TAr), length(nbAr));
    for k=1:length(nbAr)
        nb = nbAr(k);
        l  = L/(nb+1); 
        for b=1:length(TAr)
            beta = 1/(TAr(b));
            WAr(b, k) = getWorkNoInt(Estates, beta, nb, l, pnMaxAr(:,b), hb, m)/TAr(b)/log(nb+1)/Escale;
        end
    end
    %%
    figure(8)
    colormap jet
    surf(WAr); hold on;
    imagesc(WAr);
    xlabel('nbar')
    ylabel('T [10^x]')
    title('W(nbar, T)/T/Wclass/Escale')
    save('../../thermalWorkspace/WTnb')
end
return
%%
figure(8)
colormap jet
imagesc([nbAr(1), nbAr(end)],[TAr(1), TAr(end)],WAr)
colormap
xlabel('nbar')
ylabel('T')
title('W(nbar, T)') 
%set(gca, 'YScale', 'log')



%axis scaling...??
%nstates for high T
%normalize by e1 ??