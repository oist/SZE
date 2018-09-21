clear all;
clf;

%% Settings
npart = 2; %thou shalt not change the particle number
nb = 1;
L = 10;
hb = 1;
m = 1;
maxlevel = 1000;
Escale = getEnInf(1,L, hb, m); % ground state of system, use to scale
T = 10;
%% Work formula approach
% W = -beta*sum(pn*log(pn/pnrem))
TAr = [0.01, 0.1, 1, 10, 100, 1000, 1e5, 1e7, 1e9, 1e11, 1e15, 1e20];
bArr = 1./TAr;
nbAr = [1,2,3,4,5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 20, 30];


%% energy states
Estates = getElevelsDP(maxlevel, 1);

%% energy values

beta = 1/TAr(8);
l = 9.9;
nstates = 10;
systemE0= zeros(1,nstates);
i = 1;
for n =1:nstates
    thisE = 1/2*((i-1)*getEnInf(Estates(n,1), L-l,hb, m) + ... 
    (npart+1-i)*getEnInf(Estates(n,2), l, hb, m) + (i-1)*getEnInf(Estates(n,2), L-l, hb, m) + ...
    (npart+1-i)*getEnInf(Estates(n,1), l,hb, m));
    systemE0(n) = thisE;
end

systemE2= zeros(1,nstates);
i = 3;
for n =1:nstates
    thisE = 1/2*((i-1)*getEnInf(Estates(n,1), L-l,hb, m) + ... 
    (npart+1-i)*getEnInf(Estates(n,2), l, hb, m) + (i-1)*getEnInf(Estates(n,2), L-l, hb, m) + ...
    (npart+1-i)*getEnInf(Estates(n,1), l,hb, m));
    systemE2(n) = thisE;
end

figure(8)
plot(systemE0, 'b'), hold on;
plot(systemE2, 'r')

%% get pnrem 
% situation is after expansion process, only one barrier left, removal
% position needs to be determined, therefore same as 1barrier case

 npoints = 2; % 2*npoints will be calculated to determine pnrem
 lmin = 0.01*L;
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
            for n =1:length(1000)
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
        %subplot(2,1,1)
        plot(lvalues, pnremAr(i,:)); hold on; 
        xlabel('L');
        ylabel('pnrem(l)');
        %legend('0', '1', '2', '3');
        
        %{
        subplot(2,1,2)
        plot(lvalues, EAr(i,:)); hold on;
        xlabel('L')
        ylabel('E')
        %}
        pnMaxAr(i, b) = max(pnremAr(i,:));
     end
 end
 return

 %%
 npoints = 2; % values to be taken for l
 lmin = 1;
 L=10;
 hb =1;
 m= 1;
 %lmax = L - 0.01*L;
 n=1;
 lhalf = linspace(lmin, L/2, npoints);
 lhalf2 = lhalf + L/2 - lhalf(1);
 lvalues = [lhalf, lhalf2(2:end)];
 El =zeros(3,length(lvalues));
 Er =zeros(3,length(lvalues));
 table = zeros(3, 6);
% lvalues = [2];
for p =1:length(lvalues)
    l = lvalues(p);
     table = zeros(3, 4);
     for i=1:1
     %El(i,m) =  (i-1)*getEnInf(Estates(n,1), L-l,hb, m);
     %Er(i,m) =  (npart+1-i)*getEnInf(Estates(n,1), l, hb, m);
     table(i,1)=L-l;
     table(i,2)=l;
     table(i,3)= hb^2 *pi^2 * n^2/(2*m*l^2);
     table(i,4)= hb^2 *pi^2 * n^2/(2*m*(L-l)^2);
     end
     table  
end


%%
for i=1:3
   plot(lvalues, El(i,:)); hold on; plot(lvalues, Er(i,:)) 
end