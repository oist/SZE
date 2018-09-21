%thermal distribution for 2 gap well
%look at expectation value of the energy for different parts in the cycle
%take work output as difference in energy expectation values
%compute for different temperatures
%compute for different energies

clear all;
clf;

%% Settings
nbar = 1
L = 10;
hb = 1;
m = 1;
T = 1; % in Kelvin
%kB = 1.38*10^(-23);
kB = 1; 
beta = 1/(kB*T);

cutOff = 2^5;

%% Energy solutions
%E = getEnInf
%Z = sum(exp(-beta*E))
%determine cutoff
%expectation value is given by <E> = -d(lnZ)/dbeta = 1/Z*(sum(E*exp(-beta*E)))

%% Step 1: Infinite Well

Evals1 = zeros(1, cutOff);
for i=1:cutOff
   Evals1(i) = getEnInf(i, L, hb, m);
end
expFunction1 =  exp(-beta*Evals1); 
forExpVal1 = Evals1.*exp(-beta*Evals1);
Z1 = sum(expFunction1);
EZ1 = sum(forExpVal1);
plot(expFunction1); hold on;
plot(forExpVal1);

%% Step 2: Split Well
%expFunction = zeros(1, cutOff);
%forExpVal = zeros(1, cutOff);
Evals2 = zeros(1, cutOff);
for i=1:cutOff
   Evals2(i) = getEnInf(i, L/2, hb, m);
   
end
expFunction2 =  exp(-beta*Evals2); 
forExpVal2 = Evals2.*exp(-beta*Evals2);

Z3 = sum(expFunction2);
Z2 = 2*Z3;
EZ2 = 2*sum(forExpVal2);
EZ3 = sum(forExpVal2);
plot(expFunction2); hold on;
plot(forExpVal2);

%% work
%% expectation value
expectE1 = EZ1/Z1;
expectE2 = EZ2/Z2;
expectE3 = EZ3/Z3; %(no change in energy as system is measured)

Dexpect21 = expectE2 - expectE1; %insertion
Dexpect32 = expectE3 - expectE2; %measurement
Dexpect13 = expectE1 - expectE3; %expansion

%% partition function approach

wZ21 = 1/beta * (log(Z2) - log(Z1));
wZ32 = 1/beta * (log(Z3) - log(Z2));
wZ13 = 1/beta * (log(Z1) - log(Z3));

%% some plot

figure(3)
xv = [1,2,3,4];
zeroLine= zeros(1,4);
plot(xv, zeroLine, '--k'); hold on;
plot(xv,[expectE1, expectE2, expectE3, expectE1], '-.or'); hold on;
plot([1.5, 2.5, 3.5],[wZ21, wZ32, wZ13],'-.ob'); hold on;
plot([1.5, 2.5, 3.5],[Dexpect21, Dexpect32, Dexpect13],'-.og');
plot([1.5, 2.5, 3.5],[Dexpect21+wZ21, Dexpect32+wZ32, Dexpect13+wZ13],'-.om');

disp('le Cycle')
disp('total expectation value energy change')
disp(string(Dexpect21 + Dexpect32 +Dexpect13))
disp('total work')
disp(string(wZ21+ wZ32 + wZ13))
disp('total heat')
disp(string(Dexpect21+wZ21 + Dexpect32+wZ32 + Dexpect13+wZ13))
disp('--------')
disp('ignore measurement')
disp('total expectation value energy change')
disp(string(Dexpect21 + Dexpect13))
disp('total work')
disp(string(wZ21 + wZ13))
disp('total heat')
disp(string(Dexpect21+wZ21  + Dexpect13+wZ13))




%% heat calculation to confirm stuff above

%heat given by dQ = sum(En*dpn)
% step of interest: 21
%{
nx = 2^3;
DL = L - L/2;
binWidth = DL/nx;
PnInt = zeros(1, cutOff);
PInt = 0;
for i =1:cutOff
    dPnAr = zeros(1, nx);
    for j=1:nx
        thisDer = getPnDerivative(i, L - j*binWidth, hb, m, beta, cutOff)
        PInt = PInt + thisDer;
        dPnAr(j) = thisDer;
    end
    PnInt(i) = sum(dPnAr);
    PInt;
end
%}

%% 2 particles, 1 bar

%generate energy levels
nb = 1;
maxLevel = 2^8;
nStates = maxLevel^2- maxLevel*(maxLevel+1)/2;

%non-degenerate
dpLevels = zeros(nStates, 2);
count= 1;
for i=1:maxLevel
    for j = 1:(i-1)
        dpLevels(count,:) = [j,i];
        count = count+1;
        
    end
end

%degenerate dpLevels
dpLevelsDegen = zeros(nStates, 2);
count= 1;
for i=1:maxLevel
    for j = 1:i
        dpLevelsDegen(count,:) = [j,i];
        count = count+1;
    end
end
   
%% DP: Step 1: Infinite Well

dpEvals1 = zeros(1, nStates);
for i=1:nStates
   dpEvals1(i) = getEnInf(dpLevels(i,1), L, hb, m) + getEnInf(dpLevels(i,2), L, hb, m);
end
dpexpFunction1 =  exp(-beta*dpEvals1); 
dpforExpVal1 = dpEvals1.*exp(-beta*dpEvals1);
dpZ1 = sum(dpexpFunction1);
dpEZ1 = sum(dpforExpVal1);
plot(dpexpFunction1); hold on;
plot(dpforExpVal1);
%% test - passed!
[tz, tez]= getPartFkt(dpLevels,nStates, 2, L, hb,m, beta);


%% DP: Step 2: Split Well
%expFunction = zeros(1, cutOff);
%forExpVal = zeros(1, cutOff);
dpEvals2 = zeros(1, nStates);
for i=1:cutOff
   dpEvals2(i) = getEnInf(dpLevels(i,1), L/(nb+1), hb, m) + getEnInf(dpLevels(i,2), L/(nb+1), hb, m);
   
end
disp('1')
dpexpFunction2 =  exp(-beta*dpEvals2); 
disp('1')
dpforExpVal2 = dpEvals2.*exp(-beta*dpEvals2);
disp('1')
dpZ2 = sum(dpexpFunction2);
dpEZ2 = sum(dpforExpVal2);
%i know this from the density plot! pos Arrays are needed to establish the
%number of possibilities --> only 1 in this case!
dpZ3 = dpZ2;
dpZ4 = dpZ3; % i know this from the fact that barrier position is at the middle
dpEZ3 = dpEZ2;
dpEZ4 = dpEZ2;

%% DP: summing up
leCycle(dpZ1, dpZ2, dpZ3, dpZ4, dpEZ1, dpEZ2, dpEZ3, dpEZ4, beta);
%add legend at some point!

%% 2 particles, 2 barriers
nbar = 2;

%% 2P2B: step one:converging Z in 2p case for infinite well -> nstate cutoff
npoints = 50;
count = 1;
zArr = zeros(1,npoints);
ezArr = zeros(1,npoints);
for i =1:round(maxLevel/npoints):maxLevel
    [zArr(count), ezArr(count)] = getPartFkt(dpLevels,  i, 2, L, hb,m,beta);
    count= count+1;
end
xzvals = linspace(1,maxLevel,length(zArr));
plot(xzvals, zArr); hold on; plot(xzvals, ezArr);
%--> 2^8 more than enough, maybe around 50?
%density plot approaches even distribution over gaps (same probability)
%as high energy levels are spread over all gaps: pn fromthermal
%distribution
%next step: how to get partition function?
%which energies in which gaps?

%% scaling of system with L

npoints = 20;
count = 1;
zArr = zeros(1,npoints);
ezArr = zeros(1,npoints);
L = 20;
beta = 0.1;
for i =1:round(L/npoints):L
    [zArr(count), ezArr(count)] = getPartFkt(dpLevels,  50, 2, i, hb,m,beta);
    count= count+1;
end
xzvals = linspace(1,L,length(zArr));
plot(xzvals, zArr, 'b'); hold on; plot(xzvals, ezArr,'r');

%% 2P2B: partition functions
%% 2P2B: step one: infinite well
% different energy levels

[nbZ1, nbEZ1] = getPartFkt(dpLevels,  100, 2, L, hb,m,beta);

%% 2P2B: step two: superposition

[nbZbar, nbEZbar] = getPartFkt(dpLevels, 100, 2, L/(nbar+1), hb, m, beta);
nbZ2 = 6*nbZbar;
nbEZ2 = 6*nbEZbar;

%% 2P2B: step three: measure; and step four: expand
%% A) measurement destroys all coherence, degenerate energy levels in in
%different boxes

[nbZbarDegen, nbEZbarDegen] = getPartFkt(dpLevelsDegen, 100, 2, L/(nbar+1), hb, m, beta);
nbZ3a = nbZbarDegen;
nbEZ3a = nbEZbarDegen;
 
[nbZ4a, nbEZ4a]= getPartFkt(dpLevelsDegen, 100, 2, L/2, hb, m, beta);

%% B) measurement leaves occupation of different energies in tact
%no degenerate energy levels
[nbZbarDegen, nbEZbarDegen] = getPartFkt(dpLevels, 100, 2, L/(nbar+1), hb, m, beta);
nbZ3b = nbZbarDegen;
nbEZ3b = nbEZbarDegen;
 
[nbZ4b, nbEZ4b]= getPartFkt(dpLevels, 100, 2, L/2, hb, m, beta);

%% cycle

leCycle(nbZ1,nbZ2,nbZ3a,nbZ4a,nbEZ1, nbEZ2, nbEZ3b, nbEZ4b, beta);



        

