%get reduced sp density matrix for nb = 2 (try to keep general)
%cycle treatment

%has bugs

clear all

%check current folder

if not(string(pwd) == '/home/lisa')
    disp('currently in wrong directory, please change to /home/lisa')
    disp('then run script again')
    return
end

%% Settings
directorywf = 'Documents/OIST/sze_out/matlab/WF/';
directoryrspdm = 'Documents/OIST/sze_out/matlab/RSPDM/';
if not(exist(directorywf,'dir')==7)
    disp(strcat('please make directory ', directorywf))
end
if not(exist(directoryrspdm,'dir')==7)
    disp(strcat('please make directory ', directoryrspdm))
end

makenew = 1;
L = 10;
nbar = 2;
zeroCutOff = 1e-9;
%resolution, i.e. number of points
np = 2^11
%number of particles CAN ONLY DO 2 ATM
npart = 2;
%x spacing
dx = L/(np-1);
%kspacing (??)
dkx = dx;
x=-L./2:dx:L./2;

%% step 1: single well

%check if psi exists, if not, save psi to file
%parameters: np, L, nbar, n = 1 or 2
thisFilename = char(strcat(directorywf, myInfWFfilename(1, nbar, np, L)));
if and(exist(thisFilename, 'file') == 2, not(makenew)) 
    [psi1, xv] = textread(thisFilename);
    psi1 = psi1.'; 
else
    [psi1, xv] = psiInfWell(L, pi/L, np);
    fileID = fopen(thisFilename,'w');
    fprintf(fileID,'%f %f\n',[psi1;xv]);
    fclose(fileID);
  
end   


thisFilename = char(strcat(directorywf, myInfWFfilename(2, nbar, np, L)));
if and(exist(thisFilename, 'file') == 2, not(makenew))
    [psi2,xv] = textread(thisFilename);
    psi2 = psi2.'; 

else
    [psi2, xv] = psiInfWell(L, 2*pi/L, np);  
    fileID = fopen(thisFilename,'w');
    fprintf(fileID,'%f %f\n',[psi2; xv]);
    fclose(fileID);
end

psiAr = [psi1;psi2];

%%

thisFilename = char(strcat(directoryrspdm, myRSPDMfilename('infWell_1', nbar, np, L)));
if exist(thisFilename, 'file') == 2
    rspdm = textread(thisFilename);
else
    [rspdm0, b, c] = dmnk(psiAr,dx,dkx,npart,np);
    rspdm=rspdm0./trace(rspdm0); 
    fileID = fopen(thisFilename,'w');
    fprintf(fileID,'%f\n',rspdm);
    fclose(fileID);
end

S1 = getEntropy(rspdm, zeroCutOff)

%% step 2: ngap split well

thisFilename = char(strcat(directorywf, mySplitWFfilename(1, nbar, np, L)));
if exist(thisFilename, 'file') == 2
    psiS1 = textread(thisFilename);

else
    psiS1 = analyticNbar_noTun_exec(nbar, 1, np, psi1, xv); 
    psiS1 = psiS1.';
    fileID = fopen(thisFilename,'w');
    fprintf(fileID,'%f\n',psiS1);
    fclose(fileID);
end

thisFilename = char(strcat(directorywf, mySplitWFfilename(2, nbar, np, L))); 
if exist(thisFilename, 'file') == 2
    psiS2 = textread(thisFilename);

else
    psiS2 = analyticNbar_noTun_exec(nbar, 2, np, psi2, xv);
    psiS2 = psiS2.';
    fileID = fopen(thisFilename,'w');
    fprintf(fileID,'%f\n',psiS2);
    fclose(fileID);
end

psiArS = [psiS1;psiS2];
%% step 2: entropy
thisFilename = char(strcat(directoryrspdm, myRSPDMfilename('splitWell_2', nbar, np, L)));
if exist(thisFilename, 'file') == 2
    rspdmS = textread(thisFilename);
else
    [rspdmS0, b, c] = dmnk(psiArS,dx,dkx,npart,np);
    rspdmS=rspdmS0./trace(rspdmS0); 
    fileID = fopen(thisFilename,'w');
    fprintf(fileID,'%f\n',rspdmS);
    fclose(fileID);
end
S2 = getEntropy(rspdmS, zeroCutOff)

%% step 3: 2 gap split well

thisFilename = char(strcat(directorywf, mySplitWFfilename(1, 1, np, L)));
if exist(thisFilename, 'file') == 2
    psiE1 = textread(thisFilename);

else
    psiE1 = analyticNbar_noTun_exec(1, 1, np, psi1, xv); 
    psiE1 = psiE1.';
    fileID = fopen(thisFilename,'w');
    fprintf(fileID,'%f\n',psiE1);
    fclose(fileID);
end

thisFilename = char(strcat(directorywf, mySplitWFfilename(2, 1, np, L))); 
if exist(thisFilename, 'file') == 2
    psiE2 = textread(thisFilename);
else
    psiE2 = analyticNbar_noTun_exec(1, 2, np, psi2, xv);
    psiE2 = psiE2.';
    fileID = fopen(thisFilename,'w');
    fprintf(fileID,'%f\n',psiE2);
    fclose(fileID);
end

psiArE = [psiE1;psiE2];

%% step 3 entropy

thisFilename = char(strcat(directoryrspdm, myRSPDMfilename('doubleWell_3', nbar, np, L)));
if exist(thisFilename, 'file') == 2
    rspdmE = textread(thisFilename);
else
    [rspdmE0, b, c] = dmnk(psiArE,dx,dkx,npart,np);
    rspdmE=rspdmE0./trace(rspdmE0); 
    fileID = fopen(thisFilename,'w');
    fprintf(fileID,'%f\n',rspdmE);
    fclose(fileID);
end

S3 = getEntropy(rspdmE, zeroCutOff)

%calculating entropy:
%diagonalize rspdm, sum over plnp where p are diagonal elements

%% save entropy
thisFilename = char(strcat(directoryrspdm, myCycleEntropyfilename(nbar, np, L)));
fileID = fopen(thisFilename,'w');
fprintf(fileID,'%s\n','S1');
fprintf(fileID,'%f\n',S1);
fprintf(fileID,'%s\n','S2');
fprintf(fileID,'%f\n',S2);
fprintf(fileID,'%s\n','S3');
fprintf(fileID,'%f\n',S3);
fclose(fileID);

%% plotting wf
clf
figure(1)
subplot(3,2,1); plot(xv, psi1); hold on;
subplot(3,2,2); plot(xv, psi2); hold on;
subplot(3,2,3); plot(xv, psiS1); hold on;
subplot(3,2,4); plot(xv, psiS2); hold on;
subplot(3,2,5); plot(xv, psiE1); hold on;
subplot(3,2,6); plot(xv, psiE2); hold on;

%% save plot
filename = char(strcat(directorywf, 'wfCycle_noTun_nbar',string(nbar),'_np_', string(np), '.png'))
saveas(figure(1), filename);
close(figure(1))


