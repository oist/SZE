%reduced density matrix for szilard engine
clear all;
L= 10
np = 2^6
format long
%get psi arrays (single particle solutions)
psiA = psiInfWell(L, pi/L, np);
psiB = psiInfWell(L, 2*pi/L, np);
xv = linspace(0,L, np);
%plot single array
%p = plot(xv,psiA);

%plot many body WF
%case: unsplit well
probArray = zeros(np);
for i = 1:np
    for j = 1:np
    probArray(i,j) = abs(psiA(i)*psiB(j) - psiA(j)*psiB(i))^2;
    end
end
%imagesc(probArray)
%return
%case: split well
k1 = 
k2 = 
psiA2 = psiInfWell(L, k1, np);
psiB2 = psiInfWell(L, k2, np);
%psiA2 = [psiA2, psiA2]
%psiB2 = [psiB2, psiB2]
xv2 = linspace(0,L, np);

size(psiA2)
size(xv2)
%plot(xv2, psiA2)
%return
disp('le wavefunctions')
psi1 = psiA2
%psi1(end) = 0
psi2 = psiB2
if size(psi1) ~= size(psi2)
    disp('Arrays not the same size')
    return
else 
    len = size(psi1);
    np = len(2);
end

probArray = zeros(np);
probArray2 = zeros(np);
probArray3 = zeros(np);

for i = 1:np
    for j = 1:np
        probArray(j,i) = vpa(abs(psi1(i)*psi2(j) - psi1(j)*psi2(i))^2);% psi1(j)*psi2(i))^2;
        probArray2(i,j) = psi1(j)*psi2(i);    
        probArray3(i,j) = psi1(i)*psi2(j);    
        %probArray(i,j) = abs(psi1(j)*psi2(i));
    end
end

probArray
probArray2
probArray3

%pArr2 = probArray(psiA2, psiB2)
imagesc(probArray3)
%imagesc(probArray)
%compare solutions to irinas solutions
%calculate entropy as done naively
%calculate entropy from rspdm