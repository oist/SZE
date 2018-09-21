%no idea what this was supposed to do

function [probArray]= probArray(psiA,psiB)
disp('starting')
if size(psiA) ~= size(psiB)
    disp('Arrays not the same size')
    return
else
    np = length(psiA);
end

probArray = zeros(np);
for i = 1:np
    for j = 1:np
    probArray(i,j) = abs(psiA(i)*psiB(j) - psiA(j)*psiB(i))^2;
    end
end