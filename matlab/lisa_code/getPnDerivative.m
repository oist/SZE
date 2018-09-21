function[dPnEn] = getPnDerivative(n, L, hb, m, beta, cutOff)
Evals = zeros(1, cutOff);
for i=1:cutOff
   Evals(i) = getEnInf(i, L, hb, m);   
end
expFunction =  exp(-beta*Evals); 
Z = sum(expFunction);
if Z == 0
    dPnEn = 0;
else
    dPnEn = beta/Z*exp(-beta*getEnInf(n, L, hb, m))*(1/Z-1);
end

end