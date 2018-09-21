function[W] = getWorkNoInt(Estates, beta, nb, l, thisPnRem, hb, m)
    Z = 0;
    E = 0;
    nscenarios = (nb+1)*nb/2 +nb+1;
    pnAr = ones(1, nscenarios);
    for n =1:length(Estates)
        thisE = getEnInf(Estates(n,1), l,hb, m) + getEnInf(Estates(n,2), l, hb, m);
        Z = Z + exp(-beta*thisE);
        E = E + thisE*exp(-beta*thisE);
    end
    E = E/Z;
    Ztotal = nscenarios*Z;
    pnAr = pnAr*Z/Ztotal;
    %pnAr = 1/nscenarios.*pnAr;
    remAr = zeros(1, nscenarios);
    for i = 1:nb+1
        %remAr(i) = thisPnRem(1);
        remAr(i) = 1; %note: reversible process if in same well, as wall removal does not change system!

    end
    remAr(remAr ==0) = thisPnRem(2);
    W= -1/beta* sum(pnAr.*log(pnAr./remAr));
end