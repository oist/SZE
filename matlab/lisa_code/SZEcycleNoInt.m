function[] = SZEcycleNoInt(nb, T, L, hb, m, Estates)
    %% Step 1
    E1 = 0;
    Z1 = 0;
    beta = 1/T;
    for i=1:length(Estates)
        thisE = getEnInf(Estates(i,1), L, hb,m) + getEnInf(Estates(i,2), L, hb,m);
        E1 = E1 + thisE*exp(-beta*thisE);
        Z1 = Z1 + exp(-beta*thisE);

    end
    E1 = E1/Z1;
    %% Step 2 and Step 3
    nscenario = (nb+1)*nb/2 + (nb+1);
    Z3 = 0;
    E3= 0;
    l  = L/(nb+1); 
    for n =1:length(Estates)
        thisE = getEnInf(Estates(n,1), l,hb, m) + getEnInf(Estates(n,2), l, hb, m);
        Z3 = Z3 + exp(-beta*thisE);
        E3 = E3 + thisE*exp(-beta*thisE);
    end
    E3 = E3/Z3;
    E2 = nscenario*E3;
    Z2 = nscenario*Z3;

    %% Step 4, fractions should be takem from pnrem!
    [Zhalf, Ehalf] = getPartFkt(Estates,  length(Estates), 2, L/2, hb,m,beta);
    psame= 2/3;
    pdiff = 1/3;
    Z4 = psame*Z1 + pdiff*Zhalf;
    E4 = psame*E1 + pdiff*Ehalf;

    %% cycle
    leCycle(Z1, Z2, Z3, Z4, E1, E2, E3, E4, beta);
end