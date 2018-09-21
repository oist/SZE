function[energies] = getEnergiesA(l, maxlevel, estat, nstat, hb, m)
    %get Elevels
    eSingleA = zeros(1,maxlevel);
    for i = 1:maxlevel
        eSingleA(i) = getEnInf(i, l, hb, m);
    end
    eDa = zeros(1, length(estat));
    for i = 1:length(estat)
        eDa(i) = eSingleA(estat(i,1)) + eSingleA(estat(i,2));
    end
    eDa = sort(eDa);
    eDa = eDa(1:nstat);
    energies = eDa;
end