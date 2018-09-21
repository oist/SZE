% make script that reads energy states for same and different wells and
% calculates pn and work, takes pnremD
% returns one value for work done for specific input scenario
function[work] = getWorkInt(esame, ediff, premdiff, T, nb)
newlen = min(length(esame), length(ediff));
esame = esame(1:newlen);
ediff = esame(1:newlen);

beta = 1/T;
Zsame= getZ(esame, beta);
Zdiffa = getZ(ediff, beta);
nS = nb+1;
nD = (nb+1)*nb/2;
Ztot = nS*Zsame + nD*Zdiffa;
pnS= Zsame/Ztot;
pnD = Zdiffa/Ztot;
work = -1/beta*(nS*pnS*log(pnS) + nD*pnD*log(pnD/premdiff));
end
