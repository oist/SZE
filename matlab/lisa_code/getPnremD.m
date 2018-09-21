% get prem
% assumption: barrier in the middle is equilibrium position (two particle
% case), only pnrem for L/2 considered
% --> esame and ediff should be energy states for nb = 1 case
function[premdiff] = getPnremD(esame, ediff,T) 
    beta = 1/T;
    [zd, ed]= getZ(ediff, beta);
    [zs, es] = getZ(esame, beta);
    zs2 = zs;
    ztotal = zs + zd + zs2;  %0, 1, 2, particles in left well
    premdiff = zd/ztotal;
end