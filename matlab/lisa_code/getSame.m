function[zsame, esame]= getSame(beta, nx, l, nstat, g)
%get Elevels
[psiSame, eS] = sparce_magic_nd_exec_all(0, nx, l, 2, nstat, g);
%figure(3)
%plot(eS); hold on;
%get Z
[zsame,esame] = getZ(eS, beta);


end