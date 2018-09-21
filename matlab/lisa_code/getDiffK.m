function[zdiff, ediff] = getDiffK(beta, nx, l, maxlevel, estat, nstat)
    [psiDiff, eSingle] = sparce_magic_nd_exec_all(0, nx, l, 1,maxlevel,0);
    eD = zeros(1, length(estat));
    for i = 1:length(estat)
        eD(i) = eSingle(estat(i,1)) + eSingle(estat(i,2));
    end
    eD = sort(eD);
    eD = eD(1:nstat);
    %get Z
    [zdiff,ediff] = getZ(eD, beta);
end