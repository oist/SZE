function[QZ] = getAnalyticPartFkt(beta, l, hb, m)
e1 = hb^2*pi^2/(2*m*l^2);
QZ = sqrt(pi)/2*(1/beta/e1)^(1/2) + (-1/90*(e1*beta)^3 + 1/60*(e1*beta)^2 ...
    + 1/6*(e1*beta) +1/2)*exp(-e1*beta);
end