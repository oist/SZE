function[z, ez]= getPartFkt(elevels,  nstates, npart, L, hb,m,beta)

evals = zeros(1, nstates);
for i=1:nstates
    for j=1:npart
        evals(i) = evals(i) + getEnInf(elevels(i,j), L, hb, m);
    end
end
expfunc =  exp(-beta*evals); 
forexpval = evals.*exp(-beta*evals);
z = sum(expfunc);
ez = sum(forexpval)/z; 