function[limWI] = getTInfLimW(nb)
limWI = 1/(nb+1)*(2*log(nb+1)*(1+nb) - 2*nb*log(2));
end
