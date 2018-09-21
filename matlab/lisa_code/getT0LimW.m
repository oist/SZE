function[limW0] = getT0LimW(nb)
limW0 = 1/(nb+2)*(log((nb+1)*(nb+2))*(2+nb) - 2*log(2) - nb*log(6));
end
