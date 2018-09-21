function [S]= getEntropy(rspdm, zeroCutOff)
    ev1 = eig(rspdm);
    %ev1(abs(ev1)<zeroCutOff) =0;
    ev1log = zeros(size(ev1));
    for i= 1:length(ev1)
        if ev1(i) ==0
            ev1log(i) = 0;
        else
            ev1log(i) = ev1(i)*log(ev1(i));
        end
    end
    S = -sum(ev1log);
end