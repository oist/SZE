%get wave function infinite well of length l in state n 
%returns array with psi values at x positions with resolution set by
%npoints
function [psiArr, xv]= psiInfWell(L,k, npoints)

psiArr = zeros( 1, npoints);
xv = zeros( 1, npoints);
xspacing = vpa(L/(npoints-1));
for i = 0:npoints-1
    x = vpa(i)*xspacing;
    psiArr(i+1) = vpa(sin(k*x));
    xv(i+1) = x;   
end
%normalize psi
n = norm(psiArr);
psiArr = psiArr/n;

end



