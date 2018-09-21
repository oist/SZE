function[posArray] = getPosArray(psi1, psi2, nbar)

if not(length(psi1)== length(psi2))
    disp('WF have different lenghts!')
end

np = length(psi1);
probArray = zeros(np);
for i = 1:np
    for j = 1:np
        pn = abs(psi1(i)*psi2(j) - psi1(j)*psi2(i))^2;
        probArray(j,i) = pn;
    end
end 
subplot(2,1,1);
imagesc(probArray)
%%
posArray = zeros(nbar+1);
Nx = np;
gapWidth = (Nx -nbar) /(nbar+1);
if not(gapWidth == round(gapWidth))
    disp('had to round gapWidth')
    gapWidth = round(gapWidth);
end



for i =1:(nbar+1)
    for j = 1:i
        thisArr = probArray(((i-1)*gapWidth+i): (i*gapWidth), ((j-1)*gapWidth+j): (j*gapWidth));
        posArray(i,j) = sum(sum(thisArr));
    end
end

norm = sum(sum(posArray));
posArray = posArray/norm;
subplot(2,1,2)
imagesc(posArray);
end
