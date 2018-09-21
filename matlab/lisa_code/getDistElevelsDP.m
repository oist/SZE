function[stateArr] = getDistElevelsDP(maxlev)
nstat = maxlev^2- maxlev*(maxlev+1)/2;
stateArr = zeros(nstat, 2);
count= 1;
for i=1:maxlev
    for j = 1:maxlev
        stateArr(count,:) = [j,i];
        count = count+1;
    end
end

end
