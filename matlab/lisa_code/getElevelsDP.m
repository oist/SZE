function[stateArr] = getElevelsDP(maxlev, degenerate)
nstat = maxlev^2- maxlev*(maxlev+1)/2;
stateArr = zeros(nstat, 2);
count= 1;

if degenerate
    for i=1:maxlev
        for j = 1:i
            stateArr(count,:) = [j,i];
            count = count+1;

        end
    end
else
    for i=1:maxlev
        for j = 1:(i-1)
            stateArr(count,:) = [j,i];
            count = count+1;

        end
    end
end

end
