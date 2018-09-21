function[stateArr] = getElevelsTP( maxlev, degenerate)
nstat = 0;
for i=1:maxlev
    for j=1:i
        nstat = nstat + j;
    end
end


stateArr = zeros(nstat, 3);
count= 1;

if degenerate
    for i=1:maxlev
        for j = 1:i
            for k=1:j
                stateArr(count,:) = [k,j,i];
                count = count+1;
            end

        end
    end
else
    for i=1:maxlev
        for j = 1:(i-1)
            for k=1:(j-1)
                stateArr(count,:) = [k,j,i];
                count = count+1;
            end
        end
    end
end
stateArr( all(~stateArr,2), : ) = [];
end