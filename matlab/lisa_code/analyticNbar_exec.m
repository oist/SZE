%make single particle WF for well with nb barriers starting from infinite
%well
function [psiSplit]= analyticNbar_exec(nbar, n, np, psiI, xv)

    
    zeroCutOff = 1e-7;
    %np = 2^8;
    safetyDistance = 2;
    L = 10;
    xspacing = L/(np-1);
    %nbar =2;
    gamma = (nbar-1)/(nbar + 1);
    nparticles = 2;


    %WF for ininite well, nth state
    %separate well according to barrier position
    %find nodes in new wells
    %new WF: 1WF for each well, state = number of nodes, stick together
    %adjust weights of WF to have same energy in every well

    %ground state
    %n = 2;
    if psiI
        psiInf = psiI;
        xvalues = xv;
    else
        k = n*pi/L;
        [psiInf,xvalues] =  psiInfWell(L, k, np);   %array with psi values along well
        psiInf(abs(psiInf)<zeroCutOff) =0;
    end

    %find barrier positions
    %create array with gap boundaries
    gaps = zeros(nbar+1, 2);
    barPos = zeros(nbar+2,1);
    %psiGapInd = zeros(nbar+1, 2);

    %even barrier spacing
    gapWidth = L/(nbar+1);
    for i=0:(nbar)
        leftB = i*gapWidth;
        rightB = (i+1)*gapWidth;
        barPos(i+2) = rightB;
        gaps(i+1,:) = [leftB,rightB ];
       %psiGapInd(i+1,:) = [(leftB/xspacing +1), (rightB/xspacing +1)];
    end
    


    indexBar = zeros(nbar+1, 1);
    for i=1:(nbar+2)
        [c1, index1] = min(abs(xvalues-barPos(i)));
       indexBar(i) = index1;
    end

    %psiInfSplit = zeros(nbar+1,1);
    signPsi = zeros(nbar+1,1);

    for i = 1:(nbar+1)
        %psiInfSplit = psiInf(indexBar(i): indexBar(i+1))
        sgAr = psiInf((indexBar(i)+safetyDistance): (indexBar(i+1)-safetyDistance));
        if all(sgAr>=0)
            signPsi(i) = 1;
        elseif all(sgAr<=0)
            signPsi(i) = -1;
        else
            signPsi(i) =0;
        end
    end

    %create new psi based on signPsi
    psiSplit = zeros(length(psiInf),1);
    %enSplit = zeros(nbar+1,1);
    %assumption: max 1 node in gap
    enLow = energyInfWell(1*pi/gapWidth);
    enHigh = energyInfWell(2*pi/gapWidth);
    weight = enLow/enHigh;
    for i = 1:(nbar+1)
        lenThisPsi = length(psiSplit(indexBar(i):indexBar(i+1)));
        if abs(signPsi(i)) ==1
            k = 1*pi/gapWidth;
            thisPsi = signPsi(i)*psiInfWell(gapWidth, k, lenThisPsi);
        elseif signPsi(i) ==0
            if i>1
                disp('i greater 1')
                thisSign = sign(psiSplit(indexBar(i)-1));
                j = 0;
                while thisSign ==0
                    j=j-1;
                    thisSign = sign(psiSplit(indexBar(i)-j));
                end
            else
                thisSign = 1;
            end
            k = 2*pi/gapWidth;
            thisPsi = weight*thisSign*psiInfWell(gapWidth, k, lenThisPsi);
        else
            disp('something went wrong with signPsi!')
        end
        %enSplit(i) = energyInfWell(k);
        psiSplit(indexBar(i):indexBar(i+1)) = thisPsi;
    end
    psiSplit = psiSplit/norm(psiSplit);
    plot(xvalues, psiSplit);
end
    %normalize the WF(necessary?)
    %determine weight of state in each gap by making energy in each gap equal


    



