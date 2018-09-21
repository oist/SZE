function[pnRemAr] = getPnRem(tAr, npoints, lminP, L, npart, hb, m, estat)
     %npoints = 50; % 2*npoints will be calculated to determine pnrem
     lmin = lminP*L;
     
     lhalf = linspace(lmin, L/2, npoints);
     lhalf2 = lhalf + L/2 - lhalf(1);
     lvalues = [lhalf, lhalf2(2:end)];
     pnRemAr = zeros(npart+1, length(tAr));

     figure(7)
     for b =1:length(tAr)
        beta = 1/tAr(b);
        pnremAr = zeros(npart+1, length(lvalues));
        EAr = zeros(npart+1, length(lvalues));
        for p=1:length(lvalues)
            l = lvalues(p);
            Z= zeros(1, npart+1);
            E= zeros(1, npart+1);
            for i = 1:npart+1 %i particles on the left side
                %n p on left side = i-1
                for n =1:length(estat)
                    thisE = 1/2*((i-1)*getEnInf(estat(n,1), L-l,hb, m) + ... 
                    (npart+1-i)*getEnInf(estat(n,2), l, hb, m) + (i-1)*getEnInf(estat(n,2), L-l, hb, m) + ...
                    (npart+1-i)*getEnInf(estat(n,1), l,hb, m));
                    Z(i) = Z(i) + exp(-beta*thisE);
                    E(i) = E(i) + thisE*exp(-beta*thisE);
                end
            end
            Ztotal = sum(Z);
            pnremAr(:,p) = Z/Ztotal;
            EAr(:,p) = E./Z;
            EAr(isnan(EAr)) = 0;
        end
         for i = 1:npart+1
            plot(lvalues, pnremAr(i,:)); hold on; 
            xlabel('L');
            ylabel('pnrem(l)');
            pnRemAr(i, b) = max(pnremAr(i,:));
         end
     end

end