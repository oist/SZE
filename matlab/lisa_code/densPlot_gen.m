function[] = densPlot_gen(psiA, psiB, nbar, hbar)
figure(1)
subplot(2,2,1);
plot(psiA);

subplot(2,2,3);
plot(psiB);

if not(length(psiA)== length(psiB))
    disp('WF have different lenghts!')
end
np = length(psiA);
probArray = zeros(np);

for i = 1:np
    for j = 1:np
        probArray(j,i) = abs(psiA(i)*psiB(j) - psiA(j)*psiB(i))^2;

    end
end

subplot(2,2,[2,4]); imagesc(probArray); hold on;
return
cd 'Documents/OIST/sze_out/matlab/'
filename = char(strcat('densPlot_gen_nbar_',string(nbar),'_hbar_',string(hbar), '.png'))
saveas(figure(1), filename);
cd ../../../..
close(figure(1))

end