%program for density of wf in nbar case
function[] = densPlot_noTun(nbar)
np = 2^8;
figure(1)
subplot(2,2,1);
psi1 = analyticNbar_noTun(nbar,2, np);

subplot(2,2,3);
psi2 = analyticNbar_noTun(nbar,4, np);

probArray = zeros(np);
%probArray2 = zeros(np);
%probArray3 = zeros(np);

for i = 1:np
    for j = 1:np
        probArray(j,i) = abs(psi1(i)*psi2(j) - psi1(j)*psi2(i))^2;
        %probArray2(i,j) = psi1(j)*psi2(i);    
        %probArray3(i,j) = psi1(i)*psi2(j);    
        %probArray(i,j) = abs(psi1(j)*psi2(i));
    end
end

subplot(2,2,[2,4]); imagesc(probArray); hold on;
cd 'Documents/OIST/sze_out/matlab/'
filename = char(strcat('densPlot_noTun_nbar',string(nbar), '.png'))
saveas(figure(1), filename);
cd ../../../..
close(figure(1))
%subplot(3,3,6);imagesc(probArray2); subplot(3,3,9); imagesc(probArray3)
end