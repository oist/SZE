energies = eS;
states = psiSame;
number_of_states = length(energies);
symAr = zeros(1,number_of_states);
for i=1:length(energies)
    symAr(i) = checkSymmetry2D(reshape(states(:,i),Nx,Nx));       
end
bInd = find(symAr);
nbState = length(bInd);
stateSize = size(states);
bstates = zeros(stateSize(1), nbState);
benergies = zeros(1, nbState);
for i = 1:nbState
    bstates(:,i) = states(:,bInd(i));
    benergies(i) = energies(bInd(i));
end


%% obtain return wf

figure(2)
for i = 1:nbState
    subplot(subplot_x, subplot_y, i)
    imagesc(x,x,(reshape(bstates(:,i),Nx,Nx)))
end
