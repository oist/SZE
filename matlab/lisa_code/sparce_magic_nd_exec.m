function[psi1, psi2, e1, e2] = sparce_magic_nd_exec(nbar, nx, L, n1, n2)

%% UNITS DISCLAIMER
% This code is set up to work in dimensionless units.
% You can work with inferior units by adding your human-centric units to
% the length, potential,...  and multiplying all Hamiltonian terms that
% go with 1/dx^2 by hbar^2/m [ NOT hbar^2/(2m) ].

%% BOSONS/FERMIONS DISCLAIMER
% If themultiple dimensions refer to multiple indistinguishable particles,
% then one needs to be careful before using the eigenstates, as the code
% will spit both bosonic (symmetric) and fermionic (antisymmetric) states
% (or superpositions of both if they are degenerate!). Therefore, one needs
% to check what is the symmetry of the states obtained 


%% PARAMETERS %%

number_of_states = 300; % number of eigenstates to compute

ndim = 1; %number of dimensions

%nbar = 2;
%nx = double(2^9);
Nx= (nbar+1)*nx + nbar;
(Nx -nbar) /(nbar+1);
% careful here on how your systems scales with ndim and Nx...
% the limit on an iMac is around Nx^ndim ~ 2^21

gridLength = L; % grid length (x goes from -length/2 to length/2)

% define position and momentum arrays
% notice that this is abc's symmetric x version
[x,dx,px,~] = fftdef_abc(gridLength/2.,Nx); 

%% TRIPLE WELL POTENTIAL %%
% distance = 9.; % trap separation
% 
% VL = 0.5*(x+distance).^2.; % left trap
% VM = 0.5*x.^2.; % middle trap
% VR = 0.5*(x-distance).^2.; % right trap
% 
% V_1D = min(min(VL, VM), VR);


%% DOUBLE WELL POTENTIAL %%
%distance = 5.; % trap separation

%VL = 0.5*(x+distance/2.).^2.; % left trap
%VR = 0.5*(x-distance/2.).^2.; % right trap

%V_1D = min(VL, VR);


%% HARMONIC POTENTIAL %%
% V_1D = 0.5*x.^2.;


%% DOUBLE SQUARE WELL POTENTIAL %%
 large_number = 1.e300;
 barrier_strength = 30.;
 V_1D = zeros(1,Nx);

 %insert the barriers! (evenly spaced!)
 for i=(1:nbar)
     V_1D(nx*i+i) = barrier_strength/dx;     
   
 end
 %V_1D = 0.5 * min((x+1.5).^2,(x-1.5).^2);


%% INTERACTION POTENTIAL 
g = 0;
%g = large_number; % interaction strength
V_int = eye(Nx) * g/dx;


%% DEFINING THE HAMILTONIAN AS A SPARSE MATRIX %%
switch ndim
    case 1
        Hamiltonian = kraken_1d(Nx, dx, V_1D, V_int);
    
    case 2
        Hamiltonian = kraken_2d(Nx, dx, V_1D, V_int);
    
    case 3
        Hamiltonian = kraken_3d(Nx, dx, V_1D, V_int);
    
    case 4
        Hamiltonian = kraken_4d(Nx, dx, V_1D, V_int);

    otherwise
        error('unkown dimension')

end


%% OPTION 1: FIND EIGENSTATES with energies around energy_around
%energy_around = 100.; % around which energy
%[states, energies] = eigs(Hamiltonian, number_of_states, energy_around);


%% OPTION 2: FIND EIGENSTATES with smallest energies
[states, energies] = eigs(Hamiltonian, number_of_states, 'sm');

%% FLIP THE VECTORS (by default Matlab sorts from highest to lowest eigenvalue)
energies = flipud(diag(energies));
states = flip(states,2);

%% obtain return wf

psi1 = states(:,n1);
psi2 = states(:,n2);
e1 = energies(n1);
e2 = energies(n2);

return

%% SHOW THE RESULTS

% plot the energies
figure(1)
plot(energies);
%number of states to plot
number_of_states = 10;
subplot_y = ceil(sqrt(number_of_states));
subplot_x = ceil(number_of_states/subplot_y);

switch ndim
    case 1
        %% SHOW THE 1D RESULTS
        % plot the eigenstates
        figure(2)
        for i = 1:number_of_states
            subplot(subplot_x, subplot_y, i)
            plot(x,abs(states(:,i)).^2)
        end

        %plot the eigenstates separated by their energies
        scale_factor = 4.;
        for i=1:number_of_states
            states_shift(:,i) = scale_factor * states(:,i) + energies(i);
        end

        figure(3)
        plot(x,states_shift)
    
    case 2 
        %% SHOW THE 2D RESULTS
        figure(2)
        for i = 1:number_of_states
            subplot(subplot_x, subplot_y, i)
            imagesc(x,x,abs(reshape(states(:,i),Nx,Nx)).^2)
        end

    case 3
        %% SHOW THE 3D RESULTS        
        figure(2)
        [xx,yy,zz] = meshgrid(x,x,x);
        for i = 1:number_of_states
            subplot(subplot_x, subplot_y, i)
            interp3(xx,yy,zz,abs(reshape(states(:,i),Nx,Nx,Nx)).^2,x*0,x*0,x*0)
        end
    
    case 4
        %% SHOW THE 4D RESULTS

    otherwise
        error('unkown dimension')

end

end
