function[z, ez]= getZ(energies, beta)
%returns partition function as sum over all given energy levels of the
%complete system in form of an array

expfunc =  exp(-beta*energies); 
forexpval = energies.*exp(-beta*energies);
z = sum(expfunc);
ez = sum(forexpval)/z; 