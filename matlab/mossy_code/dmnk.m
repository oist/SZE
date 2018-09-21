
%%%part of crpations program for the rspdm of the tonks gas

%psi - array with single particle wf at different x values

function [rspdm,nk,Angps]=dmnk(psi,dx,dkx,Nparticles,Npts)

%calculate density matrix from other m script
%create nxn matrix with zeros
rspdm=zeros(Npts);
for i=1:Npts
    %fill diagonal
    [rspdm(i,i)]=rspdmij(psi,Nparticles,dx,i,i);
    for j=i+1:Npts
        %fill off-diagonal elements
        [rspdm(i,j)]=rspdmij(psi,Nparticles,dx,i,j);
        rspdm(j,i)=conj(rspdm(i,j));
    end
end
%density matrix has been obtained


nk=Npts*fftshift(diag(fft((ifft(rspdm)).'))*dx^2/(2*pi));
nk=nk.';

uF=zeros(size(psi));
for j=1:Nparticles
    uF(j,:)=fftshift(abs(fft(psi(j,:))).^2)/Npts;
end
Angps=ones(1,Nparticles)*uF*dx/dkx;
