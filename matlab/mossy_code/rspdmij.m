
%%part of the croation code

function [rspdmij]=rspdmij(psi,Nparticles,dx,ii,jj)

%single matrix entry
rspdmij=0;

%need not be calculated, only jj>ii considered (complex conjugate)
if (jj<ii)
    return
end

%diagonal elements
if (ii==jj)
    P=eye(Nparticles);
    %makes P matrix the size of Nparticles
else
    %define wave functions as needed from given psi
    psi_pom0=psi(1:Nparticles,ii:jj); %Determines size of matrix
    psi_pom1=psi_pom0;
    %psi_pom(1:Nparticles,ii:jj)=-psi_pom(1:Nparticles,ii:jj);
    psi_pom1(1:Nparticles,1)=0.5*psi_pom1(1:Nparticles,1);
    psi_pom1(1:Nparticles,jj-ii)=0.5*psi_pom1(1:Nparticles,jj-ii);
    %make matrix as density matrix as given for P
    P=eye(Nparticles)-2*conj(psi_pom0)*(psi_pom1).'*dx;
end
 
A=(inv(P).')*det(P);
%rspdmij=(psi(:,ii)')*psi(:,jj); %column indexing
rspdmij=(psi(:,jj).')*A*conj(psi(:,ii));
