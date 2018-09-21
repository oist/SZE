
%%%  fftdef.m --- defining the grids for the fft/split method
%%%
%%%  2000-2001, Thomas Busch    (thbusch@phys.ucc.ie)
%%%  Version 1.0
%%%
%%%  Change history:
%%%  14 Nov 01 Thomas Busch    - starting better documentation 
%%%
%%%  Notes: defining various grids for the split-operator/fft
%%%         method for solving the GPE
%%%
%%%         [X,DX,PX,DPX]=FFTDEF(XMAX,NGX) 
%%%         where XMAX is the maximum value for the X-coordinate 
%%%         and NGX the amount of points in that direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[x,dx,px,dpx]=fftdef_abc(xmax,Ngx);


% spacing in position space
dx = 2*xmax/(Ngx-1);

% maximum values
pxmax = pi/dx;

% spacing in momentum space
dpx = 2*pxmax/Ngx;

% grid vector in position space
x = (0:Ngx-1)'*dx-xmax;

% grid vector in momentum space
pxn = ((1:Ngx)*dpx-pxmax)';   

% reordination needed for the fourier transform
px(Ngx/2+2:Ngx,1) = pxn(1:Ngx/2-1,1);
px(1:Ngx/2+1,1) = pxn(Ngx/2:Ngx,1);





