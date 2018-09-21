%example code

%% imaging

v=rand(1024,1024);
imagesc(v)
colorbar
axis equal
axis xy

pcolor(v);shading interp

[ww,ee]=meshgrid(xv,xv);

subplot(2,1,1);imagesc(ww); subplot(2,1,2); imagesc(ee)
subplot(3,1,1); imagesc(ww.^2 + ee.^2); subplot(3,1,2);imagesc(ww); subplot(3,1,3); imagesc(ee)

%v=rand(1024,1024);
%imagesc(v)
%colorbar
%axis equal
%axis xy

%pcolor(v);shading interp

%% file operations

%% write
clear all
x = 0:.1:1;
A = [x; exp(x)];

fileID = fopen('exp.txt','w');
%fprintf(fileID,'%6s %12s\n','x','exp(x)');
fprintf(fileID,'%6.2f %12.8f\n',A);
fclose(fileID);

%% read
clear all
[a,b] = textread('exp.txt')



