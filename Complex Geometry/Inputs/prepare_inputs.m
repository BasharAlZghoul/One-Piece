clear
clc
close all

% INPUTS
% -----------------------------------------------------------
xc = 501; % number of cells in x direction
yc = 501; % number of cells in y direction

% velocity field
load('Ufx.mat') % Should be Ufx(x,y) NOT Ufx(y,x)
load('Ufy.mat') % Should be Ufy(x,y) NOT Ufy(y,x)

% geometry
% column1: x location % column2: y location % column3: grain radius
load('D.mat')  

% colloid initial location
load('xp.mat')  % colloids initial x position
load('yp.mat')  % colloids initial y position

% ------------------------------------------------------------



U(1:yc,1:xc)=Ufx';
V(1:yc,1:xc)=Ufy';


fileID2 = fopen('Ufx.txt','w');
fileID3 = fopen('Ufy.txt','w');

for i=1:yc
    for j=1:xc
        
        N=(j-1)*yc+i;        
        fprintf(fileID2,'%10.8f ',U(N));
        fprintf(fileID3,'%10.8f ',V(N));
     
    end      
end


fclose(fileID2);
fclose(fileID3);



fileID4 = fopen('D.txt','w');
fprintf(fileID4,'%3d',length(D));
fprintf(fileID4,'\n');

D = D';

for i = 1:3
    for j = 1:length(D(1,:))
                fprintf(fileID4,'%10.8f ',D(i,j));
    end
end

fclose(fileID4);




fileID10 = fopen('yp.txt','w');
for j = 1:length(yp)
                fprintf(fileID10,'%10.8f ',yp(j));
end
fclose(fileID10);

fileID11 = fopen('xp.txt','w');
for j = 1:length(xp)
                fprintf(fileID10,'%10.8f ',xp(j));
end
fclose(fileID11);



% histogram(yp,2000)