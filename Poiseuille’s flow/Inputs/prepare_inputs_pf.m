clear 
close all
clc


load('yp.mat')

fileID = fopen('yp.txt','w');

for j = 1:length(yp)
 %               N=(j-1)*3+i;
                fprintf(fileID,'%10.8f ',yp(j));
     
  
end
fclose(fileID);