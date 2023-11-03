clear 
close all
clc


% inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%
case_num = 0;            % just for output name
colloid_dia = '_2ap_';   % just for output name
injection_type = 'FW_';  % just for output name

Np = 192000;             % Number of colloids injected

load('results.dat')      % Fortran code output
load('log_att.dat')      % Fortran code output
load('log_pass.dat')     % Fortran code output
load('D.mat')            % geometry

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



D = D';
aux_att = find(log_att == 1);
aux_pass = find(log_pass == 1);
N_att = length(aux_att);
N_pass = length(aux_pass);
N_res = Np - N_att - N_pass;
r = 23;
attached = zeros(N_att,r);
pass = zeros(N_pass,r);
res = zeros(N_res,r);

domain = Np*r;
k = 1;
p = 1;
n = 1;
for i = 1: Np
    if (log_att(i) == 1)
    for j = 1:r
        c=(i-1)*r+j;
        
        attached(k,j) = results(c);         
    end
    k = k+1;
    end




    if (log_pass(i) == 1)
    for j = 1:r
        c=(i-1)*r+j;
        
       pass(p,j) = results(c);         
    end
    p = p+1;
    end

    if (log_pass(i) == 0 && log_att(i) == 0)
    for j = 1:r
        c=(i-1)*r+j;
        
       res(n,j) = results(c);         
    end
    n = n+1;
    end
end




attached_variable_name = ['attached',colloid_dia,injection_type, num2str(case_num)];
attached_variable_save = [attached_variable_name,'.mat'];
eval([attached_variable_name ' = attached;']);

pass_variable_name = ['pass',colloid_dia,injection_type, num2str(case_num)];
pass_variable_save = [pass_variable_name,'.mat'];
eval([pass_variable_name ' = pass;']);

save(attached_variable_save,attached_variable_name)
save(pass_variable_save,pass_variable_name)





