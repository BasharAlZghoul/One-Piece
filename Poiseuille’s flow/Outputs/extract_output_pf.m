clear 
close all
clc

% inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%
Np = 384000;            % Number of colloids injected

case_num = 0;           % just for output name
IS = 6;                 % just for output name
velocity = 16;          % just for output name
colloid_dia = '_03ap_'; % just for output name
injection_type = 'U';   % just for output name

load('results.dat')     % Fortran code output
load('log_att.dat')     % Fortran code output
load('log_pass.dat')    % Fortran code output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


aux_att = find(log_att == 1);
aux_pass = find(log_pass == 1);
N_att = length(aux_att);
N_pass = length(aux_pass);
N_res = Np - N_att - N_pass;

attached = zeros(N_att,8);
pass = zeros(N_pass,8);
res = zeros(N_res,8);


domain = Np*8;
k = 1;
p = 1;
n = 1;
for i = 1: Np
    if (log_att(i) == 1)
    for j = 1:8
        c=(i-1)*8+j;
        
        attached(k,j) = results(c);         
    end
    k = k+1;
    end




    if (log_pass(i) == 1)
    for j = 1:8
        c=(i-1)*8+j;
        
       pass(p,j) = results(c);         
    end
    p = p+1;
    end

    if (log_pass(i) == 0 && log_att(i) == 0)
    for j = 1:8
        c=(i-1)*8+j;
        
       res(n,j) = results(c);         
    end
    n = n+1;
    end
end



% inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

also = 'IS_';
also2 = 'md_';

attached_variable_name = ['attached',colloid_dia, num2str(IS),also, num2str(velocity),also2,injection_type];
attached_variable_save = [attached_variable_name,'.mat'];
eval([attached_variable_name ' = attached;']);

pass_variable_name = ['pass',colloid_dia, num2str(IS),also, num2str(velocity),also2,injection_type];
pass_variable_save = [pass_variable_name,'.mat'];
eval([pass_variable_name ' = pass;']);

save(attached_variable_save,attached_variable_name)
save(pass_variable_save,pass_variable_name)

