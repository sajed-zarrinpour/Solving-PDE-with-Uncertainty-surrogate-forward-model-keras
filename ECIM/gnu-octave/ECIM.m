clear;clc;
dataFileNameEx = '.csv';
dataFileNamePre = 'data';
files = dir('*.csv');
fileName = [dataFileNamePre,'-',num2str(length(files)+1),dataFileNameEx];

% Solving NLSE in 1D on [0,1]: 
domain = [0,1];
number_of_samples = 5000;
maxIteration = 500;
% creating mesh
h = 0.125;
mesh = domain(1):h:domain(end);
mesh_length = length(mesh);
%calculating the unformly distributed uncertainty on the mesh
lambda0 = 0.3;
lambda1 = 1.5;

xi = 1;
epsilon =0.0001;
tic
% res = [];
for k=1:number_of_samples
    %generating 
    a = (lambda1 - lambda0).*rand(mesh_length,1) + lambda0;
    % Constructing L, b 
    L_center = zeros(mesh_length,mesh_length);
    L_left_patch = zeros(mesh_length, 1); 
    L_right_patch = zeros(mesh_length, 1);
    
    b = zeros(mesh_length,1);
     for i=1:mesh_length
        if i == 1
            L_center(i,i)   = L_center(i,i) + (2*a(i)+a(end)+a(i+1))/2;
            L_left_patch(i ,1) = - (a(i+1)+a(i))/2;
            L_center(i,i+1) = L_center(i,i) - (a(end)+a(i))/2;

            b(i) = xi * ((a(i+1)-a(i))/2); 
        elseif i == mesh_length
            L_center(i,i)   = L_center(i,i) + (2*a(i)+a(i-1)+a(1))/2;
            L_center(i,i-1) = L_center(i,i) - (a(1)+a(i))/2;
            L_right_patch(end,1) = - (a(i-1)+a(i))/2;

            b(i) = xi * ((a(1)-a(i))/2);
        else
            L_center(i,i)   = L_center(i,i) + (2*a(i)+a(i-1)+a(i+1))/2;
            L_center(i,i-1) = L_center(i,i) - (a(i+1)+a(i))/2;
            L_center(i,i+1) = L_center(i,i) - (a(i-1)+a(i))/2;

            b(i) = xi * ((a(i+1)-a(i))/2);
        end
     end
    L = [L_left_patch L_center L_right_patch]; % L is n+2 * n
    L = L .* (1/(h^2));
    b = b .* (1/h);
    
    u = zeros(mesh_length,1);
    A_ex = A_exact(a);
    
    u_patched = [u(length(u)); u ; u(1)]; % n+2 *1
    A_app = A_eff(u_patched,u.',a,L,b,h,1);
    
    %error_in_loop = [];
    min_error = 1;
    current_error = 1;
    v = u;
    
    iter = 1;
    while (abs(A_app - A_ex) > epsilon) && iter<=maxIteration
        current_error = abs(A_app - A_ex);
        if min_error > current_error
            v = u;
            min_error = current_error;
        end
        %error_in_loop = [error_in_loop  , current_error];
        
        f = @(x)(L*x-b);
        u_patched = fsolve(f,u_patched);
        clc
        u = u_patched(2:length(u_patched)-1);
        

        u_patched = [u(length(u)); u ; u(1)]; % n+2 *1
        
        A_app = A_eff(u_patched,u.',a,L,b,h,1); 
        %
        iter = iter + 1;
    end

    v_patched = [v(length(v)); v ; v(1)];
    A_app = A_eff(v_patched,v.',a,L,b,h,1); 
    res = [a.', v.', A_app,A_ex, min_error]
    res_str = '';
    for r = 1:length(res)-1
        res_str = [res_str, num2str(res(r)),','];
    end
    res_str = [res_str, num2str(res(end))];
    system(['echo ', res_str, '>>',fileName]);
%     res = [res;[a.', v.', A_app,A_ex, min_error]];
    %plot(error_in_loop)
end
toc
4
% csvwrite('data.csv',res)