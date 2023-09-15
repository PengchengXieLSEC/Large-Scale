%Codes of a model-based method for solving large-scale DFO
%Copyright: Pengcheng Xie & Ya-xiang Yuan 
%Connect: xpc@lsec.cc.ac.cn

function result = max_abs(f, Delta)
% Compute the maximum absolute value of f(alpha, beta) within the disk
% of radius Delta centered at the origin

N = 100; % Number of grid points in each dimension

alpha = linspace(-Delta, Delta, N);
beta = linspace(-Delta, Delta, N);

result = -inf;

for i = 1:N
    for j = 1:N
        alpha_val = alpha(i);
        beta_val = beta(j);
        
        f_val = abs(f(alpha_val, beta_val));
        if f_val > result && alpha_val^2 + beta_val^2 <= Delta^2
            result = f_val;
        end
    end
end
end