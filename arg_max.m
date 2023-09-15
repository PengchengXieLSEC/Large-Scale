%Codes of a model-based method for solving large-scale DFO
%Copyright: Pengcheng Xie & Ya-xiang Yuan 
%Connect: xpc@lsec.cc.ac.cn


function [alpha_star, beta_star] = arg_max(f, Delta)
    % Compute the argument that maximizes f(alpha, beta) within the disk
    % of radius Delta centered at the origin
    
    N = 100; % Number of grid points in each dimension
    
    alpha = linspace(-Delta, Delta, N);
    beta = linspace(-Delta, Delta, N);
    
    max_val = -inf;
    alpha_star = 0;
    beta_star = 0;
    
    for i = 1:N
        for j = 1:N
            alpha_val = alpha(i);
            beta_val = beta(j);
            
            f_val = f(alpha_val, beta_val);
            if f_val > max_val && alpha_val^2 + beta_val^2 <= Delta^2
                max_val = f_val;
                alpha_star = alpha_val;
                beta_star = beta_val;
            end
        end
    end
end