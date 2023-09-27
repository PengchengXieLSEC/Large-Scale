%Codes of a model-based method for solving large-scale DFO
%Copyright: Pengcheng Xie & Ya-xiang Yuan 
%Connect: xpc@lsec.cc.ac.cn



function [fun, info] = evalfun(name, x, n)
    %% TestProblemF
    
    INFINITY = 1.0e308;
    f = INFINITY;
    info = 0; % info = 0 means successful evaluation
    y = zeros(n + 2, 1);
    y(0 + 1) = 0.0e0;
    y(n + 1 + 1) = 0.0e0;

    for i = 1:n
        y(i + 1) = x(i);
    end

   
    name = lower(name);
    %%
    if (strcmp(name, 'dqrtic'))
    
        f = 0.0e0;

        for i = 1:n
            f = f + (y(i + 1) - real(i)).^4;
        end
    else
        % info = 1 means unknown function name
        info = 1;
    end

    fun = f;

end
