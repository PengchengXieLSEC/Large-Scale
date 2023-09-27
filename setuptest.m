%Codes of a model-based method for solving large-scale DFO
%Copyright: Pengcheng Xie & Ya-xiang Yuan 
%Connect: xpc@lsec.cc.ac.cn



function [xb, rhobeg, fopt, info] = setuptest(fun, n)
    xb = zeros(n, 1);
    fopt = 1.0e308;
    info = 0;
    ful = lower(fun);

    if strcmp(strtrim(ful), 'dqrtic')
        xb(:)= 2.0e0;
        rhobeg = 1.0e0;
    else
        info = 1;
    end

end
