%Codes of a model-based method for solving large-scale DFO
%Copyright: Pengcheng Xie & Ya-xiang Yuan 
%Connect: xpc@lsec.cc.ac.cn

function value = func(x, dd1, dd2, xbase)
% Define the objective function here
alpha = dot(x - xbase, dd1);
beta = dot(x - xbase, dd2);
value = 9 - 2* alpha + alpha ^ 2 + 4 * beta + beta ^ 2;
end