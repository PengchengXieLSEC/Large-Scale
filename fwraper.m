%Codes of a model-based method for solving large-scale DFO
%Copyright: Pengcheng Xie & Ya-xiang Yuan 
%Connect: xpc@lsec.cc.ac.cn



function value = fwraper(x, func)
% Define the objective function here
global NF;
NF = NF + 1;
value = func(x);
end