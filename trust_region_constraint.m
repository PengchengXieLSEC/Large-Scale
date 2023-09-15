%Codes of a model-based method for solving large-scale DFO
%Copyright: Pengcheng Xie & Ya-xiang Yuan 
%Connect: xpc@lsec.cc.ac.cn

function [c, ceq] = trust_region_constraint(x, dd1, dd2, xbase, Delta)
% Trust-region constraint
alpha = dot(x - xbase, dd1);
beta = dot(x - xbase, dd2);
c = norm([alpha, beta], 2) - Delta;
ceq = [];
end