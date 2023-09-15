%Codes of a model-based method for solving large-scale DFO
%Copyright: Pengcheng Xie & Ya-xiang Yuan 
%Connect: xpc@lsec.cc.ac.cn

function [Qk] = Qkf(xbase, currx, dd1, dd2, a_k, b_k, c_k, d_k, e_k, f)
alpha = dot(currx - xbase, dd1);
beta = dot(currx - xbase, dd2);
Qk = f(xbase) + a_k * alpha + b_k * alpha ^ 2+...
    c_k * beta + d_k * beta ^ 2 + e_k * alpha * beta;
end