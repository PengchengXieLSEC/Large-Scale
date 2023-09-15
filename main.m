%Codes of a model-based method for solving large-scale DFO
%Copyright: Pengcheng Xie & Ya-xiang Yuan 
%Connect: xpc@lsec.cc.ac.cn

global NF;
NF = 0;
[x, f, tp, N] = main_algorithm();
x
% norm(x - [1; 2; 3; 4; 5], 2) < 1e-3
f
NF
tp
N