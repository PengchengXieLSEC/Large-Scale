%Codes of a model-based method for solving large-scale DFO
%Copyright: Pengcheng Xie & Ya-xiang Yuan 
%Connect: xpc@lsec.cc.ac.cn

function value = ff(x, n)
% Define the objective function here
value = 0;
for i=1:n-1
 for i=1:n
     value = value + i * (x(i, :) - i) .^ 2;
 end
end
end




