%Codes of a model-based method for solving large-scale DFO
%Copyright: Pengcheng Xie & Ya-xiang Yuan 
%Connect: xpc@lsec.cc.ac.cn


function [x, f, alg, func] = twoDMoSub_4test(func, x_0, ~, ~)
alg='twoDMoSub';
  try
    [x, f] = twoDMoSub(func, x_0, false);
  catch
    x = nan;
    f = Inf;
  end
  
end

