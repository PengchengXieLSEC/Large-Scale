%Codes of a model-based method for solving large-scale DFO
%Copyright: Pengcheng Xie & Ya-xiang Yuan 
%Connect: xpc@lsec.cc.ac.cn

function [x] = trust_region_subproblem(Qk, Delta, dd1, dd2, xbase)
% Solve the trust-region subproblem
nonlcon = @(x) trust_region_constraint(x, dd1, dd2, xbase, Delta);
options = optimoptions('fmincon', 'Display', 'off', 'OptimalityTolerance', 1e-16, 'ConstraintTolerance', 1e-16, 'StepTolerance',1e-16);
x = fmincon(Qk, xbase, [], [], [], [], [], [], nonlcon, options);
end