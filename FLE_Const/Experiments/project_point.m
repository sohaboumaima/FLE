function [x_proj] = project_point(x, lb,ub,Aeq,beq,Aineq,lbineq,ubineq)
% This function is used in case the initial point in unfeasible. Uses 
% quadprog to perform projection into the feasible set

A = [Aineq; -Aineq];

b = [ubineq; -lbineq];

x_proj = quadprog(eye(length(x)),-x,A,b,Aeq,beq,lb,ub);

end
