function b = isfeasible(x, presbounds, preslineq, lb, ub, Aeq, beq, Aineq, lbineq, ubineq, tol_feas)
% Checks the feasibility of a point x with respect to bounds and constraints.

% Inputs:
%   - x: The point to check for feasibility.
%   - presbounds: Indicates whether there are bounds present.
%   - preslineq: Indicates whether there are linear inequalities present.
%   - lb: Lower bounds on x.
%   - ub: Upper bounds on x.
%   - Aeq: Matrix defining equality constraints.
%   - beq: Vector defining equality constraints.
%   - Aineq: Matrix defining inequality constraints.
%   - lbineq: Lower bounds on inequality constraints.
%   - ubineq: Upper bounds on inequality constraints.
%   - tol_feas: Feasibility tolerance.

% Output:
%   - b: A boolean value indicating whether the point x is feasible.

if presbounds
    % Check if x satisfies the bounds
    xl = (min(x - lb) > -tol_feas);
    xu = (max(x - ub) < tol_feas);
else
    xl = 1;
    xu = 1;
end

if ~isempty(Aeq)
    % Check if x satisfies the equality constraints
    Ax = (norm(Aeq * x - beq, Inf) < tol_feas);
else
    Ax = 1;
end

if preslineq
    % Check if x satisfies the linear inequalities
    xli = (min(Aineq * x - lbineq) > -tol_feas);
    xui = (max(Aineq * x - ubineq) < tol_feas);
    Aix = (xli && xui);
else
    Aix = 1;
end

% Check if x is feasible by satisfying all conditions
b = xl && xu && Ax && Aix;

