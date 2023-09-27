function [alpha, fval, func, p, xn] = armijo_proj(fc, f, x, H, g, ds_step, pr, Z, unconst)
% Performs the Armijo backtracking line search with projection onto the feasible set.
%
% Inputs:
%   - fc: The current function value at the current point x.
%   - f: The objective function.
%   - x: The current point.
%   - H: The approximation of the inverse of the Hessian matrix.
%   - g: The gradient at the current point.
%   - ds_step: step-size for Low-Eval iteration. If not used in FLE, choose ds_step < 0.
%   - pr: The problem structure containing constraints and bounds.
%   - Z: The null space of the equality constraints.
%
% Outputs:
%   - alpha: The computed step size.
%   - fval: The number of function evaluations.
%   - func: The function value at the new point xn. 
%   - p: The search direction.
%   - xn: The new point.
%
% Author: Oumaima Sohab

%Initialization 
alpha = 1;
fval = 0;
eta = 1e-8;
tau = 0.5;
func = fc;

% Projection on the feasible set
if unconst
    xbar = x - H*g;
else
    xbar = project_on_feasible(x - Z*H*g, pr);
end

while (1)
    p = xbar - x;
    xn = x + alpha * (xbar - x);
    
    
    % Stop if the step size is too small
    if (alpha <= 1e-16)
        alpha = 0;
        func = fc;
        return
    end
    
    % Switch condition
    if (ds_step > 0 && alpha < min(1e-2, 1e-3*ds_step))
        alpha = 0;
        func = fc;
        return
    end
    
    func = f(xn, 0);
    fval = fval + 1;
    
    % Armijo condition
    if (func <= fc + eta * alpha * (Z * g)' * (xbar - x))
        return
    end
    
    alpha = tau * alpha; % Backtrack the stepsize
end

