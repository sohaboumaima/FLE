function [x, f, g, H, beta, success, func_eval] = fullevalstep(fun, x, f, H, x_old, g_old, h, Z, numfe, pr, alpha)
% Full-Eval step which corresponds to projected BFGS step
%
% Inputs:
%   - fun: Function handle of the objective function. The function should accept two arguments:
%           - x: Vector representing the current point.
%           - flag: A flag indicating whether to compute the function value (flag = 0) or the true function
%             value (flag = 1).
%   - x: Current iterate.
%   - f: Function value at x.
%   - H: Hessian matrix (approximation) at x.
%   - x_old: Previous iterate.
%   - g_old: Gradient at x_old.
%   - h: Finite difference parameter.
%   - Z: Null space matrix (for equality constraints).
%   - numfe: Number of iterations performed so far.
%   - pr: Structure containing the problem definition (constraints)
%   - alpha: Step-size for Low-Eval. It should be given a negative value (i.e -1) if .
%
% Outputs:
%   - x: Updated point after the step.
%   - f: Updated function value at the new point.
%   - g: Gradient at the previous point.
%   - H: Updated Hessian matrix.
%   - beta: Step size for the line search.
%   - success: Success flag indicating if the step was successful (1) or not (0).
%   - func_eval: Total number of function evaluations performed after the step.
%
%
% Author: Oumaima Sohab

    % Compute the gradient using finite differences
    [g, fe] = grad_finitediff(fun, x, f, h, Z);
    
    if isnan(norm(g))
        % Gradient evaluation failed
        beta = 0;
        success = 0;
        func_eval = fe;
        return
    end
    
    Im = eye(size(Z, 2));
    
    if numfe > 1
        % Compute the difference and product vectors for BFGS update
        s = Z' * (x - x_old);
        y = g - g_old;
        c = y' * s;
        r = 1 / c;
        
        if (numfe == 2) && (c > 0)
            H = (c / (y' * y)) * Im;
        end
        
        if (c > 1e-10 * norm(y) * norm(s))
            % Update the Hessian matrix using the BFGS formula
            H = (Im - r * s * y') * H * (Im - r * y * s') + r * (s * s');
        end
    end
    
    % Perform a line search using Armijo-backtracking
    [beta, fval, f, d, x] = armijo_proj(f, fun, x, H, g, alpha, pr, Z);
    func_eval = fe + fval;
    
    if (beta == 0 || norm(d) == 0)
        success = 0;
    else
        success = 1;
    end
end

