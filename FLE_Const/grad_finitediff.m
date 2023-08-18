function [g, func_eval] = grad_finitediff(fun, x, f, h, W)
% Computes the gradient of a function using finite differences.

% Inputs:
%   - fun: The function handle representing the objective function.
%   - x: The current point at which to evaluate the gradient.
%   - f: The current function value at point x.
%   - h: The step size for finite differences.
%   - W: The matrix of directions for finite differences.

% Outputs:
%   - g: The computed gradient vector.
%   - func_eval: The number of function evaluations performed.

n = size(W, 2); % Number of directions
g = zeros(n, 1); % Initialize the gradient vector
I = eye(n); % Identity matrix

for i = 1:n
    try
        % Evaluate the function 
        fval = feval(fun, x + h * W * I(:, i), 0);
    catch
        warning('grad approx has trouble computing f. exiting program.');
        return
    end
    g(i) = (fval - f) / h; % Compute the finite difference approximation of the gradient
end

func_eval = n; % Set the number of function evaluations to n (since n directions were used)

