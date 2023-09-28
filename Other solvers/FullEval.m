function [xsol, fsol, output] = FullEval(fun, x0, pr, maxevals, unconst)
% Funtion implementing pure Full-Eval
%
% Input:
%   - fun: Function handle of the objective function. The function should accept two arguments:
%           - x: Vector representing the current point.
%           - flag: A flag indicating whether to compute the function value (flag = 0) or the true function
%             value (flag = 1).
%   - x0: Initial point.
%   - pr: Structure containing the problem definition (constraints).
%   - maxevals: Maximum number of function evaluations allowed.
%
% Output:
%   - xsol: Solution point found by the optimization algorithm.
%   - fsol: Function value at the solution point.
%   - output: Structure containing additional information:
%             - numIter: Total number of iterations performed.
%             - funcCount: Total number of function evaluations performed.
%             - fvalues: Array of function values at each iteration.
%             - funcCounts: Array of function evaluation counts at each iteration.
%             - truefunc: Array of true function values at each iteration (in case of noisy function).
%             - exitflag: Exit flag indicating the reason for termination:
%               - 0: Maximum number of function evaluations reached.
%               - 1: Maximum number of iterations reached.
%               - 2: Termination condition satisfied (gradient magnitude is below threshold).
%
%
% Author: Oumaima Sohab

% Initialization
n          = length(x0);            % Dimension of the problem
h          = sqrt(eps);             % Small step size
x          = x0;                    % Initial point
f          = fun(x, 0);             % Evaluate the function at the initial point
func_eval  = 1;                     % Count of function evaluations performed

maxiter    = 2000*n;                % Maximum number of iterations based on problem dimension
fvalues    = zeros(maxiter, 1);     % Array to store function values at each iteration
nfeval     = zeros(maxiter, 1);     % Array to store function evaluation counts at each iteration
truef      = zeros(maxiter, 1);     % Array to store true function values (in case of noisy function)
stepsize   = zeros(maxiter, 1);     % Array to store step sizes at each iteration

num_success_iter   = 0;             % Counter for successful iterations
num_unsuccess_iter = 0;             % Counter for unsuccessful iterations
numIter = 0;                        % Counter for total iterations

% Checking for equality constraints
Aeq    = pr.Aeq;                    % Equality constraint matrix
beq    = pr.beq;                    % Equality constraint vector

if ~isempty(Aeq)
    % Normalize the rows of Aeq and beq
    for i = 1:size(Aeq, 1)
        auxa = norm(Aeq(i, :));
        Aeq(i, :) = Aeq(i, :) / auxa;
        beq(i) = beq(i) / auxa;
    end
    
    Z = null(Aeq);                  % Null space of Aeq
    m = size(Z, 2);                 % Dimension of the null space
else
    Z = eye(n);                      % If no equality constraints, use identity matrix
    m = n;
end

H = eye(m);                         % Hessian matrix for linear inequalities
g_old = [];
x_old = x;

while m          % The problem is already at the solution or infeasible when m= 0
    numIter = numIter + 1;          % Update the number of iterations

    fvalues(numIter) = f;           % Store the current function value
    nfeval(numIter)  = func_eval;   % Store the current function evaluation count
    truef(numIter)   = fun(x, 1);    % Store the true function value
    
    % Check termination conditions
    if func_eval >= maxevals
        exitflag = 0;               % Maximum number of function evaluations reached
        break
    elseif numIter >= maxiter
        exitflag = 1;               % Maximum number of iterations reached
        break
    end

    % Perform a full eval step 
    [xn, fn, g, H, beta, success, nf] = fullevalstep(fun, x, f, H, x_old, g_old, h, Z, num_success_iter, pr, -1, unconst);
    if (numIter == 1)
        ng0 = norm(g);
    end
    func_eval = func_eval + nf;     % Update the function evaluation count
    stepsize(numIter) = beta;
    stepsize(numIter) = beta;


    % Update iterate and function value
    if ~success
        num_unsuccess_iter = num_unsuccess_iter + 1;
    else
        x_old = x;
        g_old = g;
        x     = xn;
        f     = fn;
        num_success_iter = num_success_iter + 1;
    end

    if norm(g) <= ng0 * 1e-12
        exitflag = 2;               % Termination condition satisfied (gradient magnitude is below threshold)
        break
    end
end

xsol = x;                           % Solution point found by the optimization algorithm
fsol = f;                           % Function value at the solution point

% Output structure
output.numIter = numIter;           % Total number of iterations performed
output.funcCount = func_eval;       % Total number of function evaluations performed
output.fvalues = fvalues(1:numIter);           % Array of function values at each iteration
output.funcCounts = nfeval(1:numIter);         % Array of function evaluation counts at each iteration
output.truefunc = truef(1:numIter);            % Array of true function values at each iteration (in case of noisy function)
output.flag = exitflag;             % Exit flag indicating the reason for termination
end
