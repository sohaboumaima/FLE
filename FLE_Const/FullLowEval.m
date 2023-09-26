function [xsol, fsol, output] = FullLowEval(fun, x0, alpha0, dirgen, pr, maxevals, tol_feas)

% Code for Full-Low Evaluation methods for bound linearly constrained optimization.
%
% Inputs:
%   - fun: Objective function to be minimized. It should take the form f = fun(x, mode),
%          where x is the input vector and mode is a flag indicating the evaluation mode 
%          (0 for function value, 1 for true function value).
%   - x0: Initial iterate.
%   - alpha0: Initial step size for Low-Eval.
% 	- dirgen: Method that is used for generating the directions with
%           0: Deterministic variant - always choose the generators of an 
%           approximate tangent cone (but randomly ordered)
%           1: Proceed with a random subset selected within the 
%           generators of an approximate tangent cone
%           2: Attempt to benefit from linear subspaces included in 
%           an approximate tangent cone to reduce the size of the (random) 
%           polling set, otherwise draw directions similarly to 1
%   - pr: Structure containing problem information including bounds and constraints.
%   - maxevals: Maximum number of function evaluations allowed.
%   - tol_feas: Tolerance for feasibility in the LE step.
%
% Outputs:
%   - xsol: Solution point.
%   - fsol: Objective function value at the solution point.
%   - output: Structure containing additional output information.
%       - numIter: Number of iterations performed.
%       - funcCount: Total number of function evaluations.
%       - fvalues: Array of function values at each iteration.
%       - funcCounts: Array of function evaluation counts at each iteration.
%       - truefunc: Array of true function values at each iteration.
%       - numsfe: Number of successful FE steps.
%       - numusfe: Number of unsuccessful FE steps.
%       - numsle: Number of successful LE steps.
%       - numusle: Number of unsuccessful LE steps.
%       - flag: Termination flag indicating the reason for termination
%           0: Max evaluations reached
%           1: Max iterations reached
%           2: Gradient magnitude below threshold
%           3: Low-Eval step-size below tolerance
%
% Author : Oumaima Sohab
   
n          = length(x0);                          % Dimension of the problem
h          = sqrt(eps);                           % Smallest positive floating-point number
evl        = "FE";                                % Initial evaluation type (Full Eval)
x          = x0;                                  % Initial solution
f          = fun(x, 0);                           % Initial function value
func_eval  = 1;                                   % Function evaluation count

% Initialize arrays to store optimization information
maxiter    = maxevals;                            % Maximum number of iterations (same as maxevals)
fvalues    = zeros(maxiter, 1);                   % Array to store function values at each iteration
nfeval     = zeros(maxiter, 1);                   % Array to store function evaluation counts at each iteration
truef      = zeros(maxiter, 1);                   % Array to store true function values (in case of noisy function)
betas      = zeros(maxiter, 1);                   % Array to store step sizes at each iteration
alphas     = zeros(maxiter, 1);                   % Array to store step sizes at each iteration

% Constants and parameters for the optimization algorithm
theta      = 0.5;
gamma      = 2;
In         = eye(n);
alpha      = alpha0;                              % Initial step size
tol_alpha  = 1e-10;                               % Threshold for step size termination
p0         = log(theta)/log(theta/gamma);
nbunc      = ceil(log(1-log(theta)/log(gamma))/log(2)); % Minimum number of directions required in the unconstrained case

% Initialize counters for successful/unsuccessful iterations in both evaluation steps
num_success_iter_le   = 0;
num_unsuccess_iter_le = 0;
num_success_iter_fe   = 0;
num_unsuccess_iter_fe = 0;
numIter = 0;
nule    = 0;

% Determine the presence of bounds
lowbnd = find(pr.lb > -1e20, 1);
upbnd  = find(pr.ub < 1e20, 1);16
presbounds = ~isempty(upbnd) || ~isempty(lowbnd);

% Get the constraints from the structure pr
Aeq    = pr.Aeq;
beq    = pr.beq;
Aineq  = pr.Aineq;
lbineq = pr.lbineq;
ubineq = pr.ubineq;
lb     = pr.lb;
ub     = pr.ub;

if ~isempty(Aeq)
    presleq = 1;
    for i=1:size(Aeq,1)
        auxa = norm(Aeq(i,:));
        Aeq(i,:) = Aeq(i,:)/auxa;
        beq(i) = beq(i)/auxa;
    end
    Z = null(Aeq);
    m = size(Z,2);
else
    presleq = 0;
    Z = In;
    m = n;
end

% Gathering linear inequalities and bounds
if ~isempty(Aineq)
    preslineq = 1;
    for i=1:size(Aineq,1)
        auxa = norm(Aineq(i,:));
        Aineq(i,:) = Aineq(i,:)/auxa;
        lbineq(i) = lbineq(i)/auxa;
        ubineq(i) = ubineq(i)/auxa;
    end
    Ai = [Aineq;In];
    li = [lbineq;lb];
    ui = [ubineq;ub];
    ni = length(ui);
else
    preslineq = 0;
    Ai = In;
    li = lb;
    ui = ub;
    ni = n;
end

% Calculate the projected normal vectors
W = Z*Z'*Ai';
Im = eye(m);
H = Im;
numfe = 0;
g_old = [];
x_old = x;

while m  % The problem is already at the solution or infeasible when m = 0

    numIter = numIter + 1;                          % Update the number of iterations

    fvalues(numIter) = f;                           % Store the current function value
    nfeval(numIter)  = func_eval;                   % Store the current function evaluation count
    truef(numIter)   = fun(x, 1);                    % Store the true function value
    
    % Check termination conditions
    if func_eval >= maxevals
        exitflag = 0;                               % Maximum number of function evaluations reached
        break
    elseif numIter >= maxiter
        exitflag = 1;                               % Maximum number of iterations reached
        break
    end
    
    switch evl
        case "FE"
            % Full eval step
            [xn, fn, g, H, beta, success, nf] = fullevalstep(fun, x, f, H, x_old, g_old, h, Z, numfe, pr, alpha);
            if (numIter == 1)
                ng0 = norm(g);
            end

            func_eval = func_eval + nf;
            betas(numIter) = beta;
            
            % Update iterate and function value
            if ~success
                evl = "LE";
                num_unsuccess_iter_fe = num_unsuccess_iter_fe + 1;
            else
                x_old = x;
                g_old = g;
                x     = xn;
                f     = fn;
                num_success_iter_fe = num_success_iter_fe + 1;
                numfe = numfe + 1;
            end
            if norm(g) <= ng0 * 1e-12
                exitflag = 2;                           % Termination condition satisfied (gradient magnitude is below threshold)
                break
            end

        case "LE"
            % Low eval step
            [x, f, alpha, success, nf] = lowevalstep(fun, x, f, alpha, pr, ni, li, ui, Ai, Z, W, dirgen, preslineq, presbounds, presleq, m, p0, nbunc, tol_feas);
            
            func_eval       = func_eval + nf;
            alphas(numIter) = alpha;
            
            if ~success
                nule = nule + 1;
                num_unsuccess_iter_le = num_unsuccess_iter_le + 1;
            else
                num_success_iter_le   = num_success_iter_le + 1;
            end
            if (nule > max(2,nf - 1))
                nule = 0;
                evl = "FE";
            end

            if (alpha <= tol_alpha)
                exitflag = 3;                           % Termination condition satisfied (step size is below threshold)
                break
            end
            
    end
end

xsol = x;                                            % Final solution
fsol = f;                                            % Final function value

% Output structure
output.numIter = numIter;                            
output.funcCount = func_eval(1:numIter);             
output.fvalues = fvalues(1:numIter);                 
output.funcCounts = nfeval(1:numIter);               
output.truefunc = truef(1:numIter);                  
output.numsfe = num_success_iter_fe;                 
output.numusfe = num_unsuccess_iter_fe;              
output.numsle = num_success_iter_le;                 
output.numusle = num_unsuccess_iter_le;              
output.flag = exitflag;                              
end
