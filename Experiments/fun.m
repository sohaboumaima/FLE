function [fx] = fun(f, x, istrue, varargin)


% This function calculates the value of f at x and add noise
% based on the provided inputs. The noise can be either deterministic or 
% stochastic; additive or multiplicative.
% Inputs:
% -------
% f: A function handle that represents the function f.
% x: The input value at which to evaluate the function f.
% istrue: A flag that determines whether the returned value is the true
%           function or the noisy version.
% Outputs:
% --------
% fx: The calculated value of f at x.


probtype = '';
fx = f(x);

if (nargin<3)
    istrue = 0;
end

if (~istrue)
    switch probtype
        % Add stochastic noise
        case 'noisy'
            sigma=10^-3;
            u = sigma*(-1+2*rand(1,1));
            fx= fx*(1+u);
            % fx = fx + u;

        % Add deterministic noise
        case 'wild'
            sigma=10^-3;
            phi = 0.9*sin(100*norm(x,1))*cos(100*norm(x,inf)) + 0.1*cos(norm(x,2));
            phi = phi*(4*phi^2 -3);
            fx = (1 + sigma*phi)*fx;
            % fx = sigma*phi + fx;
    end
end
end
