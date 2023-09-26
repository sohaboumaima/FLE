function [xsol, fsol, output] = FullLowEval(fun, x0, alpha0, pr, dirgen, maxevals, tol_feas, varargin)

% Check that the entries have the right formats and values


if (nargin<2)
		error('fle: At least two inputs required');
end
%
if isempty(x0)
	error('fle: Empty initial point, please input a feasible vector');
else
	x0 = x0(:);
end
%
if (nargin >= 3)
	if (alpha0<=0)
		warning('fle: Nonpositive step size, replaced by 1');
		alpha0 = 1;
	end
else
	alpha0 = 1;
end
%
if (nargin < 5)
	dirgen = 2;
end

unconstrained = 0;

if (nargin < 4)
	pr.lb = [];
    pr.ub = [];
    pr.Aeq    = [];
    pr.beq    = [];
    pr.Aineq  = [];
    pr.lbineq = [];
    pr.ubineq = [];
    unconstrained = 1;
end
%

% Fields in the struct

fields_to_check = {'lb', 'ub', 'Aeq', 'beq', 'Aineq','lbineq', 'ubineq'};

% Loop through the fields and check their existence
for i = 1:numel(fields_to_check)
    field_name = fields_to_check{i};
    
    if ~isfield(pr, field_name)
        pr.(field_name) = [];
    end
end

if isempty(pr.ub)
    pr.ub = Inf*ones(length(x0),1);
elseif length(pr.ub)~=length(x0)
    error('fle: Incorrect upper bound vector');
end


if isempty(pr.lb)
    pr.lb = -Inf*ones(length(x0),1);
elseif length(pr.lb)~=length(x0)
    error('fle: Incorrect lower bound vector');
end
%
if sum(pr.lb<=pr.ub)~=length(x0)
	error('fle: Infeasible bounds');
end
%
%



if ((isempty(pr.Aeq) && ~isempty(pr.beq)) || (~isempty(pr.Aeq) && isempty(pr.beq)))
	error('fle: Aeq or beq is empty when it should not be');
elseif (~isempty(pr.Aeq) && ~isempty(pr.beq))
	if size(Aeq,1)~=length(beq)
		error('fle: Dimensions of Aeq and beq do not match');
	end
	if size(pr.Aeq,2)~=length(x0)
		error('fle: Dimensions of Aeq and x0 do not match');
	end
end
%

if ~isempty(pr.Aineq) && ~isempty(pr.lbineq) && ~isempty(pr.ubineq)
	if size(pr.Aineq,2)~=length(x0)
		error('fle: Dimensions of Aineq and x0 do not match');
	end
%
	if (size(pr.Aineq,1)~=length(pr.lbineq) || size(Aineq,1)~=length(pr.ubineq))
		error('fle: Dimensions of Aineq and lbineq/ubineq do not match');
	end
%
	if (pr.lbineq > pr.ubineq)
		error('fle: Infeasible inequality constraints');
	end
%
elseif ~(isempty(pr.Aineq) && isempty(pr.lbineq) && isempty(pr.ubineq))
	if isempty(pr.Aineq)
		error('fle: Matrix Aineq empty while it should not be');
	end
	if isempty(pr.lbineq) && isempty(pr.ubineq)
		error('fle: lbineq and ubineq empty but should not be');
	end
	if isempty(pr.lbineq)
		pr.lbineq = -Inf*ones(size(pr.Aineq,1),1);
		if size(Aineq,1)~=length(pr.ubineq)
			error('fle: Dimensions of Aineq and ubineq do not match');
		end
	end
	if isempty(pr.ubineq)
		pr.ubineq = Inf*ones(size(pr.Aineq,1),1);
		if size(pr.Aineq,1)~=length(pr.lbineq)
			error('fle: Dimensions of Aineq and lbineq do not match');
		end
	end
%
end

%
if (nargin < 6)
	maxevals = 2000*length(x0);
end

if (nargin < 7) || isempty(tol_feas)
	tol_feas = 1e-8;	
elseif (tol_feas<0)
	error('fle: tol_feas parameter must be nonnegative');
end

[xsol, fsol, output] = FullLowEvalalg(fun, x0, alpha0, dirgen, pr, maxevals, tol_feas, unconstrained);

end

