function [xsol,fsol,exitflag,output] = dspfd(fun,x0,alpha0,dirgen,lb,ub,Aeq,...
beq,Aineq,lbineq,ubineq,seed,maxevals,ftarget,tolf,tol_feas,varargin)
%
%	A wrapper for calling the dspfd (Direct Search based on Probabilistic 
%	Feasible Descent) method.
%
%	xsol = dspfd(fun,x0) applies the DSPFD algorithm to minimize the function 
%	fun starting at x0; and returns the best solution point obtained xsol.
%
%	xsol = dspfd(fun,x0,alpha0) applies the DSPFD algorithm to minimize the 
%	function fun starting at x0 with the initial step size alpha0 (set
%	alpha0=[] for the default value).
%
%	xsol = dspfd(fun,x0,alpha0,dirgen) applies the strategy dirgen to generate 
%	the polling directions (set dirgen=[] for the default value).
%
%	xsol = dspfd(fun,x0,alpha0,dirgen,lb,ub) applies the DSPFD algorithm to 
%	minimize the function fun subject to the bounds lb and ub. It is 
%	assumed that x0 satisfies the bound constraints (set lb=[] and ub=[] if 
%	no bounds are enforced).
%
%	xsol = dspfd(fun,x0,alpha0,dirgen,lb,ub,Aeq,beq) applies the DSPFD 
%	algorithm to minimize fun subject to the bound constraints lb <= x <= ub 
%	as well as the linear equality constraints Aeq*x=beq (set Aeq=[] and beq=[]
%	in absence of such constraints).
%
%	xsol = dspfd(fun,x0,alpha0,dirgen,lb,ub,Aeq,beq,Aineq,lbineq,ubineq) 
%	applies the DSPFD algorithm to minimize fun subject to the bound 
%	constraints lb <= x <= ub as well as the linear equality constraints 
%	Aeq*x=beq and the linear inequalty constraints lbineq <= Aineq*x <= ubineq 
%	(set Aeq=[] and beq=[] in absence of such constraints).
%
%	xsol = dspfd(fun,x0,alpha0,dirgen,lb,ub,Aeq,beq,Aineq,lbineq,ubineq,seed) 
%	applies the DSPFD method with a specified random seed (set seed=[] to use
%	the default setting).
%
%	xsol = dspfd(fun,x0,alpha0,dirgen,lb,ub,Aeq,beq,Aineq,lbineq,ubineq,seed,...
%	maxevals) applies the DSPFD method using at most maxevals calls to the 
%	objective function (set maxevals=[] for the default setting).
%
%	xsol = dspfd(fun,x0,alpha0,dirgen,lb,ub,Aeq,beq,Aineq,lbineq,ubineq,seed,...
%	maxevals,ftarget,tolf) applies the DSPFD method and uses the stopping 
%	criterion f-ftarget < tolf*(f(x0)-ftarget) (set ftarget=[] and/or tolf=[] 
%	to use the default values).
%
%	xsol = dspfd(fun,x0,alpha0,dirgen,lb,ub,Aeq,beq,Aineq,lbineq,ubineq,seed,...
%	maxevals,ftarget,tolf,tol_feas) applies the DSPFD method and tolerates a 
%	violation of the constraints within a tolerance of tol_feas (set 
%	tol_feas=[] to use the default value).
%
%	[xsol,fsol] = dspfd(fun,x0,...) returns the best solution computed xsol 
%	and the corresponding value fsol=fun(xsol).
%
%	[xsol,fsol,exitflag] = dspfd(fun,x0,...) additionally returns a flag 
%	related to the run of the method.
%
%	[xsol,fsol,exitflag,output] = dspfd(fun,x0,...) additionally returns a 
%	structure containing detailed information about the algorithmic run.
%
%	-----------------
%	Input parameters
%	-----------------
%
%	fun: Objective function handle
%
%	x0: Initial feasible point - Must be provided
%
%	alpha0: Initial step size - Can be empty, 1 by default
%
%   dirgen: Method that is used for generating the directions
%		Default value if empty: 2	
%       0: Deterministic variant - always choose the generators of an 
%		approximate tangent cone (but randomly ordered)
%       1: Proceed with a random subset selected within the 
%       generators of an approximate tangent cone
%       2: Attempt to benefit from linear subspaces included in 
%       an approximate tangent cone to reduce the size of the (random) 
%		polling set, otherwise draw directions similarly to 1
%
%	lb/ub: Bound vectors (infinite if not assigned)
%		Lower and upper bound vectors. Some components might be
%		unbounded, in which case the corresponding bound is equal to
%		+/- 1e20 or +/- Inf
%
%   Aeq: matrix (should be of full row rank) associated to the m<n
%   equality constraints - Empty if not assigned
%
%   beq: vector of size m associated to the m<n equality constraints
%		Empty if not assigned
%
%   Aineq: matrix associated to linear inequality constraints
%		Empty if not assigned
%
%   lbineq,ubineq: lower and upper bound vectors corresponding to the
%   linear inequality constraints
%		Empty if not assigned
%
%   seed: Integer used to initialize the random number generator
%		Default value if not assigned: 1
%
%   maxevals: Maximum number of function evaluations allowed
%		Default value if not assigned: 2000*length(x0)
%
%   ftarget: Target value for the method
%		Default value if not assigned: -Inf
%
%   tolf: Tolerance with respect to the target value
%		Default value if not assigned: 0
%
%   tol_feas: Feasibility tolerance to consider a point as feasible
%		Default value if not assigned: 10^(-3)
%		(as in MATLAB's patternsearch function)
%
%	-----------------
%	Output parameters
%	-----------------
%
%   xsol,fsol : Best point xsol and its value fsol found by the method
%
%   exitflag: flag indicating what happened in the method
%       0: Method stopped because the target function value was
%       reached
%       1: Method stopped because the step size was too small
%       (expected behaviour)
%       2: Method stopped because the budget of function evaluations
%       was exceeded
%       3: No feasible point was found past the initial one
%       -1: An error occurred and the method did not complete its
%       task
%
%   output:  structure containing the following information
%       funcCount: number of function evaluations performed in total
%       numIter: number of iterations
%       meanDirgen: number of directions generated per iteration
%       meanFeaspts: mean number of feasible points per iteration
%       fracDD: fraction of iterations calling the description method
%       histF: history matrix of size (numIter+1)*2
%           histF(:,1): function values at the current iterate
%           histF(:,2): number of function evaluations
%
%   C. W. Royer - March 26, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin<2)
		error('dspfd: At least two inputs required');
end
%
if isempty(x0)
	error('dspfd: Empty initial point, please input a feasible vector');
else
	x0 = x0(:);
end
%
if (nargin >= 3)
	if isempty(alpha0)
		alpha0 = 1;
	elseif (alpha0<=0)
		warning('dspfd: Nonpositive step size, replaced by 1');
		alpha0 = 1;
	end
else
	alpha0 = 1;
end
%
if (nargin < 4) || isempty(dirgen)
	dirgen = 2;
end
%
if (nargin >= 5)
	if isempty(lb)
		lb = -Inf*ones(length(x0),1);
	elseif length(lb)~=length(x0)
		error('dspfd: Incorrect lower bound vector');
	end
else
	lb = -Inf*ones(length(x0),1);
end
%
if (nargin >= 6)
	if isempty(ub)
		ub = Inf*ones(length(x0),1);
	elseif length(ub)~=length(x0)
		error('dspfd: Incorrect upper bound vector');
	end
else
	ub = Inf*ones(length(x0),1);
end
%
if sum(lb<=ub)~=length(x0)
	error('dspfd: Infeasible bounds');
end
%
if (nargin == 7)
	error('dspfd: Argument beq missing');
end
%
if (nargin >= 8)
	if ((isempty(Aeq) && ~isempty(beq)) || (~isempty(Aeq) && isempty(beq)))
		error('dspfd: Aeq or beq is empty when it should not be');
	elseif (~isempty(Aeq) && ~isempty(beq))
		if size(Aeq,1)~=length(beq)
			error('dspfd: Dimensions of Aeq and beq do not match');
		end
		if size(Aeq,2)~=length(x0)
			error('dspfd: Dimensions of Aeq and x0 do not match');
		end
	end
else
	Aeq = [];
	beq = [];
end
%
if (nargin==9) || (nargin==10)
	error('dspfd: Arguments lbineq and/or ubineq missing');
end
%
if (nargin >= 11)
	if ~isempty(Aineq) && ~isempty(lbineq) && ~isempty(ubineq)
		if size(Aineq,2)~=length(x0)
			error('dspfd: Dimensions of Aineq and x0 do not match');
		end
%
		if (size(Aineq,1)~=length(lbineq) || size(Aineq,1)~=length(ubineq))
			error('dspfd: Dimensions of Aineq and lbineq/ubineq do not match');
		end
%
		if (lbineq > ubineq)
			error('dspfd: Infeasible inequality constraints');
		end
%
	elseif ~(isempty(Aineq) && isempty(lbineq) && isempty(ubineq))
		if isempty(Aineq)
			error('dspfd: Matrix Aineq empty while it should not be');
		end
		if isempty(lbineq) && isempty(ubineq)
			error('dspfd: lbineq and ubineq empty but should not be');
		end
		if isempty(lbineq)
			lbineq = -Inf*ones(size(Aineq,1),1);
			if size(Aineq,1)~=length(ubineq)
				error('dspfd: Dimensions of Aineq and ubineq do not match');
			end
		end
		if isempty(ubineq)
			ubineq = Inf*ones(size(Aineq,1),1);
			if size(Aineq,1)~=length(lbineq)
				error('dspfd: Dimensions of Aineq and lbineq do not match');
			end
		end
%
	end
%
else
	Aineq = [];
	lbineq = [];
	ubineq = [];
end
%
if (nargin < 12) || isempty(seed)
	seed=1;
end
%
if (nargin < 13) || isempty(maxevals)
	maxevals = 2000*length(x0);
end
%
if (nargin < 14) || isempty(ftarget)
	ftarget = -Inf;
% 	ftarget = 1e-9;
end
%
if (nargin < 15) || isempty(tolf)
	tolf = 0;
elseif (tolf<0)
	error('dspfd: tolf parameter must be nonnegative');
end
%
if (nargin < 16) || isempty(tol_feas)
	tol_feas = 1e-8;	
elseif (tol_feas<0)
	error('dspfd: tol_feas parameter must be nonnegative');
end
%
% Calling the dspfd method
[xsol,fsol,exitflag,output] = dspfdalgo(fun,x0,alpha0,dirgen,lb,ub,Aeq,beq,...
Aineq,lbineq,ubineq,seed,maxevals,ftarget,tolf,tol_feas);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
