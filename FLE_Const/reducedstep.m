function D = reducedstep(alpha,x,lb,ub,lbineq,ubineq,Aineq,V,tolfeas)
%
%	Computation of a reduced step in infeasible directions.
%
%	D = reducedstep(alpha,x,lb,ub,lbineq,ubineq,Aineq,V,tolfeas) returns 
%	(if possible) directions corresponding to those in V such that for 
%	every d in D, one has lb <= x+alpha*d <= ub and
%	lbineq <= Aineq*(x+alpha*d) <= ubineq, up to a tolerance tolfeas.
%
% 	Inputs:
%
%	alpha: Current step size
%
%	x: Current point
%
%	lb: Vector of lower bounds
%
%	ub: Vector of upper bounds
%
%	lbineq: Vector of lower bounds
%
%	ubineq: Vector of upper bounds
%
%	Aineq: Matrix of inequality constraints
%
%	V: Matrix with columns consisting of unitary vectors
%
%	tolfeas: Feasibility tolerance for constraint violation
%
%	Outputs:
%
%	D: Directions resulting from scaling of those in V
%
%	C. W. Royer - March 26, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
sigmatol = 1e-3; % Minimum fraction of step size allowed
D = [];
v = zeros(length(x),1);
tolfeas = max(1e-15,tolfeas);
tolb = tolfeas;
%
for i=1:size(V,2)
	v = V(:,i);
%
%
	kl = find(x+alpha*v-lb <= -tolb);
	ku = find(x+alpha*v-ub >= tolb);
	if ~isempty(kl)
		zl = min((x(kl)-lb(kl)+tolb)./(-alpha*v(kl)));
	else
		zl = Inf;
	end
%
	if ~isempty(ku)	
		zu = min((ub(ku)-x(ku)+tolb)./(alpha*v(ku)));
	else
		zu = Inf;
	end
%
%
	if ~isempty(Aineq)
		vai = Aineq*v;
		jl = find(Aineq*x+alpha*vai-lbineq <= -tolfeas);
		ju = find(Aineq*x+alpha*vai-ubineq >= tolfeas);
		if ~isempty(jl)
			al = min((Aineq(jl,:)*x-lbineq(jl)+tolfeas)./...
			(-alpha*vai(jl)));
		else
			al = Inf;
		end
		if ~isempty(ju)	
			au = min((ubineq(ju)-Aineq(ju,:)*x+tolfeas)./...
			(alpha*vai(ju)));
		else
			au = Inf;
		end
	else
		al = Inf;
		au = Inf;
	end
%
	atemp = min([al au zl zu]);
	if ~isinf(atemp) && (atemp >= sigmatol)
		D = [D atemp*v];
	end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
