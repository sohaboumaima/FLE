function [xsol,fsol,exitflag,output]= dspfdalgo(fun,x0,alpha0,dirgen,lb,ub,...
Aeq,beq,Aineq,lbineq,ubineq,seed,maxevals,ftarget,tolf,tol_feas)
%
%	Direct Search based on Probabilistic Feasible Descent.
%
%	Implementation of direct-search methods based upon active-set 
%	considerations, with the aim of solving minimization problems with 
%	variables subject to bound and linear equality constraints.
%
%	[xsol,fsol,exitflag,output] = dspfd(fun,x0,alpha0,dirgen,lb,ub,...
%	Aeq,beq,Aineq,lbineq,ubineq,seed,maxevals,ftarget,tolf,tol_feas) 
%	attempts to solve the optimization problem
% 		minimize f(x) subject to 	Aeq*x=beq
%									lbi <= Aineq*x <= ubi 
%									lb <= x <= ub.
%
%
%	Inputs:
%
%	fun: The objective function
%
%	x0: Initial (feasible) point
%
%	alpha0: Initial step size
%
%	dirgen: Method that is used for generating the directions with
%       0: Deterministic variant - always choose the generators of an 
%       approximate tangent cone (but randomly ordered)
%       1: Proceed with a random subset selected within the 
%       generators of an approximate tangent cone
%       2: Attempt to benefit from linear subspaces included in 
%       an approximate tangent cone to reduce the size of the (random) 
%       polling set, otherwise draw directions similarly to 1
%
%	lb,ub: Lower and upper bound vectors. Some components might be
%	unbounded, in which case the corresponding bound is equal to
%	+/- 1e20 or +/- Inf
%
%	Aeq: matrix (should be of full row rank) associated to the m<n
%	equality constraints
%
%	beq: vector of size m associated to the m<n equality constraints
%
%	Aineq: matrix associated to linear inequality constraints
%
%	lbineq,ubineq: lower and upper bound vectors corresponding to the
%	linear inequality constraints
%
%	seed: Integer used to initialize the random number generator
%
%	maxevals: Maximum number of function evaluations allowed
%
%	ftarget: Target value for the method
%
%	tolf: Tolerance with respect to the target value
%
%	tol_feas: Feasibility tolerance to consider a point as feasible
%
%	Outputs:
%
%	xsol,fsol : Best point xsol and its value fsol found by the method
%
%	exitflag: flag indicating what happened in the method
%		0: Method stopped because the target function value was
%		reached
%		1: Method stopped because the step size was too small
%		(expected behaviour)
%		2: Method stopped because the budget of function evaluations
%		was exceeded
%		3: No feasible point was found past the initial one
%		-1: An error occurred and the method did not complete its
%		task	
%
%	output:  structure containing the following information
%		funcCount: number of function evaluations performed in total
%		numIter: number of iterations
%		meanDirgen: number of directions generated per iteration
%		meanFeaspts: mean number of feasible points per iteration
%		fracDD: fraction of iterations calling the description method
%		histF: history matrix of size (numIter+1)*2
%			histF(:,1): function values at the current iterate
%			histF(:,2): number of function evaluations
%
%	C. W. Royer - March 26, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	Initialization and setting of auxiliary parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format compact;
format long;
n = length(x0);
% rand('state',seed);
In = eye(n);
exitflag = -1;
tol_alpha = alpha0*(1e-10);
%tol_feas = 1e-8;
%tol_feas = 1e-3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Step size update parameters and corresponding probability of descent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dirgen > 0
	gamma = 2;
else
	gamma = 1;
%	Other possible value
%	gamma = 2;
end
theta = 0.5; 
p0 = log(theta)/log(theta/gamma);
% Minimum number of directions required in the unconstrained case
nbunc = ceil(log(1-log(theta)/log(gamma))/log(2));
% nbunc = 2*n;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	Counters and indexes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
funcCount = 0;
numIter	= 0;
numFeaspts = 0;
numDirgen = 0;
sumdd = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Before the loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
x = x0;
f = feval(fun,x, 0);
diff0 = tolf*(f-ftarget);
funcCount = funcCount+1;
histF = [f funcCount];
tf = [fun(x, 1)];
alpha = alpha0;
alphamax = 20*alpha0;
%
stopcrit = 0;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%% Determination of the presence of bounds %%%%%%%%%%%%%%
%
lowbnd = find(lb > -1e20);
upbnd = find(ub < 1e20);
presbounds = ~isempty(upbnd) || ~isempty(lowbnd);
%
%%%%%%%%%%%%%%%%%%%%%%%% Checking for equality constraints %%%%%%%%%%%%%%%%%%%%
%
if ~isempty(Aeq)
	presleq = 1;
	for i=1:size(Aeq,1)
		auxa = norm(Aeq(i,:));
		Aeq(i,:) = Aeq(i,:)/auxa;
		beq(i) = beq(i)/auxa;
	end
	Z = null(Aeq);
	m = size(Z,2);
%
        if m==0
%               We are already at the solution, or the problem is infeasible
                stopcrit = 1;
        end
else
	presleq = 0;
	Z = In;
        m = n;
end
%
%%%%%%%%%%%%%%%%%%%% Gathering linear inequalities and bounds %%%%%%%%%%%%%%%%%
%
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
%
% The projected normal vectors
W = Z*Z'*Ai';
%
% 
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% fprintf( '  it      neval        fvalue       success      alpha  \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Main loop 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while ~stopcrit
%
	numIter = numIter + 1;
%
%	Computation of the forcing function at the current iterate
%
%	rho = 1e-4*(alpha^2);
%	Alternative
	rho = min(1e-4,1e-4*(alpha^2));
%
%	Computation of the alpha-active and inactive constraints
	[Iap,Iam,nbgen] = activeconsAiZ(x,alpha,li,ui,Ai,Z);
%
	Iip = setdiff(1:ni,Iap);
	Iim = setdiff(1:ni,Iam);
%
	Ya = [];
%
%%%%%%%%%%%%%%%%%%% Computation of the direction sets %%%%%%%%%%%%%%%%%%%%%%%%%
%
	if (nbgen==0)
%
%		No active bounds - situation similar to linear 
%		equality-constrained (or unconstrained) case
%	
		switch dirgen
			case 0
%				Positive spanning set corresponding to the null space of the 
%				linear equality constraints (if any)
				if ~preslineq && ~presbounds
					D = [Z -Z];
				else
					D = [W -W];
                end
                
%
            case 1
                
%				Random sample of the previous one
				if ~preslineq && ~presbounds
					nbdir = ceil(2*m*p0);
					D = randcolZsearch(Z,1:m,1:m,nbdir);
				else
					nbdir = ceil(2*ni*p0);
					D = randcolZsearch(W,1:ni,1:ni,nbdir);
				end
%
			case 2
%               Using uniform distribution in the null space
%				of the equality constraints
%
				D = 2*rand(n,max(ceil(nbunc/2))) - 1;    %(1,ceil(nbunc/2)))-1;
				if ~preslineq && ~presbounds
					D = Z'*D;
				else
					D = W'*D;
				end
%
				for i=1:size(D,2)
					D(:,i) = D(:,i)/norm(D(:,i));
				end
%
				if ~preslineq && ~presbounds
					D = Z*D;
				else
					D = W*D;
				end
%
				D = [D -D];
%
                otherwise
                        error('dirgen value not supported');
                end
%
	else
%		Case nbgen >=1 :
%
		if (preslineq || (presleq && presbounds))
%
			[Ys_R,Yc_R,calldd] = gentech(Aeq,Ai,Z,Iap,Iam);
            if (size(Ys_R,1) ~= size(Z, 2) || size(Yc_R,1) ~= size(Z, 2))
                 if (size(Ys_R,1)>0 || size(Yc_R,1)>0)
                    [Ys_R,Yc_R,calldd] = gentech(Aeq,Ai,Z,Iap,Iam);
                 end
            end
            if (numIter == 38)
                [Ys_R,Yc_R,calldd] = gentech(Aeq,Ai,Z,Iap,Iam);
            end
            if isempty(Ys_R)
                Ys = [];
            else
                Ys = Z*Ys_R;
            end
            if isempty(Yc_R)
                Yc = [];
            else
                Yc = Z*Yc_R;
            end			
            sumdd = sumdd + calldd;
%
%			Orthonormalization process
%
			if ~isempty(Ys)
				for i=1:size(Ys,2)
					wsi = Ys(:,i);
					Ys(:,i) = wsi/norm(wsi);
				end
%
				if ~isempty(Yc)
%
					Yc = Yc - Ys*Ys'*Yc;
					for i=1:size(Yc,2)
						wci = Yc(:,i);
						Yc(:,i) = wci/norm(wci);
					end
				end
			end
%
%			Additional directions - currently empty			
%				Possibilities: coordinate vectors, more generators, etc.
%
			Ya = [];
%
%
		else
%			Simply select the complement of the active coordinate
%			vectors
			Is = intersect(Iip,Iim);
			Icp = setdiff(Iip,Is);
			Icm = setdiff(Iim,Is);
			Ys = W(:,Is);
			Yc = [W(:,Icp) -W(:,Icm)];
		end
%
		dyy = 2*size(Ys,2)+size(Yc,2);

%		
		if dyy>0
%
			switch dirgen
				case 0
					D = [Ys -Ys Yc];
					if ~isempty(D)
						pDD = randperm(size(D,2));
						D = D(:,pDD);
					end
%						
				case 1
%					Ensuring a p>p0_probability of using a 
%					feasible descent direction among the 
%					columns of Z and their opposite
					nbdir = ceil(dyy*p0);
					D = randcolZsearch([Ys -Ys Yc],...
					1:dyy,[],nbdir);
%
                case 2
%                   Using uniform distribution in subspaces
%					if possible, in cones otherwise
% 
					D = [];
%
%					Subspace vectors
					nbdir = ceil(nbgen*p0);
%
					if ~isempty(Ys)
						Dsb = 2*rand(size(Ys,2),...
						max(1,ceil(nbunc/2)))-1;
						Dsb = Ys*Dsb;
						for i=1:size(Dsb,2)
							Dsb(:,i) = ...
							Dsb(:,i)/norm(...
							Dsb(:,i));
						end
						D = [Dsb -Dsb];	
					end
%					Cone vectors
					if ~isempty(Yc)
						nyc = size(Yc,2);
						nbdir = ceil(nyc*p0);
%
						Dco = randcolZsearch(Yc,...
						1:nyc,[],nbdir);
						D = [D Dco];
					end
%
			otherwise
				error('dirgen value not supported');
			end
%
		else
			D = [];
		end
%
%		Adding additional directions if any
		if ~isempty(Ya)
			if dirgen>0
				ca = size(Ya,2);
				dda = ceil(ca*p0);
				Ya = Ya(:,dda);
			end
			D = [D Ya];
		end
	end
%
%	Additional directions II - Infeasible constraints normals
%	Set the boolean addII to 1 to activate this option
	addII = 1;
%
	if (nbgen > 0) && addII
%
%
		Ias = intersect(Iap,Iam);
		Iacp = setdiff(Iap,Ias);
		Iacm = setdiff(Iam,Ias);
%
%		Use same ideas than for "tangent" directions
		switch dirgen
			case 0
				Dn = reducedstep(alpha,x,lb,ub,lbineq,...
				ubineq,Aineq,[W(:,Iap) -W(:,Iam)],tol_feas);
				if ~isempty(Dn)
					pDD = randperm(size(Dn,2));
					Dn = Dn(:,pDD);
					D = [D Dn];
				end
%
			case 1
				if ~isempty(Iap)
					Dnp = reducedstep(alpha,x,...
					lb,ub,lbineq,ubineq,Aineq,...
					W(:,Iap),tol_feas);
					if ~isempty(Dnp)
						cardDn  = size(Dnp,2);
						nbnrml = min(1,...
						ceil(cardDn*p0/2));
						tpm = randperm(cardDn);
						D = [D Dnp(:,...
						tpm(1:nbnrml))];
					end
				end
%
				if ~isempty(Iam)
					Dnm = reducedstep(alpha,x,...
					lb,ub,lbineq,ubineq,Aineq,...
					-W(:,Iam),tol_feas);
					if ~isempty(Dnm)
						cardDn  = size(Dnm,2);
						nbnrml = min(1,...
						ceil(cardDn*p0/2));
						tpm = randperm(cardDn);
						D = [D Dnm(:,...
						tpm(1:nbnrml))];
					end
				end
%
			case 2
				nsbd = length(Ias);
				if nsbd~=0
					Was = W(:,Ias);
					dnn = 2*rand(nsbd,1)-1;
					dnn = Was*dnn;
					Dnsp = reducedstep(alpha,x,...
					lb,ub,lbineq,ubineq,Aineq,...
					[dnn -dnn],tol_feas);
					D = [D Dnsp];
				end
%
				nncone = length(Iacp)+length(Iacm);
				if nncone~=0
					nbnormaldir = max(1,...
					ceil(nncone*p0));
%
					Wac = [W(:,Iacp) -W(:,Iacm)];
					t = randperm(nncone);
					Dncoco = reducedstep(alpha,x,...
					lb,ub,lbineq,ubineq,Aineq,...
					Wac(:,t(1:nbnormaldir)),...
					tol_feas);
					D = [D Dncoco];
				end
%
		otherwise
			error('dirgen value not supported');
		end	
	end			
%
	if isempty(D)
		cardD = 0;
		warning('Iteration %d: Empty Polling set %d',numIter,...
		dirgen);
	else
		cardD = size(D,2);
	end
	numDirgen = numDirgen+cardD;

    

%
%%%%%%%%%%%%%%%%%%% Loop on the directions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
	count_dir = 1;
	success = 0;
%
	oldnumFeaspts = numFeaspts;
%
	while (~success) && (count_dir <= cardD) && (funcCount < maxevals)
%
		d = D(:,count_dir);
		xtemp = x + alpha*d;
%
%		Test for the feasibility of the point
%
		if presbounds
			xl = (min(xtemp-lb) > -tol_feas);
			xu = (max(xtemp-ub) < tol_feas);
		else
			xl = 1;
			xu = 1;
		end
%
		if ~isempty(Aeq)
			Ax = (norm(Aeq*xtemp-beq,Inf) < tol_feas);
		else
			Ax = 1;
		end
%
		if preslineq
			xli = (min(Aineq*xtemp-lbineq) > -tol_feas);
			xui = (max(Aineq*xtemp-ubineq) < tol_feas);
			Aix = (xli && xui);
		else
			Aix = 1;
		end
% %
		if (~(xl && xu && Ax && Aix))
% 			fprintf('Unfeasible tentative point %d\n',dirgen);
			count_dir = count_dir +1;
		else
%
			numFeaspts = numFeaspts +1;
%
%			Test for sufficient decrease at xtemp
%
			ftemp = feval(fun,xtemp, 0);
			funcCount = funcCount + 1;
			if (ftemp < f - rho*(norm(d)^2))
				success = 1;
			else
				count_dir = count_dir +1;
			end
%
		end
	end
%
%%%%%%%%%%%%%%%%%%% End of the iteration  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	Have we sufficiently improved the current function value ?
%
	if success
		f = ftemp;
		x = xtemp;
		alpha = min(gamma*alpha,alphamax);
	else
		alpha = theta*alpha;
	end
%
	histF = [histF;f funcCount];
    tf = [tf, fun(x, 1)];
%
%	Is the stopping criterion satisfied ?
%
	stopcrit = ((f-ftarget < diff0) ||  (alpha < tol_alpha) || ...
	(funcCount >= maxevals));

    

%     fprintf( '%5d %+.4e %+.4e  %5d %+.4e\n', numIter, funcCount, f, success, alpha );


%
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	END OF THE MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Final adjustments
%
xsol = x;
fsol = f;
output.numIter = numIter;
output.funcCount = funcCount;
output.meanFeaspts = round(numFeaspts/numIter);
output.meanDirgen = round(numDirgen/numIter);
output.fracDD = sumdd/numIter;
output.histF = histF;
output.truefunc = tf;
% 
%
if (alpha >=tol_alpha)
	if (funcCount==1)
%		No initial feasible point was found
		exitflag = 3;
	else 
%		Budget of function evaluations exceeded
		exitflag = 2;
	end
else
	if ( f-ftarget < diff0)
		exitflag = 0;
	else
%		Step size shrunk below tolerance
		exitflag = 1;
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
