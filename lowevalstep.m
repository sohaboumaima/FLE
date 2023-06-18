function [x, f, alpha, success, func_eval] = lowevalstep(fun, x, f, alpha, pr, ni, li, ui, Ai, Z, W, dirgen, preslineq, presbounds, presleq, m, p0, nbunc, tol_feas)
% Inputs:
%   - fun: The objective function.
%   - x: The current point.
%   - f: The current function value at point x.
%   - alpha: The current step size.
%   - pr: Structure containing problem data.
%   - ni: Number of inequality constraints.
%   - li: Lower bounds for inequality constraints.
%   - ui: Upper bounds for inequality constraints.
%   - Ai: Matrix of inequality constraint coefficients.
%   - Z: Matrix of equality constraint coefficients.
%   - W: Basis for the null space of A (equality constraint matrix).
%   - dirgen: Method that is used for generating the directions with
%       0: Deterministic variant - always choose the generators of an 
%       approximate tangent cone (but randomly ordered)
%       1: Proceed with a random subset selected within the 
%       generators of an approximate tangent cone
%       2: Attempt to benefit from linear subspaces included in 
%       an approximate tangent cone to reduce the size of the (random) 
%       polling set, otherwise draw directions similarly to 1
%   - preslineq: Flag indicating the presence of linear equality constraints.
%   - presbounds: Flag indicating the presence of bounds.
%   - presleq: Flag indicating the presence of inequality constraints.
%   - m: Number of linear equality constraints.
%   - p0: Probability of selecting feasible directions.
%   - nbunc: Number of unconstrained directions.
%   - tol_feas: Tolerance for feasibility checks.

% Outputs:
%   - x: Updated point after the step.
%   - f: Updated function value at the new point.
%   - alpha: The updated step size.
%   - success: Flag indicating the success of the step.
%   - func_eval: Number of function evaluations performed.

   
    Aineq  = pr.Aineq;
    lbineq = pr.lbineq;
    ubineq = pr.ubineq;
    lb     = pr.lb;
    ub     = pr.ub;
    gamma      = 2;
    theta      = 0.5;
    func_eval  = 0;
    n = length(x);
    D = [];
    [Iap,Iam,nbgen] = activeconsAiZ(x,alpha,li,ui,Ai,Z);
%
    Iip = setdiff(1:ni,Iap);
    Iim = setdiff(1:ni,Iam);
    Ya  = [];
%
    if (nbgen==0)
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
                D = 2*rand(n,max(ceil(nbunc/2))) - 1;
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
			[Ys_R,Yc_R,~] = gentech(pr.Aeq,Ai,Z,Iap,Iam);
%             if (size(Ys_R,1) ~= size(Z, 2) || size(Yc_R,1) ~= size(Z, 2))
%                  if (size(Ys_R,1)>0 || size(Yc_R,1)>0)
%                     [Ys_R,Yc_R,calldd] = gentech(Aeq,Ai,Z,Iap,Iam);
%                  end
%             end
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
		warning('Iteration %d: Empty Polling set %d', dirgen);
	else
		cardD = size(D,2);
	end

    count_dir = 1;
    success = 0;

    while (~success) && (count_dir <= cardD)

	d = D(:,count_dir);
	xtemp = x + alpha * d;
        rho = min(1e-5, 1e-5*alpha^2);
%\
%	Test for the feasibility of the point

        feas = isfeasible(xtemp,presbounds, preslineq, pr.lb, pr.ub, pr.Aeq, pr.beq, pr.Aineq, pr.lbineq, pr.ubineq, tol_feas);

	if (~feas)
		count_dir = count_dir +1;
	else
		ftemp = feval(fun,xtemp, 0); 
		func_eval = func_eval + 1;
		if (ftemp < f - rho*(norm(d)^2)) % Test for sufficient decrease at xtemp
			success = 1;
		else
			count_dir = count_dir +1;
		end
	end
    end

    if success
        f = ftemp;
        x = xtemp;
        alpha = gamma * alpha;
    else
        alpha = theta * alpha;
    end
