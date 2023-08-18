function [Iap,Iam,nbgen] = activeconsAiZ(x,alpha,li,ui,Ai,Z)
%
%	Computation of the index of active linear constraints at a given point.
%
%	[Iap,Iam,nbgen] = activeboundsZ(x,alpha,li,ui,Ai,Z) returns the
%	indexes of columns of [Ai -Ai] for which a displacement from a step 
%	length ALPHA is sufficiently feasible.
%
%	Inputs:
%
%	x: current point
%
%	alpha: current step length
%
%	li: vector of lower bounds (LI(i) = -1e20 if unbounded)
%
%	ui: vector of upper bounds (UI(i) = +1e20 if unbounded)
%
%	Ai: matrix of linear inequality constraints and/or bounds
%
%	Z: orthonormal basis for a subspace of R^n, with n=dim(X)
%
%	Outputs:
%
%	Iap: indexes of columns of AI that cannot be used
%	
%	Iam: Indexes of columns of -AI that cannot be used
%
%	nbgen: Number of generators
%
%	C. W. Royer - March 26, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
n = length(x);
W = Z*Z'*Ai';
nz = sqrt(diag(Ai*Z*Z'*Ai'));
iz = find(nz < 1e-15);
%
ni = size(Ai,1);
Eni = eye(ni);
%
%
% Feasibility tolerance for defining the approximately active constraints
%tolfeas=0;%No tolerance
tolfeas=1e-15;
%tolfeas=1e-3;% A mild tolerance - initial choice of MATLAB patternsearch
%
for i=1:ni
	Vp(i) = (-(Ai(i,:)*x-ui(i)-tolfeas))/nz(i);
	Vm(i) = ((Ai(i,:)*x-li(i)+tolfeas))/nz(i);
end
%
if ~isempty(iz)
	for i=1:length(iz)
		j = iz(i);
		if (abs(Ai(j,:)*x-ui(j))<tolfeas)
			Vp(j) = 0;
		else
			Vp(j) = Inf;
		end
%
		if (abs(Ai(j,:)*x-li(j))<tolfeas)
			Vm(j) = 0;
		else
			Vm(j) = Inf;
		end
	end
end
%
% alpha = min(alpha,1e-3);
%	
Iap = find(Vp <= alpha);
Iam = find(Vm <= alpha);
nbgen = length(Iap)+length(Iam);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
