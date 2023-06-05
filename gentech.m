function [Ys,Yc,calldd] = gentech(Ae,Ai,Z,Ap,Am)
%
%	Computation of the generators of the approximate tangent cone as
%	defined in Griffin et al (SISC, 2008).
%
%	[Ys,Yc] = lisgen(A,Z,Ap,Am) returns the generators of the 
%	polar of a cone spanned by a linearly independent set of the columns of 
%	[Ae' Ai(Ap,:) -Ai(Am,:)] that violate the constraints at (x,alpha).
%
%	Inputs:
%
%	Ae: equality constraint matrix
%
%	Ai: inequality constraint matrix
%
%	Z: orthonormal basis for the null space of A
%
%	Ap: upper inequalities violated for the current iterate/step 
%		size
%
%	Am: lower inequalities violated for the current iterate/step 
%		size
%
%	Outputs:
%
%	Ys: Linear generators for the polar of a cone defined by a
%		subset of the columns of [A' I -I]
%	Yc: Positive generators for the same polar
%
%	calldd: boolean indicating a call to the double description
%		method
%
%	C. W. Royer - March 26, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
n=size(Z,1);
En = eye(n);
calldd = 0;
%
r = rank(Ae);
if r<size(Ae,1)
%	Modifying Ae to obtain a full rank matrix
	Q = orth(Ae');
	Ae = Q';
end
%
Ib = intersect(Ap,Am);
Iu = setdiff(Ap,Ib);
Il = setdiff(Am,Ib);
%
W = Z*Z'*Ai';
Vp = [];
if ~isempty(Iu)
	Vp = [Vp W(:,Iu)];
end
if ~isempty(Il)
	Vp = [Vp -W(:,Il)];
end
Vl = [Ae'];
if ~isempty(Ib)
	Vl = [Vl W(:,Ib)];
	Zl = null(Vl');
else
	Zl = Z;
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine if a call to the double description is necessary
%
if isempty(Vp)
	Ys = Z'*Zl;
	Yc = [];
elseif isempty(Zl)
	Ys = [];
	Yc = [];
else
%
	Q = Zl'*Vp;
	rp=rank(Q);
	sp=size(Q);
	if rp==min(sp)
%		Non-degenerate case - Direct computation of the generators
		[Qi,Ri]=qr(Q);
		if (sp(1)>sp(2))
%             Ys = Qi(:,sp(2)+1:end);

			Ys = Zl*Qi(:,sp(2)+1:end);
            Ys = Z'*Ys;
		else
			Ys = [];
		end
		r2=min(sp(2),size(Qi,2));
		Qi=Qi(:,1:r2);
		Ri=Ri(1:r2,1:r2);
		warning('off');
		Yc=-Zl*(Qi/Ri');
        Yc = Z'*Yc;
% 		Yc=-(Qi/Ri');

		warning('on');
	else
%		Degenerate case - calling the double description method
		[Ys,Yc] = calldoubledescLI(Ae,Ai,Ap,Am,Z);
		calldd = 1;
	end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
