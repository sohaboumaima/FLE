function [Ys,Yc] = calldoubledescLI(Ae,Ai,Ipos,Ineg,Z)

%	Computation of the generators of a cone obtained as the intersection
%	of a subspace and (possibly) some half-spaces.%%
%
%	[Ys,Yc] = calldoubledescLI(Ae,Ai,Ipos,Ineg,Z) returns a positive spanning 
%	set for the cone 
%	{x | Ae*x=0, [Ai]j*x>=0 for all j in Ipos,[Ai]_j*x<=0 for all j in Ineg}.
%
%	Inputs:
%
%	Ae: matrix of full row rank m < n defining linear equalities
%
%	Ai: matrix defining linear equalities
%
%	Ipos: index of components of the elements of the cone that must
%	be nonnegative
%
%	Ineg: index of components of the elements of the cone that must
%	be nonpositive
%
%	Z: orthonormal basis for the null space of Ae
%
%	Outputs: 
%
%	Ys: set of linear generators for the cone
%
%	Yc: set of positive generators for the cone
%
%	C. W. Royer - March 26, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
n= size(Ai,2);
r = rank(Ae);
if r<size(Ae,1)
%       Turning Ae into a full rank matrix
        Q = orth(Ae');
        Ae = Q';
end
%
W = Z*Z'*Ai';
Iss = intersect(Ipos,Ineg);
Aeplus = Ae;
Aiplus = [-W(:,Ipos)';W(:,Ineg)'];
nb = 0;
G = [];
%
% Call the double description method to compute the generators of an 
% appropriate cone
[G,nb] = doubledescLI(Aeplus,Aiplus);

G = Z'*G;
%
for i=1:size(G,2)
	G(:,i) = G(:,i)/norm(G(:,i));
end
%
% Identification of the lineality space
%
Is = [];
Ir = [];
V = G'*G;
%
Iv = 1:size(V,1);
%
while ~isempty(Iv)
	i = Iv(1);
	V = G(:,i)'*G;	
	Jopp = find(V==-1);
	if ~isempty(Jopp)
		Is = [Is i];
		Ir = [Ir Jopp];
		Iv = setdiff(Iv,Jopp);
	end
	Jsame = find(V==1);
	Jsame = setdiff(Jsame,[i]);
	if ~isempty(Jsame) 
		Iv = setdiff(Iv,Jsame);	
		Ir = [Ir Jsame];
	end
	Iv = setdiff(Iv,[i]);
end
%
Ic = setdiff(1:size(G,2),[Is Ir]);
Ys = G(:,Is);
Yc = G(:,Ic);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
