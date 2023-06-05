function [Gb,nb] = doubledescLI(B,A)
%
%	Computation of generators of a pointed cone by the double description 
%	method revisited (Fukuda, Prodon - 1997).
%
%	GB = doubledescLI(B,A) computes a minimal generating set for the cone in 
%	R^d defined as {x | B*x = 0, A*x >= 0}.
% 
%	Inputs:
%
%	B: matrix of full row rank m < d representing linear equalities
%
%	A: matrix of linear inequalities
%
%	Outputs:
%
%	Gb: Generating set of the cone
%
%	nb: number of generators
%
%	C. W. Royer - March 26, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
d = size(A,2);
Ed = eye(d);
%
% Initial set
if ~isempty(B)
	ZB = null(B);
	kb = rank(B);
else
	kb = 0;
	ZB = Ed;
end
R = [ZB -ZB]; % This is a generating set for null(B);
nb = size(R,2);
%
for j=1:size(A,1)
%
	v = A(j,:)*R;
	if j==1
		ka = kb;
		J = [];
	else
		J = 1:(j-1);
		ka = kb + rank(A(J,:));
	end
%
	Iplus = find(v>0);
	Izero = find(v==0);
	Iminus = setdiff(1:nb,[Iplus Izero]);
%
%	Finding the adjacent vectors
%
	Adj = 0;
	Rnew = [];
%
	for ip=1:length(Iplus)
		for im=1:length(Iminus)
			r1 = R(:,Iplus(ip));
			r2 = R(:,Iminus(im));
			ar1 = A(1:j,:)*r1;
			ar2 = A(1:j,:)*r2;
			if ~isempty(J)
				z1 = find(ar1(J,:)==0);
				if ~isempty(z1)
					k1 = rank([B;A(z1,:)]);
				else
					k1 = kb;
				end
				z2 = find(ar2(J,:)==0);
				if ~isempty(z1)
					k2 = rank([B;A(z2,:)]);
				else
					k2 = kb;
				end
				iz = intersect(z1,z2);
				if ~isempty(iz)
					kz = rank([B;A(iz,:)]);
				else
					kz = kb;
				end
			else
				k1=kb;
				k2=kb;
				kz=kb;
			end
			if (k1==(ka-1)) && (k2==(ka-1)) && (kz==(ka-2))
%				Adjacent directions
				r3 = ar1(j)*r2-ar2(j)*r1;
				if (norm(r3)~=0)
					Rnew = [Rnew r3];
					Adj = Adj + 1;
				end
			end
		end
	end
%
%
%	Building the new set
%
	R = [R(:,Iplus) R(:,Izero) Rnew];
	nb = size(R,2);	
%
end
%
Gb = R;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
