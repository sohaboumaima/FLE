function D = randcolZsearch(Z,Ip,Im,nbgen)
%
%	Random choice of columns of a given matrix Z to obtain a set of 
%	directions.
%
%	D = randcolZsearch(Z,Ip,Im,nbgen) selects randomly nbgen elements among 
%	a subset of vectors in Z and -Z.
%
%	Inputs:
%
%	Z: Matrix/Basis of the vectors to be used
%
%   Ip: Indexes of the inactive columns of the matrix Z
%
%   Im: Indexes of the inactive columns of the opposite matrix -Z
%
%	nbgen: Number of directions to be generated
%
%	Output:
%
%	D: The set of generated directions
%
%	C. W. Royer - March 26, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
Y = [Z(:,Ip) -Z(:,Im)];
p = randperm(size(Y,2));
D = Y(:,p(1:nbgen));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
