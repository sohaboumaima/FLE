function [T,cl,prob_max, hl] = perf_profile_dic(H,N,gate,logplot)
%     This subroutine produces a performance profile as described in:
%
%     Benchmarking Derivative-Free Optimization Algorithms
%     Jorge J. More' and Stefan M. Wild
%     SIAM J. Optimization, Vol. 20 (1), pp.172-191, 2009.
%
%     The latest version of this subroutine is always available at
%     http://www.mcs.anl.gov/~more/dfo/
%     The authors would appreciate feedback and experiences from numerical
%     studies conducted using this subroutine.
%
%     Performance profiles were originally introduced in
%     Benchmarking optimization software with performance profiles,
%     E.D. Dolan and J.J. More', 
%     Mathematical Programming, 91 (2002), 201--213.
%
%     The subroutine returns a handle to lines in a performance profile.
%
%       H contains a three dimensional array of function values.
%         H(f,p,s) = function value # f for problem p and solver s.
%       gate is a positive constant reflecting the convergence tolerance.
%       logplot=1 is used to indicate that a log (base 2) plot is desired.
%
%     Argonne National Laboratory
%     Jorge More' and Stefan Wild. January 2008.
%
% Modified by Oumaima Sohab to use dictionaries instead of matrices

ns = length(H);
firstsolv = H(1);
firstsolv = firstsolv{:};
problems= firstsolv.keys;

np = length(problems);

prob_min = Inf * ones(np, 1);
prob_max = -Inf * ones(np, 1);

for solver=H 
    dict = solver{:};
    i = 1;
    while i<=np
            
        array = dict(problems(i));
        if (isempty(array{:}))
            array ={[prob_max(i)]};
        end
        prob_min(i) = min(prob_min(i), min(array{:}));
        prob_max(i) = max(prob_max(i), max(array{:}));

        i = i+1;

    end
end

% For each problem and solver, determine the number of evaluations
% required to reach the cutoff value
T = zeros(np,ns);

for s = 1:ns
  dicH = H(s);
  dicH = dicH{:};
  dicN = N(s);
  dicN = dicN{:};

  cl = [];

  for p = 1:np
      
      cutoff = prob_min(p) + gate*(prob_max(p) - prob_min(p));
      cl = [cl, cutoff];
      arr = dicH(problems(p));
      ind_f = find(arr{:} <= cutoff,1);

    if (isempty(ind_f))
      T(p,s) = NaN;
    else
      prob = dicN(problems(p));
      prob = prob{:};
      nfval  = prob(ind_f);
      T(p,s) = nfval;
    end
  end
end


% Other colors, lines, and markers are easily possible:
colors  = ['b' 'r' 'k' 'm' 'c' 'g' 'y'];   lines   = {'-' '-.' '--'};
markers = [ 's' 'o' '^' 'v' 'p' '<' 'x' 'h' '+' 'd' '*' '<' ];

if (nargin < 3); logplot = 0; end

% Compute ratios and divide by smallest element in each row.
r = T./repmat(min(T,[],2),1,ns);

% Replace all NaN's with twice the max_ratio and sort.
max_ratio = max(max(r));
r(isnan(r)) = 2*max_ratio;


r = sort(r);


% Plot stair graphs with markers.
hl = zeros(ns,1);

for s = 1:ns

    [xs,ys] = stairs(r(:,s),(1:np)/np);
    
    % Only plot one marker at the intercept
    if (xs(1)==1)
        vv = find(xs==1,1,'last');
        xs = xs(vv:end);   ys = ys(vv:end);
    end

    sl = mod(s-1,3) + 1; sc = mod(s-1,7) + 1; sm = mod(s-1,12) + 1;
    option1 = [char(lines(sl)) colors(sc) markers(sm)];
    if (logplot)
        hl(s) = semilogx(xs,ys,option1);
    else
        hl(s) = plot(xs,ys,option1);
    end
    hold on;
end

% Axis properties are set so that failures are not shown, but with the
% max_ratio data points shown. This highlights the "flatline" effect.
if (logplot) 
  axis([1 1.1*max_ratio 0 1]);
  twop = floor(log2(1.1*max_ratio));
  set(gca,'XTick',2.^[0:twop])
else
  axis([1 1.1*max_ratio 0 1]);
end