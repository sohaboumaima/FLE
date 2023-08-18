% Include all the files in the parent directiory to the path
currentPath = pwd;
parts = strsplit(currentPath, filesep);
addpath(strjoin(parts(1:end-1), filesep))

solvers = ["dspfd", "fle", "nomad", "fl" ]; % List of Solvers
problems =  ["bound_constrained.csv","linearequ_const_copy.csv", "linearinequ_const.csv"]; % List of problems
set  = ["_boundconst_10", "_boundconst_10", "_boundconst_10", "_boundconst_10"];

fval = {};
neval = {};
j = 1;
for i=solvers
    run_function(problems, char(i), char(set(j))) % Solve problems with Solver i
    load([char(i) char(set(j)) '.mat'])
    fval{end+1} = functions_values;
    neval{end+1} =  nfun_eval;
    j = j+1;
end

% Run performance profiles 
figure(1)
gate = 1e-3;
logplot = 1;
perf_profile_dic(fval,neval,gate,logplot)
title('\tau = 10^{-3}')
xlabel('Ratio of function calls (log scale)')
ylabel('Fraction of problems solved')
legend([ "dspfd",  "constFLE", "NOMAD","constBFGS"])
% 
figure(2)
% 

gate = 1e-5;
logplot = 1;
perf_profile_dic(fval,neval,gate,logplot)
xlabel('Ratio of function calls (log scale)')
ylabel('Fraction of problems solved')
title('\tau = 10^{-5}')
legend([ "dspfd",  "constFLE", "NOMAD","constBFGS"])
