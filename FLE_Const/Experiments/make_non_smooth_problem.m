function [p]= make_non_smooth_problem(pr, to_change)
% 
% pr = macup(problem);
p.x0 = pr.x0;
lambda = 1e-4;

switch to_change
    case 'bounds'
         p.Aeq = pr.Aeq;
         p.beq = pr.beq;
         p.Aineq = pr.Aineq;
         p.bineq = pr.bineq;
         p.objective = @(x)pr.objective(x);
         if ~isempty(pr.lb)
             p.objective = @(x)(p.objective(x) + lambda*norm(x-pr.lb, 1));
        
         end
         
         if ~isempty(pr.ub)
             p.objective = @(x)(p.objective(x) + lambda*norm(pr.ub-x, 1));
         end
         p.lb = [];
         p.ub = [];

    case 'linearequ'
         p.lb = pr.lb;
         p.ub = pr.ub;
         p.Aineq = pr.Aineq;
         p.bineq = pr.bineq;
         p.objective = @(x)pr.objective(x);
         if ~isempty(pr.Aeq)
             p.objective = @(x)(p.objective(x) + lambda*norm(pr.Aeq*x-pr.beq, 1));
         end
         p.Aeq = [];
         p.beq = [];

    case 'linearinequ'
         p.lb = pr.lb;
         p.ub = pr.ub;
         p.Aeq = pr.Aeq;
         p.beq = pr.beq;
         p.objective = @(x)pr.objective(x);
         if ~isempty(pr.Aineq)
             p.objective = @(x)(p.objective(x) + lambda*norm(pr.Aineq*x-pr.bineq, 1));
         end
         p.Aineq = [];
         p.bineq = [];
end