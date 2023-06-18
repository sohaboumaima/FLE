function [functions_values, nfun_eval]= run_function(problems, solver, set)


[list_problems, ~]  = get_list_problems(problems);

num_runs  = 1;
nfun_eval = dictionary;
exitflags = dictionary;
functions_values = dictionary;

for i= 1:length(list_problems)
    
    prob          = list_problems(i);  
    p             = prob{:};
    disp(['problem ' int2str(i) ': ' char(p)])
    pr = macup(p);
%     if pr.numlb || pr.numub
%         pr = make_non_smooth_problem(pr, 'linearequ');
%         pr = make_non_smooth_problem(pr, 'linearinequ');
        if (isempty(pr.lb))
            pr.lb = -Inf *ones(size(x0));
        end
    
        if (isempty(pr.ub))
            pr.ub = Inf *ones(size(x0));
        end
        x0 = pr.x0;
        lb= pr.lb;
        ub = pr.ub;
        Aeq = pr.Aeq;
        beq = pr.beq;
        Aineq = pr.Aineq;
        pr.ubineq = pr.bineq;
        ubineq = pr.bineq;
        if isempty(Aineq)
            pr.lbineq = [];
        else
            pr.lbineq = -Inf *ones(size(Aineq, 1), 1);
        end
        lbineq = pr.lbineq;
        
    
        presbounds = ~isempty(find(ub < 1e20, 1)) || ~isempty(find(lb > -1e20, 1));
        if (isempty(Aineq))
            preslineq = 0;
        else
            preslineq = 1;
        end
        
        tol_feas = 1e-8;
        if (~isfeasible(x0,presbounds, preslineq, lb, ub, Aeq, beq, Aineq, lbineq, ubineq, tol_feas))
        
            x0 = project_on_feasible(x0, pr);
        end

%         if i==17
%             x0=[0.5; 1e-12; 1e-12; 0.5];
%         end
    
        max_eval = 2000*length(x0);

        for j=1:num_runs
    
        switch solver
    
            case 'nomad'
                if (j==1)
    % 
    %             bad = [15,17,18,19, 26, 34, 37, 38, 39];
                bad = [];
    
                if any(bad==i)
                    func_val = [];
                    nfeval = [];
                    exitflag = -1;
                else
    
                    [fx0, bbtype] = fun_nomad(@(x)fun(@pr.objective,x), x0, Aeq, beq, Aineq, lbineq, ubineq);
        
                    if (length(bbtype) > 1)
        
                        opts = nomadset('bb_output_type',bbtype);
                    else
                        opts = nomadset('bb_output_type',['OBJ']);
                    end
    %                 opts.display_stats = 1;
                    opts.max_bb_eval = max_eval;
                    opts.model_search = 0;
                    opts.stats_file = './output.txt';
                    opts.history_file = './iterates.txt';
                    [~,~,exitflag,~,~] = nomad(@(x)fun_nomad(@(x)fun(@pr.objective,x), x, Aeq, beq, Aineq, lbineq, ubineq), x0, lb, ub,opts);
                    fileID = fopen('./output.0.txt','r');
                    formatSpec = '%d %f';
                    sizeA = [2 Inf];
                    feval = fscanf(fileID,formatSpec,sizeA);
                    fclose(fileID);
                    if (isempty(feval))
                        func_val = [];
                        nfeval = [];
                    else
                        func_val = feval(2,:);
                        nfeval   = feval(1,:);
                        fileID_iterates = fopen('/home/ouma/Documents/rBFGS/constrained_LF/New_code/iterates.0.txt','r');
                        n = length(pr.x0);
                        iterates = fscanf(fileID_iterates, repmat('%f ', [1, length(fx0)+n]), [length(fx0)+n Inf]);
                        fvalues = [];
                        for k=1:length(feval(1,:))
                            val = func_val(k);
                            closest = min(abs(iterates(n+1, :) - val));
                            ind = find(abs(iterates(n+1, :) - val) == closest(1));
                            x = iterates(1:n, ind(1));
                            fvalues = [fvalues, pr.objective(x)];
                        end
                    end
                   
    
                end
                end
    
            case 'dspfd'
                [~,~,exitflag, output_ds] = dspfd(@(x,y)fun(@pr.objective,x,y),x0,1,2,lb,ub, Aeq, beq, Aineq, lbineq, ubineq);
                hist = output_ds.histF;
                func_val = output_ds.truefunc;
                nfeval = hist(:,2);
    
            case 'fle'
                try       
                    [~, ~, output_fe]= FullLowEval(@(x,y)fun(@pr.objective,x,y),x0,1,2,pr,max_eval,tol_feas);
                    func_val = output_fe.truefunc;
                    nfeval = output_fe.funcCounts;
                    exitflag = output_fe.flag;
                catch
                    func_val = [];
                    nfeval = [];
                    exitflag = -1;
                end
            case 'fl'
                if (j==1)
                    try
                        [~, ~, output_fe] = FullEval(@(x,y)fun(@pr.objective,x,y),x0, pr, max_eval);
                        func_val = output_fe.truefunc;
                        nfeval   = output_fe.funcCounts;
                        exitflag = output_fe.flag;
                    catch
                        func_val = [];
                        nfeval = [];
                        exitflag = -1;
                    end
                end
        end
    
        functions_values([p '_' int2str(j)]) = {func_val};
        nfun_eval([p '_' int2str(j)]) = {nfeval};
        exitflags([p '_' int2str(j)]) = exitflag;
        end
%     end
end


file = [solver set];
save(file,'functions_values','nfun_eval', "exitflags")
