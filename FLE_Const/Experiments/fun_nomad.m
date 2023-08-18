function [out, bbtype] = fun_nomad(f, x, Aeq, beq, Aineq, lbineq, ubineq)

out = [f(x)];
if ~isempty(Aeq)
        c1 = -Aeq*x + beq;
        c2 = -beq + Aeq*x;
        out = [out, c1', c2'];
end
if ~isempty(Aineq)
    c3 = Aineq*x - ubineq;
    c4 = lbineq - Aineq*x;
    ind3 = find(c3 <= -Inf);
    c3(ind3) = [];
    ind4 = find(c4 <= -Inf);
    c4(ind4) = [];
    out = [out, c3', c4']; 
end

bbtype = ['OBJ '];

for ind=1:length(out)-1

    bbtype = [bbtype, 'PB '];

end

