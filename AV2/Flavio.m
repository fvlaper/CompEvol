function [ind_best, aval_best, IGD_best, IGD_m, ind_pior, aval_pior, IGD_pior] = Flavio(n_aval, problem, n_obj, n_exec)

if problem == 1
    func = @dtlz1;
    k = 5;
    precisao = 0.5;
else
    func = @dtlz2;
    k = 10;
    precisao = 0.5;
end

nvar = n_obj + k - 1;

if n_obj == 3
    npop = 92;
else
    npop = 212;
end

xi = 1; xf = xi+nvar-1; L.COLX  = xi:xf;
fi = xf+1; ff = fi+n_obj-1; L.COLF  = fi:ff;

name = sprintf('dtlz%d_%dd',problem,n_obj);
pareto = load(name);
pareto = pareto.fronteiraReal;

IGD_best = realmax;
IGD_pior = realmin;
IGD_m = 0;

for i = 1:n_exec
    display(['Avaliando ', num2str(i)]);
    ps = moea_3w(n_aval, npop, nvar, n_obj, precisao, func);
    ind = ps(:,L.COLX);
    aval = ps(:,L.COLF);
    igd = calcIgd(aval,pareto);
    IGD_m = IGD_m + igd;
    if igd < IGD_best
        IGD_best = igd;
        ind_best = ind;
        aval_best = aval;
    end
    if igd > IGD_pior
        IGD_pior = igd;
        ind_pior = ind;
        aval_pior = aval;
    end
    
    %----------------
    outname = sprintf('%s_%d',name,i);
    med = IGD_m / i;
    save(outname,'ind_best','aval_best','IGD_best','med','ind_pior','aval_pior','IGD_pior');
    %----------------
end

IGD_m = IGD_m / n_exec;

end
