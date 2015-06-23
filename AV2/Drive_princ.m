n_aval = 500;
problem = 2; % 1 = DTLZ1; 2 = DTLZ2;
n_obj = 3;
n_exec = 3;

[ind_best, aval_best, IGD_best, IGD_m, ind_pior, aval_pior, IGD_pior] = Flavio(n_aval, problem, n_obj, n_exec);

display(n_aval);
display(ind_best);
display(aval_best);
display(IGD_best);
display(IGD_m);
display(ind_pior);
display(aval_pior);
display(IGD_pior);