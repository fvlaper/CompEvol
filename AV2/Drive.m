% dtlz1
%{
ncal = 13;
npop = 10;
nvar = 7;
no = 3;
%nvar = 9;
%no = 5;
resolucao = 5;

ps = moea_3w(ncal,npop,nvar,no,resolucao,@dtlz1);

display(ps);
%}

% dtlz2
ncal = 13;
npop = 10;
nvar = 12;
no = 3;
%nvar = 14;
%no = 5;
resolucao = 5;

ps = moea_3w(ncal,npop,nvar,no,resolucao,@dtlz2);

display(ps);