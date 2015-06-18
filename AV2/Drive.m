% dtlz1
%{
ncal = 13;
npop = 10;
nvar = 7;
no = 3;
%nvar = 9;
%no = 5;
precisao = 0.5;

ps = moea_3w(ncal,npop,nvar,no,precisao,@dtlz1);

display(ps);
%}

% dtlz2
ncal = 13;
npop = 10;
nvar = 12;
no = 3;
%nvar = 14;
%no = 5;
precisao = 0.5;

ps = moea_3w(ncal,npop,nvar,no,precisao,@dtlz2);

display(ps);