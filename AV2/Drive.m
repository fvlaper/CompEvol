    % dtlz1
ncal = 3680;
npop = 92;
nvar = 7;
no = 3;
%nvar = 9;
%no = 5;
precisao = 0.5;

ps = moea_3w(ncal,npop,nvar,no,precisao,@dtlz1);

display(ps);

% dtlz2
%{
ncal = 23000;
npop = 92;
nvar = 12;
no = 3;
%nvar = 14;
%no = 5;
precisao = 0.5;

ps = moea_3w(ncal,npop,nvar,no,precisao,@dtlz2);

display(ps);
%}
