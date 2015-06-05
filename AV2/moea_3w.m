function [ ps ] = moea_3w( ncal, nvar, no )
%MOEA_3W Multi Objective Evolutionary Algorithm - 3 Way
%   Esta fun��o implementa um algoritmo para otimiza��o de fun��es
%   Multi-objetivo. XXX Explicar o algoritmo.
%
%   Par�metros de entrada:
%     - ncal: n�mero m�ximo de c�lculos das fun��es objetivo;
%     - nvar: n�mero de vari�veis;
%     - no: n�mero de fun��es objetivo.
%
%   Par�metros de sa�da:
%     - ps: conjunto de Pareto.

% Layout de array de popula��o: array contendo uma linha por
% indiv�duo. Para cada linha:
%   - Vari�veis de decis�o: nvar colunas.
%   - Fun��es objetivo: no colunas.
%   - Valor agregado: uma coluna.
%   - Fronteira de Pareto: uma coluna.
xi = 1; xf = xi+nvar-1;
fi = xf+1; ff = fi+no-1;
ag = ff+1;
pt = ag+1;
nc = pt;

L.COLX  = xi:xf;
L.COLF  = fi:ff;
L.COLAG = ag;
L.COLPT = pt;
L.NC = nc;

% Dados iniciais da popula��o: faixa das vari�veis e 
% quantidade de indiv�duos.
xmin = 0; xmax = 1;
npop = 10;

% Estabelecimento da polula��o inicial.
pop = popinit(npop,xmin,xmax,L);

% Teste
pop = dtlz1(pop,nvar,no,L);
%pop = dtlz2(pop,nvar,no,L);

pop = pareto(pop,L);

% Retorno do resultado
ps = pop;
end

function [pop] = popinit (npop, xmin, xmax, L)
%POPINIT Gera��o da popula��o inicial.
%   Gera a popula��o inicial.
%
%   Par�metros de entrada:
%     - npop: n�mero de indiv�duos da popula��o;
%     - xmin: valor m�nimo de uma vari�vel;
%     - xmax: valor m�ximo de uma vari�vel;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - pop: array contendo a popula��o inicial.

% Aloca o array.
pop = zeros(npop, L.NC);

% Inicializa as vari�veis.
pop(:,L.COLX) = xmin * ones(npop,max(L.COLX)) + ...
                (xmax - xmin) * rand(npop,max(L.COLX));
end

function pop = pareto(pop,L)
%PARETO Calcula a fronteira de Pareto.
%   Para cada indiv�duo da polula��o, calcula a que fronteira de
%   Pareto ele pertence.
%   Inspirado em XXX.
%
%   Par�metros de entrada:
%     - pop: array contendo a popula��o inicial;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - pop: popula��o inicial com fronteira de Pareto atualizada.

npop = size(pop,1); % n�mero de indiv�duos
no = max(L.COLF) - min(L.COLF) + 1; % n�mero de fun��es objetivo
f = 1;  % f-�sima fronteira de Pareto

% Fronteira: vetor com uma estrutura por fronteira de Pareto:
%   n: n�mero de indiv�duos na fronteira;
%   ids: �ndices dos indiv�duos na fronteira.
indfront = struct('n',0, 'ids',zeros(1,npop));
front = repmat(indfront,1,npop);

% Domina��o: vetor com uma estrutura por indiv�duo i:
%   n: n�mero de indiv�duos que dominam i;
%   nd: n�mero de indiv�duos dominados por i;
%   d: �ndices dos indiv�duos dominados por i.
elemento = struct('n',0, 'nd',0, 'd',zeros(1,npop));
dominacao = repmat(elemento,1,npop);

% Compara cada indiv�duo com todos os demais para determinar
% a domina��o.
for i = 1:npop
    for j = 1:npop
        if i ~= j
            % Compara indiv�duos i e j
            difs = pop(i,L.COLF) - pop(j,L.COLF);
            menores = sum(difs <  0);
            iguais  = sum(difs == 0);
            maiores = sum(difs >  0);
            %display(['menor = ', num2str(menores), ' maior =  ', num2str(maiores), ' igual = ', num2str(iguais)]);
            
            if maiores == 0 && iguais ~= no
                % i domina j: nenhum maior, pelo menos um menor
                dominacao(i).nd = dominacao(i).nd + 1;
                dominacao(i).d(dominacao(i).nd) = j;
            elseif menores == 0 && iguais ~= no
                % j domina i: nenhum menor, pelo menos um diferente
                dominacao(i).n = dominacao(i).n + 1;
            end
        end        
    end
    
    % Trata os indiv�duos na primeira fronteira.
    if dominacao(i).n == 0
        front(f).n = front(f).n + 1;
        front(f).ids(front(f).n) = i;
        pop(i,L.COLPT) = 1;
    end
end

%for i = 1:npop
%    display([num2str(i), ': n = ', num2str(dominacao(i).n), ' nd = ', num2str(dominacao(i).nd)]);
%    display(dominacao(i).d);
%end

% Tratamento das demais fronteiras
while front(f).n ~= 0
    % Percorre os elementos da f-�sima fronteira, removendo-os
    % e fazendo os c�lculos para a (f+1)-�sima fronteira.
    for i = 1:front(f).n
        individuo = front(f).ids(i);
%        display(['Tratando indiv�duo ', num2str(individuo)]);
        % Trata os indiv�duos dominados pelo indiv�duo atual.
        for j = 1:dominacao(individuo).nd
            dominado = dominacao(individuo).d(j);
%            display([' Domina ', num2str(dominado)]);
            dominacao(dominado).n = dominacao(dominado).n - 1;
            if dominacao(dominado).n == 0
                front(f+1).n = front(f+1).n + 1;
                front(f+1).ids(front(f+1).n) = dominado;
                pop(dominado,L.COLPT) = f+1;
            end
        end
    end
    
%    display(['fronteira ', num2str(f), ': ', num2str(front(f).n)]);
%    display(front(f).ids);
    f = f + 1;
end

end

function pop = dtlz1 (pop,nvar,no,L)
%DTLZ1 Fun��o dtlz1.
%   Calcula as fun��es objetivo DTLZ1 para todos os indiv�duos
%   da popula��o.
%
%   Par�metros de entrada:
%     - pop: popula��o;
%     - nvar: n�mero de vari�veis;
%     - no: n�mero de fun��es objetivo;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - pop: popula��o com os valores das fun��es objetivo.

if no > nvar
    ME = MException('dtlz1:nvarError','Invalid number of variables');
    throw(ME);
end

% Colunas das vari�veis.
xi = min(L.COLX);
fi = min(L.COLF);
ff = max(L.COLF);
k = nvar - no + 1;

xg = pop(:,no:nvar);

% C�lculo
gx = 100 * (k + sum(((xg-0.5) .^ 2) - cos(20*pi*(xg-0.5)),2));

pop(:,fi) = 0.5 * (prod(pop(:,xi:(no-1)),2) .* (1+gx));

for i = 2:(no-1)
  pop(:,(fi+i-1)) = 0.5 * (prod((pop(:,xi:(no-i))),2) .* ...
                           (1 - pop(:,(no-i+1))) .* (1+gx));
end

pop(:,ff) = 0.5 * ((1 - pop(:,xi)) .* (1+gx));

end

function pop = dtlz2 (pop,nvar,no,L)
%DTLZ1 Fun��o dtlz1.
%   Calcula as fun��es objetivo DTLZ2 para todos os indiv�duos
%   da popula��o.
%
%   Par�metros de entrada:
%     - pop: popula��o;
%     - nvar: n�mero de vari�veis;
%     - no: n�emro de fun��es objetivo;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - pop: popula��o com os valores das fun��es objetivo.

if no > nvar
    ME = MException('dtlz2:nvarError','Invalid number of variables');
    throw(ME);
end

% Colunas das vari�veis.
xi = min(L.COLX);
fi = min(L.COLF);
ff = max(L.COLF);

xg = pop(:,no:nvar);

% C�lculo
gx = sum(((xg-0.5) .^ 2),2);

pop(:,fi) = (1+gx) .* prod(cos(0.5*pi*pop(:,xi:(no-1))),2);

for i = 2:(no-1)
    pop(:,(fi+i-1)) = (1+gx) .* prod(cos(0.5*pi*pop(:,xi:(no-i))),2) .* ...
                                sin(0.5*pi*pop(:,(no-i+1)));
end

pop(:,ff) = (1+gx) .* sin(0.5*pi*pop(:,xi));

end
