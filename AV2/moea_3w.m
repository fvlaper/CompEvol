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
xi = 1; xf = xi+nvar-1;
fi = xf+1; ff = fi+no-1;
nc = ff;

L.COLX = xi:xf;
L.COLF = fi:ff;
L.NC = nc;

% Dados iniciais da popula��o: faixa das vari�veis e 
% quantidade de indiv�duos.
xmin = 0; xmax = 1;
npop = 10;

% Estabelecimento da polula��o inicial.
pop = popinit(npop,xmin,xmax,L);

% Teste
pop = dtlz1(pop,nvar,no,L);

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

function pop = dtlz1 (pop,nvar,no,L)
%DTLZ1 Fun��o dtlz1.
%   Calcula as fun��es objetivo DTLZ1 para todos os indiv�duos
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
    ME = MException('dtlz1:nvarError','Invalid number of variables');
    throw(ME);
end

% Colunas das vari�veis.
xi = min(L.COLX);
xf = max(L.COLX);
k = nvar - no + 1;
xg = pop(:,no:nvar);

gx = 100 * (k + sum(((xg-0.5) .^ 2) - cos(20*pi*(xg-0.5)),2));
display(gx);

pop(:,L.COLF) = 1;
end

