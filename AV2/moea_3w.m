function [ ps ] = moea_3w( ncal, nvar, no )
%MOEA_3W Multi Objective Evolutionary Algorithm - 3 Way
%   Esta função implementa um algoritmo para otimização de funções
%   Multi-objetivo. XXX Explicar o algoritmo.
%
%   Parâmetros de entrada:
%     - ncal: número máximo de cálculos das funções objetivo;
%     - nvar: número de variáveis;
%     - no: número de funções objetivo.
%
%   Parâmetros de saída:
%     - ps: conjunto de Pareto.

% Layout de array de população: array contendo uma linha por
% indivíduo. Para cada linha:
%   - Variáveis de decisão: nvar colunas.
xi = 1; xf = xi+nvar-1;
fi = xf+1; ff = fi+no-1;
nc = ff;

L.COLX = xi:xf;
L.COLF = fi:ff;
L.NC = nc;

% Dados iniciais da população: faixa das variáveis e 
% quantidade de indivíduos.
xmin = 0; xmax = 1;
npop = 10;

% Estabelecimento da polulação inicial.
pop = popinit(npop,xmin,xmax,L);

% Retorno do resultado
ps = pop;
end

function [pop] = popinit (npop, xmin, xmax, L)
%POPINIT Geração da população inicial.
%   Gera a população inicial.
%
%   Parâmetros de entrada:
%     - npop: número de indivíduos da população;
%     - xmin: valor mínimo de uma variável;
%     - xmax: valor máximo de uma variável;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - pop: arrau contendo a população inicial.

% Aloca o array.
pop = zeros(npop, L.NC);

% Inicializa as variáveis.
pop(:,L.COLX) = xmin * ones(npop,max(L.COLX)) + ...
                (xmax - xmin) * rand(npop,max(L.COLX));
end

