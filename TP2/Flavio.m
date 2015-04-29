function [ x, f, g, h ] = Flavio( ncal, nvar )
%FLAVIO Trabalho 2 de Computação Evolucionária - Flávio Laper.
%   Esta função implementa um algoritmo genético para minimizar a
%   função de Rastrigin para nvar variáveis com restrições.
%   
%   Parâmetros de entrada:
%       - ncal: número máximo de chamadas da função de cálculo da
%         fitness do problema;
%       - nvar: número de variáveis.
%
%   Parâmetros de saída:
%       - x: vetor das variáveis de decisão do melhor indivíduo;
%       - f: melhor função objetivo;
%       - g: vetor restrição de desigualdade avaliado no ponto x;
%       - h: vetor restrição de igualdade avaliado no ponto x.

pop = popinit(3,3,-5.12,5.12);

x = [1,1,1];
f = 2;
g = pop;
%h = rastrigin(pop);
%h = restr_desig(pop);
h = restr_ig(pop);

end

function [pop] = popinit ( npop, nvar, xmin, xmax )
%POPINIT Geração da população inicial.
%   Gera a população com npop indivíduos e nvar variáveis. Cada
%   variável tem valores limitados à faixa [xmin,xmax].
%
%   Parâmetros de entrada:
%       - npop: número de indivíduos da população;
%       - nvar: número de variáveis por indivíduo;
%       - xmin: valor mínimo de uma variável;
%       - xmax: valor máximo de uma variável.
%
%   Parâmetros de saída:
%       - pop: array contendo a população (npop linhas, nvar colunas). 

pop = xmin * ones(npop,nvar) + (xmax - xmin) * rand(npop,nvar);

end

function [f] = rastrigin (pop)
%RASTRIGIN Função de Rastrigin.
%   Calcula a função de Rastrigin para todos os indivíduos da população.
%
%   Parâmetros de entrada:
%       - pop: array contendo a população (unm indivíduo por linha)
%
%   Parâmetros de saída:
%       - f: valor da função de Rastrigin para cada indivíduo (vetor
%         coluna com npop linhas).

f = 10 * size(pop,2) + sum(pop .^ 2 - 10 * cos(2 * pi * pop), 2);

end

function [g] = restr_desig (pop)
%RESTR_DESIG Função de restrição de desigualdade.
%   Calcula o valor da restrição de desigualdade para cada variável
%   de todos os indivíduos da população.
%
%   Parâmetros de entrada:
%       - pop: array contendo a população (unm indivíduo por linha)
%
%   Parâmetros de saída:
%       - g: array contendo a função de restrição (npop linhas,
%         nvar colunas). 

g = sin(2 * pi * pop) + 0.5;

end

function [h] = restr_ig (pop)
%RESTR_IG Função de restrição de igualdade.
%   Calcula o valor da restrição de igualdade para cada variável
%   de todos os indivíduos da população.
%
%   Parâmetros de entrada:
%       - pop: array contendo a população (unm indivíduo por linha)
%
%   Parâmetros de saída:
%       - g: array contendo a função de restrição (npop linhas,
%         nvar colunas). 

h = cos(2 * pi * pop) + 0.5;

end