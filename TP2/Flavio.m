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
f = pop;
%g = rastrigin(pop);
g = fmod(@rastrigin, @rastrigin_le, @rastrigin_eq, 10, 10, pop);
h = fitness_desc(@rastrigin, @rastrigin_le, @rastrigin_eq, 10, 10, 1.2, pop);
%h = restr_desig(pop);
%h = rastrigin_eq(pop);

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

function [d] = fitness_inv(f, gle, geq, r, s, n, pop)
%FITNESS_INV Função de fitness pelo método de inversão.
%   Calcula a fitness dos indivíduos da população pelo método de
%   inversão (transformando um problema de minimização em maximização).
%
%   Parâmetros de entrada:
%       - f: handle da função de otimização original;
%       - gle: handle da função de restrição de desigualdade;
%       - geq: handle da função de restrição de igualdade;
%       - r: parâmetro de penalidade para a restrição de desigualdade;
%       - s: parâmetro de penalidade para a restrição de igualdade;
%       - n: parâmetro de escalonamento para a inversão;
%       - pop: array contendo a polulação (um indivíduo por linha).
%
%   Parâmetros de saída:
%       - d: valor da fitness para cada indivíduo (vetor coluna).

h = fmod(f, gle, geq, r, s, pop);
d = 1 ./ (h - (min(h) - (10 ^ -n)));
end

function [d] = fitness_desc(f, gle, geq, r, s, k, pop)
%FITNESS_DESC Função de fitness pelo método de deslocamento.
%   Calcula a fitness dos indivíduos da população pelo método de
%   deslocamento (transformando um problema de minimização em maximização).
%
%   Parâmetros de entrada:
%       - f: handle da função de otimização original;
%       - gle: handle da função de restrição de desigualdade;
%       - geq: handle da função de restrição de igualdade;
%       - r: parâmetro de penalidade para a restrição de desigualdade;
%       - s: parâmetro de penalidade para a restrição de igualdade;
%       - k: parâmetro de escalonamento para o deslocamento;
%       - pop: array contendo a polulação (um indivíduo por linha).
%
%   Parâmetros de saída:
%       - d: valor da fitness para cada indivíduo (vetor coluna).

h = fmod(f, gle, geq, r, s, pop);
cmax = k * sum(h) / size(h,1);
display(cmax);
d = max(zeros(size(h)), (cmax-h));
end

function [h] = fmod (f, gle, geq, r, s, pop)
%FMOD Função de otimização modificada.
%   Calcula o valor da função de otimização (por indivíduo) modificada
%   para um problema irrestrito.
%
%   Parâmetros de entrada:
%       - f: handle da função de otimização original;
%       - gle: handle da função de restrição de desigualdade;
%       - geq: handle da função de restrição de igualdade;
%       - r: parâmetro de penalidade para a restrição de desigualdade;
%       - s: parâmetro de penalidade para a restrição de igualdade;
%       - pop: array contendo a polulação (um indivíduo por linha).
%
%   Parâmetros de saída:
%       - h: valor da função de otimização modificada
%         para cada indivíduo (vetor coluna).

h = f(pop) ...                                             % func. original
  + r * sum((max(zeros(size(pop)), gle(pop)) .^ 2) ,2) ... % desigualdades
  + s * sum(geq(pop) .^ 2, 2);                             % igualdades
end

%------------------------------------------------------------------------%
% Função de Rastrigin e restrições
%------------------------------------------------------------------------%

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

function [g] = rastrigin_le (pop)
%RASTRIGIN_LE Função de restrição de desigualdade - Rastrigin.
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

function [h] = rastrigin_eq (pop)
%RASTRIGIN_EQ Função de restrição de igualdade - Rastrigin.
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
