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

pop = popinit(10,3,-5.12,5.12);
%pop = [pop fitness(@fitness_desc, @rastrigin, @rastrigin_le, @rastrigin_eq, 10, 10, 1.2, pop)];

[x, f, g, h] = ga(pop,1,1, @fitness_desc, @rastrigin, @rastrigin_le, @rastrigin_eq, 10, 10, 1.2);

%g = rastrigin(pop);
%g = fmod(@rastrigin, @rastrigin_le, @rastrigin_eq, 10, 10, pop);
%g = pop;
%h = [pop fitness(@fitness_desc, @rastrigin, @rastrigin_le, @rastrigin_eq, 10, 10, 1.2, pop)];
%h = restr_desig(pop);
%h = rastrigin_eq(pop);
%h = roleta(pop);

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

function [d,fo,gleo,geqo] = fitness(ft, f, gle, geq, r, s, n, pop)
%FITNESS Função de fitness.
%   Calcula a fitness dos indivíduos da população. Esta função é
%   apenas um "stub" que recebe a função de fitness específica a
%   utilizar no cálculo.
%
%   Parâmetros de entrada:
%       - ft: handle da função de fitness;
%       - f: handle da função de otimização original;
%       - gle: handle da função de restrição de desigualdade;
%       - geq: handle da função de restrição de igualdade;
%       - r: parâmetro de penalidade para a restrição de desigualdade;
%       - s: parâmetro de penalidade para a restrição de igualdade;
%       - n: parâmetro de escalonamento da função de fitness;
%       - pop: array contendo a polulação (um indivíduo por linha).
%
%   Parâmetros de saída:
%       - d: valor da fitness para cada indivíduo (vetor coluna);
%       - fo: valor da função de otimização original
%         para cada indivíduo (vetor coluna);
%       - gleo: valores da função de restrição de desigualdade
%         (npop linhas, nvar colunas);
%       - geqo: valores da função de restrição de igualdade
%         (npop linhas, nvar colunas). 

[d,fo,gleo,geqo] = ft(f, gle, geq, r, s, n, pop);
end

function [d,fo,gleo,geqo] = fitness_inv(f, gle, geq, r, s, n, pop)
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
%       - d: valor da fitness para cada indivíduo (vetor coluna);
%       - fo: valor da função de otimização original
%         para cada indivíduo (vetor coluna);
%       - gleo: valores da função de restrição de desigualdade
%         (npop linhas, nvar colunas);
%       - geqo: valores da função de restrição de igualdade
%         (npop linhas, nvar colunas). 

[h,fo,gleo,geqo] = fmod(f, gle, geq, r, s, pop);
d = 1 ./ (h - (min(h) - (10 ^ -n)));
end

function [d,fo,gleo,geqo] = fitness_desc(f, gle, geq, r, s, k, pop)
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
%       - d: valor da fitness para cada indivíduo (vetor coluna);
%       - fo: valor da função de otimização original
%         para cada indivíduo (vetor coluna);
%       - gleo: valores da função de restrição de desigualdade
%         (npop linhas, nvar colunas);
%       - geqo: valores da função de restrição de desigualdade
%         (npop linhas, nvar colunas). 

[h,fo,gleo,geqo] = fmod(f, gle, geq, r, s, pop);
cmax = k * sum(h) / size(h,1);
d = max(zeros(size(h)), (cmax-h));
end

function [h,fo,gleo,geqo] = fmod (f, gle, geq, r, s, pop)
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
%         para cada indivíduo (vetor coluna);
%       - fo: valor da função de otimização original
%         para cada indivíduo (vetor coluna).
%       - gleo: valores da função de restrição de desigualdade
%         (npop linhas, nvar colunas);
%       - geqo: valores da função de restrição de igualdade
%         (npop linhas, nvar colunas). 

fo = f(pop);
gleo = gle(pop);
geqo = geq(pop);

h = fo ...                                             % func. original
  + r * sum((max(zeros(size(pop)), gleo) .^ 2) ,2) ... % desigualdades
  + s * sum(geqo .^ 2, 2);                             % igualdades
end

function [ x, f, g, h ] = ga(pop,ngen,nelite,ft,f,gle,geq,r,s,n)
%GA Função principal do algoritmo genético.
%   Executa o algoritmo genético (ga) para otimização da função objetivo.
%
%   Parâmetros de entrada:
%       - pop: array contendo a população inicial;
%       - ngen: número de gerações a processar;
%       - nelite: número de elementos da elite (preservados para a
%         próxima geração);
%       - ft: handle da função de fitness;
%       - f: handle da função de otimização original;
%       - gle: handle da função de restrição de desigualdade;
%       - geq: handle da função de restrição de igualdade;
%       - r: parâmetro de penalidade para a restrição de desigualdade;
%       - s: parâmetro de penalidade para a restrição de igualdade;
%       - n: parâmetro de escalonamento da função de fitness.
%
%   Parâmetros de saída:
%       - x: vetor das variáveis de decisão do melhor indivíduo;
%       - f: melhor função objetivo;
%       - g: vetor restrição de desigualdade avaliado no ponto x;
%       - h: vetor restrição de igualdade avaliado no ponto x.

% População de pais para a próxima geração. Array contendo uma linha
% por indivíduo da população. Para cada linha:
%   - colunas 1 - nvar: variáveis de decisão;
%   - colunas (nvar+1) - (2*nvar): restrições de desigualdade;
%   - colunas (2*nvar+1) - (3*nvar): restrições de igualdade;
%   - coluna (3*nvar+1): função objetivo;
%   - coluna (3*nvar+2): fitness.
columnsX = 1:size(pop,2);
columnsG = size(pop,2)+1 : 2*size(pop,2);
columnsH = 2*size(pop,2)+1 : 3*size(pop,2);
columnF = 3*size(pop,2)+1;
columnFit = 3*size(pop,2)+2;
pais = [pop zeros(size(pop)) zeros(size(pop)) zeros(size(pop,1),2)];

% Cálculos iniciais para a primeira população de pais
[fit,fo,gleo,geqo] = fitness(ft, f, gle, geq, r, s, n, pop);
pais(:,columnsG) = gleo;
pais(:,columnsH) = geqo;
pais(:,columnF) = fo;
pais(:,columnFit) = fit;
display(pais);

x = pais(1,columnsX);
f = pais(1,columnF);
g = pais(1,columnsG);
h = pais(1,columnsH);

end

%function [pais] = roleta(pop)
%ROLETA Seleção pelo método da roleta.
%   Seleciona os pais para a próxima geração pelo método da roleta.
%
%   Parâmetros de entrada:
%       - pop: array contendo a população (supõe-se que a última
%         coluna contém a fitness do indivíduo).
%
%   Parâmetros de saída:
%       - pais: população de pais selecionados para os procedimentos
%         genéticos seguintes.

%npop = size(pop,1);                    % tamanho da população
%fitcol = size(pop,2);                  % coluna com o valor da fitness
%totfitness = sum(pop(:,fitcol));       % fitness total
%pais = zeros(size(pop));               % pais selecionados

% Loop de seleção
%for n = 1:npop
%    r = totfitness * rand(); % valor do lance    
%    i = 1;                   % índice do indivíduo selecionado
%    s = pop(1,fitcol);       % soma acumulada
%    while s < r
%        i = i + 1;
%        s = s + pop(i,fitcol);
%    end
    
%    pais(n,:) = pop(i,:);
%    display(i);
        
%end

%end

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
