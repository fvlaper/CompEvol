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

% Estabelecimento dos parametros iniciais
npop = 10;                % número de indivíduos na população.
ngen = (ncal / npop) - 1; % número de gerações. Calculado a partir
                          % do número máximo de cálculos da função de
                          % fitness, considerando npop cálculos por
                          % geração, mais npop cálculos para a população
                          % inicial.
                          
f   = @rastrigin;         % Handle da função objetivo.
gle = @rastrigin_le;      % Handle da função de restrição de desigualdade.
geq = @rastrigin_eq;      % Handle da função de restrição de igualdade.
xmin = -5.12;             % Valores de contorno para as variáveis
xmax = 5.12;              % de decisão.
                          
% Estabelecimento da população inicial.
pop = popinit(npop,nvar,xmin,xmax);

% Execução do algoritmo genético para otimização.
[x, f, g, h] = ga(pop, ngen, f, gle, geq);

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

function [ x, f, g, h ] = ga(pop,ngen,f,gle,geq)
%GA Função principal do algoritmo genético.
%   Executa o algoritmo genético (ga) para otimização da função objetivo.
%
%   Parâmetros de entrada:
%       - pop: array contendo a população inicial;
%       - ngen: número de gerações a processar;
%       - f: handle da função de otimização original;
%       - gle: handle da função de restrição de desigualdade;
%       - geq: handle da função de restrição de igualdade.
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
nvar = size(pop,2);
columnsX  = 1:nvar;
columnsG  = nvar+1 : 2*nvar;
columnsH  = 2*nvar+1 : 3*nvar;
columnF   = 3*nvar+1;
%columnFit = 3*nvar+2;
pais = [pop zeros(size(pop)) zeros(size(pop)) zeros(size(pop,1),2)];

%Parâmetros do algoritmo
ft = @fitness_desc; % handle da função de fitness.
fp = @fitpar_desc;  % handle da função de parâmetros da fitness.

nelite = 1;         % número de elementos da elite (preservados para a
                    % próxima geração).

% Parâmetros de penalidade para a função objetivo modificada para
% problema sem restrições. A penalidade pelo desrespeito às restrições
% aumenta a cada geração.
r = 1 : ((10^8 - 1) / ngen) : 10^8;
s = 1 : ((10^8 - 1) / ngen) : 10^8;

% Parâmetros da função de fitness
fitpar = fp(ngen);

% Cálculos iniciais para a primeira população de pais
[pais] = fitness(ft, f, gle, geq, r(1), s(1), fitpar(1), pais, nvar);
display(pais);
pais = popSort(pais,nvar);
display(pais);

% Inicialização da população de filhos
filhos = zeros(size(pais,1)-nelite, size(pais,2));
display(filhos);

% Loop de evolução
for g = 1 : ngen
    % Guarda os indivíduos da elite
    elite = pais(1 : nelite, :);
    display(elite);
end

x = pais(1,columnsX);
f = pais(1,columnF);
g = pais(1,columnsG);
h = pais(1,columnsH);

end

function [pop] = popSort(pop,nvar)
%POPSORT Ordena o array de populaçao.
%   Ordena o array de população do melhor para o pior indivíduo.
%
%   Parâmetros de entrada:
%       - pop: array de população;
%       - nvar: número de variáveis.
%
%   Parâmetros de saída:
%       - pop: array de população ordenado.

columnFit = 3*nvar+2; % Coluna do valor de fitness.

pop = sortrows(pop, -columnFit);
end

function [pop] = fitness(ft, f, gle, geq, r, s, n, pop, nvar)
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
%       - pop: array contendo a população (um indivíduo por linha);
%       - nvar: número de variáveis.
%
%   Parâmetros de saída:
%       - pop: array de população com valores calculados.

[pop] = ft(f, gle, geq, r, s, n, pop, nvar);
end

function [pop] = fitness_inv(f, gle, geq, r, s, n, pop, nvar)
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
%       - pop: array contendo a polulação (um indivíduo por linha);
%       - nvar: número de variáveis.
%
%   Parâmetros de saída:
%       - pop: array de população com valores calculados.

columnFit = 3*nvar+2;

[h,pop] = fmod(f, gle, geq, r, s, pop, nvar);
d = 1 ./ (h - (min(h) - (10 ^ -n)));
pop(:,columnFit) = d;
end

function [par] = fitpar_inv(ngen)
%FITPAR_DESC Parâmetro para cálculo da fitness pelo método de inversão.
%   Retorna, para cada geração, o valor do parâmetro n para o cálculo
%   da fitness pelo método de deslocamento.
%
%   Parâmetros de entrada:
%       - ngen: número de gerações.
%
%   Parâmetros de saída:
%       - vetor com (ngen+1) valores do parâmetro (um para cada geração,
%         mais a população inicial).

par = 0.1 : ((1.0 - 0.1) / ngen) : 1.0;
end

function [pop] = fitness_desc(f, gle, geq, r, s, k, pop, nvar)
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
%       - pop: array contendo a polulação (um indivíduo por linha);
%       - nvar: número de variáveis.
%
%   Parâmetros de saída:
%       - pop: array de população com valores calculados.

columnFit = 3*nvar+2;

[h,pop] = fmod(f, gle, geq, r, s, pop, nvar);
cmax = k * sum(h) / size(h,1);
d = max(zeros(size(h)), (cmax-h));
pop(:,columnFit) = d;
end

function [par] = fitpar_desc(ngen)
%FITPAR_DESC Parâmetro para cálculo da fitness pelo método de deslocamento.
%   Retorna, para cada geração, o valor do parâmetro K para o cálculo
%   da fitness pelo método de deslocamento. Será utilizado um valor
%   fixo para todas as gerações.
%
%   Parâmetros de entrada:
%       - ngen: número de gerações.
%
%   Parâmetros de saída:
%       - vetor com (ngen+1) valores do parâmetro (um para cada geração,
%         mais a população inicial).

k = 1.2; % Valor fixo
par = k * ones(1,ngen+1);
end

function [h,pop] = fmod (f, gle, geq, r, s, pop, nvar)
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
%       - pop: array contendo a polulação (um indivíduo por linha);
%       - nvar: número de variáveis.
%
%   Parâmetros de saída:
%       - h: valor da função de otimização modificada
%         para cada indivíduo (vetor coluna);
%       - pop: array de população com valores calculados.

columnsG  = nvar+1 : 2*nvar;
columnsH  = 2*nvar+1 : 3*nvar;
columnF   = 3*nvar+1;

vpop = pop(:, 1:nvar);

fo = f(vpop);
gleo = gle(vpop);
geqo = geq(vpop);

h = fo ...                                              % func. original
  + r * sum((max(zeros(size(vpop)), gleo) .^ 2) ,2) ... % desigualdades
  + s * sum(geqo .^ 2, 2);                              % igualdades

pop(:,columnF) = fo;
pop(:,columnsG) = gleo;
pop(:,columnsH) = geqo;
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
%       - pop: array contendo a população (unm indivíduo por linha).
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
%       - pop: array contendo a população (um indivíduo por linha).
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
%       - pop: array contendo a população (um indivíduo por linha).
%
%   Parâmetros de saída:
%       - g: array contendo a função de restrição (npop linhas,
%         nvar colunas). 

h = cos(2 * pi * pop) + 0.5;

end
