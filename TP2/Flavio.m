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
%npop = 50;                % número de indivíduos na população.
npop = 100;                % número de indivíduos na população.
ngen = floor((ncal / npop) - 1); % número de gerações. Calculado a partir
                                 % do número máximo de cálculos da função
                                 % de fitness, considerando npop cálculos 
                                 % por geração, mais npop cálculos para 
                                 % a população inicial.
                          
f   = @rastrigin;         % Handle da função objetivo.
gle = @rastrigin_le;      % Handle da função de restrição de desigualdade.
geq = @rastrigin_eq;      % Handle da função de restrição de igualdade.
xmin = -5.12;             % Valores de contorno para as variáveis
xmax = 5.12;              % de decisão.
                          
% Estabelecimento da população inicial.
pop = popinit(npop,nvar,xmin,xmax);

% Execução do algoritmo genético para otimização.
[x, f, g, h] = ga(pop, ngen, f, gle, geq, xmin, xmax);

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

function [ x, f, g, h ] = ga(pop,ngen,f,gle,geq,xmin,xmax)
%GA Função principal do algoritmo genético.
%   Executa o algoritmo genético (ga) para otimização da função objetivo.
%
%   Parâmetros de entrada:
%       - pop: array contendo a população inicial;
%       - ngen: número de gerações a processar;
%       - f: handle da função de otimização original;
%       - gle: handle da função de restrição de desigualdade;
%       - geq: handle da função de restrição de igualdade;
%       - xmin: valor mínimo de uma variável;
%       - xmax: valor máximo de uma variável.
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
%   - coluna (3*nvar+2): fitness;
%   - coluna (3*nvar+3): número de restrições violadas.
nvar = size(pop,2);
columnsX  = 1:nvar;
columnsG  = nvar+1 : 2*nvar;
columnsH  = 2*nvar+1 : 3*nvar;
columnF   = 3*nvar+1;
%columnFit = 3*nvar+2;
%columnV   = 3*nvar+3;
pais = [pop zeros(size(pop)) zeros(size(pop)) zeros(size(pop,1),3)];

%Parâmetros do algoritmo
ft = @fitness_desc; % handle da função de fitness.
fp = @fitpar_desc;  % handle da função de parâmetros da fitness.
%ft = @fitness_inv; % handle da função de fitness.
%fp = @fitpar_inv;  % handle da função de parâmetros da fitness.

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
pais = fitness(ft, f, gle, geq, r(1), s(1), fitpar(1), pais, nvar);
%display(pais);
pais = popSort(pais,nvar);
%display(pais);

% Inicialização da população de filhos
filhos = zeros(size(pais,1)-nelite, size(pais,2));

% Probabilidades iniciais de cruzamento e mutação.
pc = 0.6;
pm = 0.05;

% Loop de evolução
for g = 1 : ngen
%for g = 1 : 1
    % Guarda os indivíduos da elite
    elite = pais(1 : nelite, :);
    
    % Selecão de pais
    selecionados = selecao(pais,size(filhos,1));
    
    % Calcula as probabilidade de cruzamento e mutação
    [pc, pm] = calcProbs(pais,nvar,pc,pm);
    
    % Efetua os cruzamentos
    lsup = size(filhos,1);
    if rem(lsup,2) == 1; lsup = lsup - 1; end
    
    for fl = 1 : 2 : lsup
        filhos = cruzamento2(filhos,pais, ...
                             selecionados(fl),selecionados(fl+1), ...
                             nvar, fl, pc);
%        filhos = cruzamento(filhos,pais, ...
%                            selecionados(fl),selecionados(fl+1), ...
%                            nvar, fl, pc);
    end
    if rem(size(filhos,1),2) == 1
        filhos(size(filhos,1),:) = pais(selecionados(size(selecionados,2)),:);
    end
    
    % Efetua as mutações
    gamma = perturbacao(filhos,nvar,pm,xmin,xmax,g,ngen);
    filhos(:,columnsX) = filhos(:,columnsX) + gamma;

    % Garante que as condições de contorno sejam respeitadas
    filhos(:,columnsX) = min(xmax, max(xmin,filhos(:,columnsX)));
    
    % Gera a nova população de pais, incluindo a elite
    pais = [elite ; filhos];
    
    % Efetua cálculos e classifica para a nova geração
    %display(pais);
    pais = fitness(ft, f, gle, geq, r(g+1), s(g+1), fitpar(g+1), ...
                   pais, nvar);
    pais = popSort(pais,nvar);
end

x = pais(1,columnsX);
f = pais(1,columnF);
g = pais(1,columnsG);
h = pais(1,columnsH);

end

function [filhos] = cruzamento2(filhos,pais,pai1,pai2,nvar,ind,pc)
%CRUZAMENTO Faz o cruzamento entre dois indivíduos.
%   Realiza o cruzamento estre os dois indivíduos cujos índices são
%   fornecidos e coloca os filhos sequencialmente no vetor de filhos.
%   Se o cruzamento não for realizado, apenas copia os pais para o
%   vetor de filhos.
%   Para as variáveis que participarão do cruzamento, determina
%   individualmente qual é o melhor pai (em relação às restrições), no
%   lugar de utilizar um melhor pai global baseado em um único critério.
%
%   Parâmetros de entrada:
%       - filhos: array dos filhos;
%       - pais: array contendo os pais;
%       - pai1: índice do primeiro pai;
%       - pai2: índice do segundo pai;
%       - nvar: número de variáveis;
%       - ind: índice a primeira linha vaga no array de filhos;
%       - pc: probabilidade de realização do cruzamento.
%
%   Parâmetros de saída:
%       - filhos: array de filhos (preenchido com novos indivíduos).


% Verifica se o cruzamento deve ser realizado; caso contrário, copia.
if rand > pc
    filhos(ind,:) = pais(pai1,:);
    filhos(ind+1,:) = pais(pai2,:);
else
    % Determina o melhor e o pior pai em relação ao valor da fitness.
    if pai1 < pai2
        melhorfit = pai1;
        piorfit = pai2;
    else
        melhorfit = pai2;
        piorfit = pai1;
    end

    % Determina os coeficientes
    kcross = randi(nvar);
    apol = 0.5 * rand + 0.5;
    a = 1.2 * rand - 0.1;
    dir = randi(2) - 1;
    
    if dir == 0
        faixa = 1:kcross;
        filhos(ind,(kcross+1):nvar) = pais(melhorfit,(kcross+1):nvar);
        filhos(ind+1,(kcross+1):nvar) = pais(piorfit,(kcross+1):nvar);
    else
        faixa = kcross:nvar;
        filhos(ind,1:(kcross-1)) = pais(melhorfit,1:(kcross-1));
        filhos(ind+1,1:(kcross-1)) = pais(piorfit,1:(kcross-1));
    end
    
    for i = faixa
        [melhor, pior] = compara(pais,pai1,pai2,i,nvar);
        filhos(ind,i) = apol * pais(melhor,i) + ...
                        (1-apol) * pais(pior,i);
        filhos(ind+1,i) = (1-a) * pais(melhor,i) + a * (pais(pior,i));
    end
    
end
end

function [melhor, pior] = compara(pop,ind1,ind2,var,nvar)
%COMPARA Compara dois indivíduos.
%   Determina qual dentre dois indivíduos é o melhor em relação
%   às restrições de uma variável.
%
%   Parâmetros de entrada:
%       - pop: array de população;
%       - ind1: índice do primeiro indivíduo;
%       - ind2: índice do segundo indivíduo;
%       - var: variável a comparar;
%       - nvar: número de variáveis.
%
%   Parâmetros de saída:
%       - melhor: índice do melhor indivíduo;
%       - pior: índice do pior indivíduo.

colg = nvar + var;     % Colunas das violações da
colh = 2*nvar + var;   % variável em questão.
     
g1 = pop(ind1,colg) <= 0;
g2 = pop(ind2,colg) <= 0;
h1 = pop(ind1,colh);
h2 = pop(ind2,colh);
     
if g1 == g2
    % Ambos respeitam ou desrespeitam a restrição de desigualdade.
    % Seleciona pela restrição de igualdade.
    if abs(h1) < abs(h2)
        melhor = ind1;
        pior = ind2;
    else
        melhor = ind2;
        pior = ind1;
    end
elseif g1 == 1
    % Indivíduo 1 na faixa, 2 fora.
    melhor = ind1;
    pior = ind2;
else
    % Indivíduo 2 na faixa, 1 fora.
    melhor = ind2;
    pior = ind1;
end
     
end

function [gamma] = perturbacao(pop,nvar,pm,xmin,xmax,g, ngen)
%PERTURBACAO Calcula a matriz de perturbações.
%   A matriz de perturbações corresponde aos valores a serem
%   adicionados às variáveis dos indivíduos para a implementação
%   da mutação.
%
%   Parâmetros de entrada:
%       - pop: array de indivíduos;
%       - nvar: número de variáveis;
%       - pm: probabilidade de mutação;
%       - xmin: valor mínimo de uma variável;
%       - xmax: valor máximo de uma variável;
%       - g: número da geração;
%       - ngen: total de gerações.
%
%   Parâmetros de saída:
%       - gamma; matriz de perturbações.

npop = size(pop,1);
gamma = zeros(npop,nvar);
maxg = floor(0.8 * ngen);  % geração máxima para cálculo por range

for f = 1 : npop
    if rand < pm
        % Parâmetros.
        kmut = randi(nvar);
        range = xmax - xmin;
        dir = randi(2) - 1;
             
        % Colunas a mutar.
        if dir == 0
            tomut = 1 : kmut;
        else
            tomut = kmut : nvar;
        end
        
        for v = tomut
            if g <= maxg
                val = range;
            else
                val = sum(pop(:,v)) / npop;
            end
            gamma(f,v) = 0.05 * (2 * rand - 1) * val;
        end
    end
end

end

function [pc, pm] = calcProbs(pop,nvar,pc,pm)
%CALCPROBS Calcula as probabilidades de cruzamento e mutação.
%   Calcula as probabilidades de cruzamento e mutação da presente geração.
%   Essas probabilidades são influenciadas pelo mdg.
%
%   Parâmetros de entrada:
%       - pop: array com os indivíduos da população;
%       - nvar: número de variáveis;
%       -8 pc: probabilidade de cruzamento atual;
%       - pm: probabilidade de mutação atual.
%
%   Parâmetros de saída:
%       - pc: nova probabilidade de cruzamento;
%       - pm: nova probabilidade de mutação.


vinf = 0.3;
vmax = 0.7;
kc = 1.2;
km = 1.2;

mdg = calcMdg(pop,nvar);

if mdg < vinf
  pc = pc / kc;
  pm = pm * km;
elseif mdg > vmax
  pc = pc * kc;
  pm = pm / km;
end

end

function [mdg] = calcMdg(pop, nvar)
%CALCMDG Calcula a medida da diversidade genética.
%   A medida da diversidade genética da população influencia nas
%   probabilidades de cruzamento e mutação.
%
%   Parâmetros de entrada:
%       - pop: array com os indivíduos da população;
%       - nvar: número de variáveis.
%
%   Parâmetros de saída:
%       - mdg: medida de diversidade genética da população.

columnFit = 3*nvar+2;

fit = pop(:,columnFit);
fmed = sum(fit) / size(pop,1);
fmax = max(fit);

mdg = (fmax - fmed) / fmax;
end

function [filhos] = cruzamento(filhos,pais,pai1,pai2,nvar,ind,pc)
%CRUZAMENTO Faz o cruzamento entre dois indivíduos.
%   Realiza o cruzamento estre os dois indivíduos cujos índices são
%   fornecidos e coloca os filhos sequencialmente no vetor de filhos.
%   Se o cruzamento não for realizado, apenas copia os pais para o
%   vetor de filhos.
%
%   Parâmetros de entrada:
%       - filhos: array dos filhos;
%       - pais: array contendo os pais;
%       - pai1: índice do primeiro pai;
%       - pai2: índice do segundo pai;
%       - nvar: número de variáveis;
%       - ind: índice a primeira linha vaga no array de filhos;
%       - pc: probabilidade de realização do cruzamento.
%
%   Parâmetros de saída:
%       - filhos: array de filhos (preenchido com novos indivíduos).


% Verifica se o cruzamento deve ser realizado; caso contrário, copia.
if rand > pc
    filhos(ind,:) = pais(pai1,:);
    filhos(ind+1,:) = pais(pai2,:);
else
    % Determina o melhor e o pior pai.
    if pai1 < pai2
        melhor = pai1;
        pior = pai2;
    else
        melhor = pai2;
        pior = pai1;
    end

    % Determina os coeficientes
    kcross = randi(nvar);
    apol = 0.5 * rand + 0.5;
    a = 1.2 * rand - 0.1;
    dir = randi(2) - 1;
    
    if dir == 0
        % Cruzamento para a esquerda
        filhos(ind,1:kcross) = apol * pais(melhor,1:kcross) + ...
                               (1-apol) * pais(pior,1:kcross);
        filhos(ind,(kcross+1):nvar) = pais(melhor,(kcross+1):nvar);
        filhos(ind+1,1:kcross) = (1-a)*pais(melhor,1:kcross) + ...
                                 a*pais(pior,1:kcross);
        filhos(ind+1,(kcross+1):nvar) = pais(pior,(kcross+1):nvar);
    else
        % Cruzamento para a direita
        filhos(ind,1:(kcross-1)) = pais(melhor,1:(kcross-1));
        filhos(ind,kcross:nvar) = apol*pais(melhor,kcross:nvar) + ...
                                  (1-apol) * pais(pior,kcross:nvar);
        filhos(ind+1,1:(kcross-1)) = pais(pior,1:(kcross-1));
        filhos(ind+1,kcross:nvar) = (1-a)*pais(melhor,kcross:nvar) + ...
                                    a*pais(pior,kcross:nvar); 
    end
end
end

function [selecionados] = selecao(pais,nfilhos)
%SELECTION Seleciona os pais que participarão dos próximos processos.
%   Utiliza torneio binário.
%
%   Parâmetros de entrada:
%       - pais: array contendo os pais;
%       - nfilhos: número de filhos.
%
%   Parâmetros de saída:
%       - filhos: vetor com índices dos pais selecionados.

npais = size(pais,1);   % quantidade de pais.
ultimo = -1;            % índice do último pai selecionado.
selecionados = zeros(1,nfilhos);

f = 1;
while f <= nfilhos
    % Seleciona dois pais para o torneio.
    p1 = randi(npais);
    p2 = randi(npais);
    while p1 == p2
        p2 = randi(npais);
    end
    
    % Determina o melhor e o pior
    if p1 < p2
        melhor = p1;
        pior = p2;
    else
        melhor = p2;
        pior = p1;
    end
    
    % Torneio
    if rand < 0.75
        selecionado = melhor;
    else
        selecionado = pior;
    end
    
    % Evita o cruzamento de un indivíduo consigo mesmo
    if rem(f,2) == 0 && ultimo == selecionado
        continue;
    end
    
    % Registra o pai selecionado
    ultimo = selecionado;
    selecionados(f) = selecionado;
    f = f + 1;
end
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
columnV = 3*nvar+3;   % Coluna da quantidade de violações de restrição.

%Ordena pela função de fitness apenas (decrescente).
pop = sortrows(pop, -columnFit);

%Ordena pela quantidade de violações e pela fitness (decrescente).
%pop = sortrows(pop,[columnV -columnFit]);
end

function [v] = violacoes(pop, nvar)
%VIOLACOES Calcula o número restrições violadas por um indivíduo.
%   Retorna, para cada indivíduo, o número de restrições (igualdade
%   e desigualdade) violadas por suas variáveis).
%
%   Parâmetros de entrada:
%       - pop: array de população.
%       - nvar: número de variáveis.
%
%   Parâmetros de saída:
%       - v: número de restrições violadas por indivíduo(vetor coluna).

columnsG  = nvar+1 : 2*nvar;
columnsH  = 2*nvar+1 : 3*nvar;

% Variação admissível no teste de violação da restrição de igualdade
% (testes de igualdade com valores em ponto flutuante não são precisos).
delta = 0.0000001;
%delta = 0.1;

%nv = sum((ceil(10 .^ (max(0,pop(:,columnsG))))-1),2) + sum((ceil(10 .^ abs(pop(:,columnsH)))-1),2);
%display([pop(:,columnsG) pop(:,columnsH)]);
%display(nv)

%v = sum((pop(:,columnsG) > 0),2) + sum((pop(:,columnsH) ~= 0),2);
v = sum((pop(:,columnsG) > 0),2) + sum(abs(pop(:,columnsH)) > delta,2);
%v = nv;
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

columnFit = 3*nvar+2;
columnV = 3*nvar+3;

pop = ft(f, gle, geq, r, s, n, pop, nvar);
pop(:,columnFit) = escalonamento(pop,nvar);
pop(:,columnV) = violacoes(pop,nvar);
end

function [e] = escalonamento(pop,nvar)
%ESCALONAMENTO Faz o escalonamento da fitness dos indivíduos.
%   Escalona a fitness dos indivíduos para evitar que um "superindivíduo"
%   domine a população. Utiliza escalonamento linear.
%
%   Parâmetros de entrada:
%       - pop: array de população.
%       - nvar: número de variáveis.
%
%   Parâmetros de saída:
%       - e: fitness escalonada (vetor coluna).

columnFit = 3*nvar+2; % Coluna da fitness
Cmax = 2; % número máximo de cópias do melhor indivíduo.

fit = pop(:,columnFit);
fmed = sum(fit) / size(pop,1);
fmax = max(fit);
a = (Cmax-1)*fmed / (fmax-fmed);
b = (1-a) / fmed;
e = max((a * fit + b), zeros(size(fit)));
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
