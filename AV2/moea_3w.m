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
%   - Funções objetivo: no colunas.
%   - Valor agregado: uma coluna.
%   - Fronteira de Pareto: uma coluna;
%   - Hyperbox: uma coluna;
%   - Fator de aglomeração (squeeze factor): uma coluna.
xi = 1; xf = xi+nvar-1;
fi = xf+1; ff = fi+no-1;
ag = ff+1;
pt = ag+1;
hy = pt+1;
sq = hy+1;
nc = sq;

L.COLX  = xi:xf;
L.COLF  = fi:ff;
L.COLAG = ag;
L.COLPT = pt;
L.COLHY = hy;
L.COLSQ = sq;
L.NC = nc;

% Dados iniciais da população: faixa das variáveis e 
% quantidade de indivíduos.
xmin = 0; xmax = 1;
npop = 10;

% Estabelecimento da polulação inicial.
pop = popinit(npop,xmin,xmax,L);

% Teste
%pop = dtlz1(pop,nvar,no,L);
pop = dtlz2(pop,nvar,no,L);

resolucao = 5;

pop = agregacao(pop,L);
pop = pareto(pop,L);
pop = hyperbox(pop,resolucao,L);

arqnd = naodominado(pop,L);
arqsq = naoaglomerado(pop,L);
display(arqsq);

%for i = 1:npop
%    display(decodifica(pop(i,L.COLHY),no,resolucao));
%end

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
%     - pop: array contendo a população inicial.

% Aloca o array.
pop = zeros(npop, L.NC);

% Inicializa as variáveis.
pop(:,L.COLX) = xmin * ones(npop,max(L.COLX)) + ...
                (xmax - xmin) * rand(npop,max(L.COLX));
end

function pop = agregacao(pop,L)
%AGREGACAO Calcula o valor agregado das funções objetivo.
%   Calcula o valor agregado das funções objetivo utilizando
%   o método da soma. Considera o mesmo pelo para todas as funções.
%
%   Parâmetros de entrada:
%     - pop: array contendo a população inicial;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - pop: população inicial com valor agregado atualizado.

no = max(L.COLF) - min(L.COLF) + 1; % número de funções objetivo
pop(:,L.COLAG) = sum(pop(:,L.COLF),2) ./ no;

end

function pop = pareto(pop,L)
%PARETO Calcula a fronteira de Pareto.
%   Para cada indivíduo da polulação, calcula a que fronteira de
%   Pareto ele pertence.
%   Inspirado em XXX.
%
%   Parâmetros de entrada:
%     - pop: array contendo a população inicial;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - pop: população inicial com fronteira de Pareto atualizada.

npop = size(pop,1); % número de indivíduos
no = max(L.COLF) - min(L.COLF) + 1; % número de funções objetivo
f = 1;  % f-ésima fronteira de Pareto

% Fronteira: vetor com uma estrutura por fronteira de Pareto:
%   n: número de indivíduos na fronteira;
%   ids: índices dos indivíduos na fronteira.
indfront = struct('n',0, 'ids',zeros(1,npop));
front = repmat(indfront,1,npop);

% Dominação: vetor com uma estrutura por indivíduo i:
%   n: número de indivíduos que dominam i;
%   nd: número de indivíduos dominados por i;
%   d: índices dos indivíduos dominados por i.
elemento = struct('n',0, 'nd',0, 'd',zeros(1,npop));
dominacao = repmat(elemento,1,npop);

% Compara cada indivíduo com todos os demais para determinar
% a dominação.
for i = 1:npop
    for j = 1:npop
        if i ~= j
            % Compara indivíduos i e j
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
    
    % Trata os indivíduos na primeira fronteira.
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
    % Percorre os elementos da f-ésima fronteira, removendo-os
    % e fazendo os cálculos para a (f+1)-ésima fronteira.
    for i = 1:front(f).n
        individuo = front(f).ids(i);
%        display(['Tratando indivíduo ', num2str(individuo)]);
        % Trata os indivíduos dominados pelo indivíduo atual.
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

function pop = hyperbox(pop,resolucao,L)
%HYPERBOX Calcula o hyperbox de cada indivíduo.
%   Para cada indivíduo da população, calcula o hyperbox ao qual
%   pertence e o fator de aglomeração.
%
%   Parâmetros de entrada:
%     - pop: array contendo a população inicial;
%     - resolucao: tamanho do grid (o mesmo para todas as dimensões).
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - pop: população inicial hyperbox e aglomeração atualizados.

% Encontra os valores máximos e mínimos das funções objetivo.
zmax = max(pop(:,L.COLF),[],1);
zmin = min(pop(:,L.COLF),[],1);

% Encontra os limites superiores e inferiores e o tamanho da caixa
delta = (zmax - zmin) / (2 * resolucao);
inf = zmin - delta;
sup = zmax + delta;
boxsize = (sup - inf) / resolucao;

% Calcula o hyperbox (para cada dimensão) e o fator de aglomeracao
npop = size(pop,1); % número de indivíduos
no = max(L.COLF) - min(L.COLF) + 1; % número de funções objetivo
hbox = zeros(npop,no);
indbox = struct('hbox',ones(1,no) * -1, 'n', 0, 'ids',zeros(1,npop));
pbox = repmat(indbox,1,npop);
p = 0;
for i = 1:npop
    hbox(i,:) = floor((pop(i,L.COLF) - inf) ./ boxsize);
    pop(i,L.COLHY) = codifica(hbox(i,:),resolucao);
    encontrado = 0;
    for k = 1:p
        if isequal(hbox(i,:), pbox(k).hbox)
            encontrado = 1;
            pbox(k).n = pbox(k).n + 1;
            pbox(k).ids(pbox(k).n) = i;
            break;
        end
    end
    
    if ~encontrado
        p = p + 1;
        pbox(p).hbox = hbox(i,:);
        pbox(p).n = 1;
        pbox(p).ids(pbox(p).n) = i;
    end
end
%display(hbox);

% Registra o valor de aglomeração para cada indivíduo.
for k = 1:p
%    display(pbox(k));
    for j = 1:pbox(k).n
        pop(pbox(k).ids(j),L.COLSQ) = pbox(k).n;
    end
end

end

function box = codifica(hbox,resolucao)
%CODIFICA Codifica um hyperbox como um valor numérico.
%   Codifica as coordenadas do hyperbox em um valor numérico único.
%
%   Parâmetros de entrada:
%     - hbox: hyperbox a codificar;
%     - resolucao: resolucao utilizada para o cálculo do hyperbox.
%
%   Parâmetros de saída:
%     - box: hyperbox codificado como inteiro.

no = size(hbox,2);
box = 0;
for i = 1:no
    box = (box * resolucao) + hbox(i);
end

end

function hbox = decodifica(box,no,resolucao)
%DECODIFICA Decodifica um hyperbox como um vetor de dimensões.
%   Obtém o vetor de dimensões de um hyperbox a aprtir do valor
%   codificado.
%
%   Parâmetros de entrada:
%     - box: hyperbox a decodificar;
%     - no: número de funções objetivo;
%     - resolucao: resolucao utilizada para o cálculo do hyperbox.
%
%   Parâmetros de saída:
%     - box: hyperbox decodificado.

hbox = zeros(1,no);
for i = 0:no-1
    hbox(no-i) = rem(box,resolucao);
    box = fix(box/resolucao);
end

end

function arqnd = naodominado(pop,L)
%NAODOMINADO Encontra os elementos não dominados de uma populacao.
%   Encontra os elementos não dominados de uma população e retorna-os
%   em uma subpopulação separada.
%
%   Parâmetros de entrada:
%     - pop: array contendo a população inicial;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - arqnd: subpopulação de elementos não dominados.

idx = (pop(:,L.COLPT) == 1);
arqnd = pop(idx,:);

end

function arqsq = naoaglomerado(pop,L)
%NAOAGLOMERADO Encontra os elementos não aglomerados de uma populacao.
%   Encontra os elementos não aglomerados de uma população (ou seja,
%   aqueles sozinhos em um hyperbox) e retorna-os
%   em uma subpopulação separada.
%
%   Parâmetros de entrada:
%     - pop: array contendo a população inicial;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - arqnd: subpopulação de elementos não aglomerados.

idx = (pop(:,L.COLSQ) == 1);
arqsq = pop(idx,:);

end

function pop = dtlz1 (pop,nvar,no,L)
%DTLZ1 Função dtlz1.
%   Calcula as funções objetivo DTLZ1 para todos os indivíduos
%   da população.
%
%   Parâmetros de entrada:
%     - pop: população;
%     - nvar: número de variáveis;
%     - no: número de funções objetivo;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - pop: população com os valores das funções objetivo.

if no > nvar
    ME = MException('dtlz1:nvarError','Invalid number of variables');
    throw(ME);
end

% Colunas das variáveis.
xi = min(L.COLX);
fi = min(L.COLF);
ff = max(L.COLF);
k = nvar - no + 1;

xg = pop(:,no:nvar);

% Cálculo
gx = 100 * (k + sum(((xg-0.5) .^ 2) - cos(20*pi*(xg-0.5)),2));

pop(:,fi) = 0.5 * (prod(pop(:,xi:(no-1)),2) .* (1+gx));

for i = 2:(no-1)
  pop(:,(fi+i-1)) = 0.5 * (prod((pop(:,xi:(no-i))),2) .* ...
                           (1 - pop(:,(no-i+1))) .* (1+gx));
end

pop(:,ff) = 0.5 * ((1 - pop(:,xi)) .* (1+gx));

end

function pop = dtlz2 (pop,nvar,no,L)
%DTLZ1 Função dtlz1.
%   Calcula as funções objetivo DTLZ2 para todos os indivíduos
%   da população.
%
%   Parâmetros de entrada:
%     - pop: população;
%     - nvar: número de variáveis;
%     - no: núemro de funções objetivo;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - pop: população com os valores das funções objetivo.

if no > nvar
    ME = MException('dtlz2:nvarError','Invalid number of variables');
    throw(ME);
end

% Colunas das variáveis.
xi = min(L.COLX);
fi = min(L.COLF);
ff = max(L.COLF);

xg = pop(:,no:nvar);

% Cálculo
gx = sum(((xg-0.5) .^ 2),2);

pop(:,fi) = (1+gx) .* prod(cos(0.5*pi*pop(:,xi:(no-1))),2);

for i = 2:(no-1)
    pop(:,(fi+i-1)) = (1+gx) .* prod(cos(0.5*pi*pop(:,xi:(no-i))),2) .* ...
                                sin(0.5*pi*pop(:,(no-i+1)));
end

pop(:,ff) = (1+gx) .* sin(0.5*pi*pop(:,xi));

end
