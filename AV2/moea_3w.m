function [ ps ] = moea_3w( ncal, npop, nvar, no, precisao, objfunc )
%MOEA_3W Multi Objective Evolutionary Algorithm - 3 Way
%   Esta função implementa um algoritmo para otimização de funções
%   Multi-objetivo. XXX Explicar o algoritmo.
%
%   Parâmetros de entrada:
%     - ncal: número máximo de cálculos das funções objetivo;
%     - npop: tamanho da população;
%     - nvar: número de variáveis;
%     - no: número de funções objetivo;
%     - precisao: XXX
%     - objfunc: handle da função objetivo
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
%   - Fator de aglomeração (squeeze factor): uma coluna;
%   - Ponto de referência: uma coluna;
%   - Fator de nicho: uma coluna;
%   - Distância ao ponto de referência: uma coluna.
xi = 1; xf = xi+nvar-1;
fi = xf+1; ff = fi+no-1;
ag = ff+1;
pt = ag+1;
hy = pt+1;
sq = hy+1;
rf = sq+1;
nh = rf+1;
dt = nh+1;
nc = dt;

L.COLX  = xi:xf;
L.COLF  = fi:ff;
L.COLAG = ag;
L.COLPT = pt;
L.COLHY = hy;
L.COLSQ = sq;
L.COLRF = rf;
L.COLNH = nh;
L.COLDT = dt;
L.NC = nc;

% Dados iniciais da população: faixa das variáveis e 
% quantidade de indivíduos.
xmin = 0; xmax = 1;

% Estabelecimento da população inicial.
pop = popinit(npop,xmin,xmax,L);

% Avaliação dos elementos da população inicial.
pop = calcFos(pop,objfunc,L);
ncal = ncal - npop;

pop = agregacao(pop,L); % cálculo do valor agregado
pop = pareto(pop,L);    % cálculo da fronteira de pareto
[inf,sup,boxsize] = calParEspaco(pop,precisao,L); % XXX
[pop, pbox] = hyperbox(pop,inf,boxsize,L); % cálculo do hyperbox
[pop, pref] = calcRefPontos(inf,sup,npop,no,pop,L); % cálculo dos pontos de referência

% Cópia dos elementos não dominados da população para o arquivo arqnd
arqnd = naodominado(pop,L);
[arqnd,ndbox] = hyperbox(arqnd,inf,boxsize,L);
[arqnd,ndref] = calcRefPontos(inf,sup,npop,no,arqnd,L); % cálculo dos pontos de referência

% Cópia dos elementos não agregados da população para o arquivo arqsq
arqsq = selecionaDistribuicao(pop,L);
[arqsq,sqbox] = hyperbox(arqsq,inf,boxsize,L);
[arqsq,sqref] = calcRefPontos(inf,sup,npop,no,arqsq,L); % cálculo dos pontos de referência

% Probabilidades iniciais de cruzamento e mutação
pc = 0.6;
pm = 0.05;

% Laço principal
while ncal >= 3
    % Seleciona os pais para os cruzamentos (um de cada populaçao/arquivo)
    p1 = seleciona(pop,@compara,L);
    p2 = seleciona(arqnd,@compara,L);
    p3 = seleciona(arqsq,@compdist,L);
    
    % Calcula as probabilidades de cruzamento e mutação
    [pc,pm] = calcProbs(pop,pc,pm,L);

    % Realiza os cruzamentos para gerar três filhos.
    c1 = cruzamento(p1,p2,pc,L);
    c2 = cruzamento(p1,p3,pc,L);
    c3 = cruzamento(p2,p3,pc,L);
    filhos = [c1;c2;c3];
    
    % Tratamentos dos filhos
    for i = 1:size(filhos,1)
        
        c = filhos(i,:);
        
        % Realiza mutações
        gamma = perturbacao(c,pm,xmin,xmax,L);
        c(L.COLX) = c(L.COLX) + gamma;

        % Garante que as condições de contorno sejam respeitadas
        c(:,L.COLX) = min(xmax, max(xmin,c(:,L.COLX)));
        
        % Calcula as funções objetivo e o valor agregado
        c = calcFos(c,objfunc,L);
        c = agregacao(c,L);
        
        % Tenta inserir c na população e nos arquivos
        [pop,pbox,pref,pextrapola] = atualizapop(c,pop,inf,sup,pbox,pref,boxsize,npop,L);

        [arqnd,ndbox,ndref,ndextrapola] = atualizadom(c,arqnd,inf,sup,ndbox,ndref,...
                                          npop,boxsize,npop,L);
        [arqsq,sqbox,sqref,sqextrap] = atualizadis(c,arqsq,inf,sup,sqbox,sqref,...
                                          npop,boxsize,npop,L);
        
        % Algum novo ponto extrapolou o espaço; recalcula
        if pextrapola || ndextrapola || sqextrap
            [inf,sup,boxsize] = calParEspaco(pop,precisao,L);
            [pop, pbox] = hyperbox(pop,inf,boxsize,L);
            [pop, pref] = calcRefPontos(inf,sup,npop,no,pop,L);
            pop = pareto(pop,L);
            [arqnd,ndbox] = hyperbox(arqnd,inf,boxsize,L);
            [arqnd,ndref] = calcRefPontos(inf,sup,npop,no,arqnd,L);
            arqnd = pareto(arqnd,L);
            [arqsq,sqbox] = hyperbox(arqsq,inf,boxsize,L);
            [arqsq,sqref] = calcRefPontos(inf,sup,npop,no,arqsq,L);
            arqsq = selecionaDistribuicao(arqsq,L);
            arqsq = pareto(arqsq,L);
        end
        
    end
  
    % Foram avaliadas as funções objetivo de três indivíduos
    ncal = ncal - 3;
    display(ncal);
end

% Retorno do resultado (elementos na primeira fronteira de Pareto)
ps = arqnd;
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
%   Para cada indivíduo da população, calcula a que fronteira de
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

% Tratamento das demais fronteiras
while front(f).n ~= 0
    % Percorre os elementos da f-ésima fronteira, removendo-os
    % e fazendo os cálculos para a (f+1)-ésima fronteira.
    for i = 1:front(f).n
        individuo = front(f).ids(i);
        % Trata os indivíduos dominados pelo indivíduo atual.
        for j = 1:dominacao(individuo).nd
            dominado = dominacao(individuo).d(j);
            dominacao(dominado).n = dominacao(dominado).n - 1;
            if dominacao(dominado).n == 0
                front(f+1).n = front(f+1).n + 1;
                front(f+1).ids(front(f+1).n) = dominado;
                pop(dominado,L.COLPT) = f+1;
            end
        end
    end
    
    f = f + 1;
end

end

function [inf, sup, boxsize] = calParEspaco(pop,precisao,L)
%CALCPARESPACO Calcula parâmetros do espaço de pontos.
%   XXX
%
%   Parâmetros de entrada:
%     - pop: array contendo a população;
%     - precisao: XXX
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - inf: limites inferiores do grid para cada dimensão;
%     - sup: limites superiores do grid para cada dimensão;
%     - boxsize: XXX

% Encontra os valores máximos e mínimos das funções objetivo.
zmax = max(pop(:,L.COLF),[],1);
zmin = min(pop(:,L.COLF),[],1);

resolucao = calcResolucao(pop,precisao,L);

% Encontra os limites superiores e inferiores e o tamanho da caixa
delta = (zmax - zmin) ./ (2 * resolucao);
inf = zmin - delta;
sup = zmax + delta;
boxsize = (sup - inf) ./ resolucao;

end

function [pop, pbox] = hyperbox(pop,inf,boxsize,L)
%HYPERBOX Calcula o hyperbox de cada indivíduo.
%   Para cada indivíduo da população, calcula o hyperbox ao qual
%   pertence e o fator de aglomeração.
%
%   Parâmetros de entrada:
%     - pop: array contendo a população;
%     - inf: XXX;
%     - boxsize: XXX;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - pop: população atualizados;
%     - pbox: vetor de estruturas com as informações dos boxes:
%             hbox: coordenadas do box;
%             n: número de elementos do box;
%             ids: índices dos indivíduos do box.

% Calcula o hyperbox (para cada dimensão) e o fator de aglomeracao
npop = size(pop,1); % número de indivíduos
no = max(L.COLF) - min(L.COLF) + 1; % número de funções objetivo
hbox = zeros(npop,no);
indbox = struct('hbox',ones(1,no) * -1, 'n', 0, 'ids',zeros(1,npop));
pbox = repmat(indbox,1,npop);
p = 0;
for i = 1:npop
    hbox(i,:) = floor((pop(i,L.COLF) - inf) ./ boxsize);
    encontrado = 0;
    for k = 1:p
        if isequal(hbox(i,:), pbox(k).hbox)
            encontrado = 1;
            pop(i,L.COLHY) = k;
            pbox(k).n = pbox(k).n + 1;
            pbox(k).ids(pbox(k).n) = i;
            break;
        end
    end
    
    if ~encontrado
        p = p + 1;
        pop(i,L.COLHY) = p;
        pbox(p).hbox = hbox(i,:);
        pbox(p).n = 1;
        pbox(p).ids(pbox(p).n) = i;
    end
end

% Registra o valor de aglomeração para cada indivíduo.
for k = 1:p
    %display(pbox(k));
    for j = 1:pbox(k).n
        pop(pbox(k).ids(j),L.COLSQ) = pbox(k).n;
    end
end

end

function [pop,refpontos] = calcRefPontos(inf,sup,npontos,no,pop,L)
%REFPONTOS XXX
%   XXX
%
%   Parâmetros de entrada:
%     - inf: XXX
%     - sup: XXX
%     - npontos: XXX
%     - no: XXX
%     - pop: XXX
%     - L: XXX
%
%   Parâmetros de saída:
%     - pop: XXX
%     - refpontos: XXX

% Observação: o cálculo do número de subdivisões para cada
% dimensão abaixo pode gerar um número muito grande de pontos
% de referência para problemas com muitos objetivos. Deve ser
% estudada uma forma de estabelecer um limite superior para esse
% número para não comprometer a eficiência do algoritmo.
n = max(2,round(nthroot(npontos,no))); % nro de subdivisoes por dimensão

npop = size(pop,1);

delta = (sup - inf) / n; % espaçamento dos pontos (por dimensão)
npontos = (n+1) ^ no; % nro total de pontos de referência

strucponto = struct ('ponto',zeros(1,no), 'n', 0, 'ids',zeros(1,npop));
refpontos = repmat(strucponto,1,npontos);

% Preenche os pontos de referência
ponto = inf; % Começa pelo "canto inferior esquerdo"
for k = 1:npontos
    refpontos(k).ponto = ponto; % grava o ponto
    
    % Faz incremento em uma das dimensões para o próximo ponto
    tol = 1e-6;
    d = 1; % dimensão a incrementar
    incrementado = 0; % incremento bem sucedido?
    while(~incrementado)
        ponto(d) = ponto(d) + delta(d); % incrementa
        if ponto(d) <= sup(d) + tol % se superou máximo, muda dimensão
            incrementado = 1;
        else
            ponto(d) = inf(d);
            d = d+1;
        end
        
        if d > no  % segurança
            break;
        end
    end
end

% Encontra o ponto de referência mais próximo de cada indivíduo
for i = 1:npop
    prox = 1;
    distprox = norm(pop(i,L.COLF) - refpontos(1).ponto);
    
    for j = 2:npontos
        dist = norm(pop(i,L.COLF) - refpontos(j).ponto);
        if dist < distprox
            prox = j;
            distprox = dist;
        end
    end
    
    % XXX Gravar aqui informações no indivíduo
    pop(i,L.COLDT) = distprox;
    refpontos(prox).n = refpontos(prox).n + 1;
    refpontos(prox).ids(refpontos(prox).n) = i;
end

% Grava as informações nos indivíduos
for k = 1:npontos
    for i = 1:refpontos(k).n
        pop(refpontos(k).ids(i), L.COLRF) = k;
        pop(refpontos(k).ids(i), L.COLNH) = refpontos(k).n;        
    end
end

end

function p = seleciona(pop,comp,L)
%SELECIONA Seleciona um indivíduo da população.
%   Seleciona um indivíduo para participar dos cruzamentos.
%   Utiliza torneio binário.
%
%   Parâmetros de entrada:
%     - pop: população;
%     - comp: handle da função de comparação;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - p: indivíduo selecionado.

npop = size(pop,1); % número de indivíduos

% Seleciona dois pais para o torneio
p1 = randi(npop);
p2 = randi(npop);
while p1 == p2
    p2 = randi(npop);
end

% Compara os dois pais
[melhor, pior] = comp(p1, p2, pop, L);

% Torneio
prob = 0.75;
if rand <= prob
    p = pop(melhor,:);
else
    p = pop(pior,:);
end

end

function [melhor, pior] = compara(i1, i2, pop, L)
%COMPARA Compara dois indivíduos e retorna-os ordenados.
%   Utiliza critérios normais de comparação:
%   pareto > fator de nicho> fator de agregação > valor agregado.
%
%   Parâmetros de entrada:
%     - i1: índice do indivíduo 1;
%     - i2: índice do indivíduo 2;
%     - pop: população;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - melhor: índice do melhor indivíduo;
%     - pior: índice do pior indivíduo.

ind1 = pop(i1,:);
ind2 = pop(i2,:);

if ind1(L.COLPT) < ind2(L.COLPT)
    melhor = i1;
    pior = i2;
elseif ind1(L.COLPT) > ind2(L.COLPT)
    melhor = i2;
    pior = i1;
elseif ind1(L.COLNH) < ind2(L.COLNH)
    melhor = i1;
    pior = i2;
elseif ind2(L.COLNH) < ind1(L.COLNH)
    melhor = i2;
    pior = i1;
elseif ind1(L.COLSQ) < ind2(L.COLSQ)
    melhor = i1;
    pior = i2;
elseif ind1(L.COLSQ) > ind2(L.COLSQ)
    melhor = i2;
    pior = i1;
elseif ind1(L.COLAG) < ind2(L.COLAG)
    melhor = i1;
    pior = i2;
else
    melhor = i2;
    pior = i1;
end

end

function [melhor, pior] = compdist(i1, i2, pop, L)
%COMPDIST Compara dois indivíduos e retorna-os ordenados.
%   Utiliza critério alternativo (prioriza distribuição):
%   fator de nicho > fator de agregação > pareto > valor agregado.
%
%   Parâmetros de entrada:
%     - i1: índice do indivíduo 1;
%     - i2: índice do indivúduo 2;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - melhor: melhor indivíduo;
%     - pior: pior indivíduo.

ind1 = pop(i1,:);
ind2 = pop(i2,:);

if ind1(L.COLNH) < ind2(L.COLNH)
    melhor = i1;
    pior = i2;
elseif ind2(L.COLNH) < ind1(L.COLNH)
    melhor = i2;
    pior = i1;
elseif ind1(L.COLSQ) < ind2(L.COLSQ)
    melhor = i1;
    pior = i2;
elseif ind1(L.COLSQ) > ind2(L.COLSQ)
    melhor = i2;
    pior = i1;
elseif ind1(L.COLPT) < ind2(L.COLPT)
    melhor = i1;
    pior = i2;
elseif ind1(L.COLPT) > ind2(L.COLPT)
    melhor = i2;
    pior = i1;
elseif ind1(L.COLAG) < ind2(L.COLAG)
    melhor = i1;
    pior = i2;
else
    melhor = i2;
    pior = i1;
end

end

function [melhor, pior] = comprefs(i1, i2, pop, L)
%COMPDIST Compara dois indivíduos e retorna-os ordenados.
%   Utiliza distância a ponto de referência para comparação.
%
%   Parâmetros de entrada:
%     - i1: índice do indivíduo 1;
%     - i2: índice do indivúduo 2;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - melhor: melhor indivíduo;
%     - pior: pior indivíduo.

ind1 = pop(i1,:);
ind2 = pop(i2,:);

if ind1(L.COLDT) < ind2(L.COLDT)
    melhor = i1;
    pior = i2;
else
    melhor = i2;
    pior = i1;
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

function arqsq = selecionaDistribuicao(pop,L)
%SelecionaDistribuicao Seleciona os elementos do arquivo de distribuição.
%   Examina os elementos de uma população e seleciona os melhores de
%   cada nicho para compor o arquivo de distribuição.
%
%   Parâmetros de entrada:
%     - pop: array contendo a população inicial;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - arqsq: subpopulação de elementos não aglomerados.

popsort = sortrows(pop,[L.COLRF L.COLDT L.COLHY L.COLPT L.COLAG]);

numnichos = length(unique(popsort(:,L.COLRF)));
arqsq = zeros(numnichos,L.NC);

lastnicho = -1;
ind = 1;

for i = 1:size(popsort,1)
    if popsort(i,L.COLRF) ~= lastnicho
        arqsq(ind,:) = popsort(i,:);
        arqsq(ind,L.COLNH) = 1;
        ind = ind + 1;
        lastnicho = popsort(i,L.COLRF);
    end
end

end

function [pop,pbox,pref,extrap] = atualizapop(c,pop,inf,sup,pbox,pref,boxsize,npontos,L)
%ATUALIZAPOP Atualiza a população com um indivíduo.
%   Tenta inserir um indivíduo em uma população. XXX
%
%   Parâmetros de entrada:
%     - c: indivíduo a inserir;
%     - pop: população;
%     - inf: limites inferiores do grid para cada dimensão;
%     - sup: limites superiores do grid para cada dimensão;
%     - pbox: box da população;
%     - pref: XXX
%     - boxsize: XXX
%     - npontos: XXX
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - pop: população atualizada.
%     - pbox: XXX
%     - pbof: XXX
%     - extrapola: XXX

% Verifica se indivíduo não existe na população
%idx = indice(c,pop,L);
%if idx > 0
%    display('c já existe: descartado');
%    return;
%end

% Calcula os indivíduos de pop dominados e que dominam c.
dom = dominacao(c,pop,L);

extrap = extrapola(c,inf,sup,L);
descartado = 0;

if dom.n > 0
    % c domina um grupo de indivíduos:
    % substitui o pior deles por c.
    inds = dom.idn(1:dom.n);
    pr = pior(pop,inds,@compara,L);
    % XXX atualizar pbox e pref
    [pop,pbox] = removebox(pop,pr,pbox,L);
    [pop,pref] = removeref(pop,pr,pref,L);
    pop(pr,:) = c;
elseif dom.m > 0
    % c é dominado e descartado.
    descartado = 1;
else
    % c não domina nem é dominado:
    % substitui o pior indivíduo da população.
    inds = 1:size(pop,1);
    pr = pior(pop,inds,@compara,L);
    % XXX atualizar pbox e pref
    [pop,pbox] = removebox(pop,pr,pbox,L);
    [pop,pref] = removeref(pop,pr,pref,L);
    pop(pr,:) = c;
    % XXX atualiza pbox e pref
end

if ~descartado
    % c foi inserido na população: refaz os cálculos de dominação
    % e distribuição.
    pop = pareto(pop,L);
    if ~extrap
        ibox = indbox(pop(pr,:),inf,pbox,boxsize,L);
        if ibox == 0
            % Caiu em novo box; recalcula
            [pop, pbox] = hyperbox(pop,inf,boxsize,L);
        else
            % Caiu em box existente; atualiza
            [pop,pbox] = inserebox(pop,pr,inf,pbox,boxsize,L);
            pop(pr,L.COLHY) = ibox;
            pop(pr,L.COLSQ) = pbox(ibox).n;
        end
        
        [iref,dist] = indref(pop(pr,:),pref,L);
        if iref == 0
            % Caiu em novo ponto de referência; recalcula
            no = max(L.COLF) - min(L.COLF) + 1;
            [pop, pref] = calcRefPontos(inf,sup,npontos,no,pop,L);
        else
            % Caiu em ponto existente; atualiza
            [pop,pref] = insereref(pop,pr,pref,L);
            pop(pr,L.COLRF) = iref;
            pop(pr,L.COLNH) = pref(iref).n;
            pop(pr,L.COLDT) = dist;
        end
    end
else
    extrap = 0; % Não considera extrapolação para ponto descartado
end

end

function [arqnd,ndbox,ndref,extrap] = atualizadom(c,arqnd,inf,sup,ndbox,ndref, ...
                                           maxpontos,boxsize,npontos,L)
%ATUALIZAPOP Atualiza o arquivo de dominação com um indivíduo.
%   Tenta inserir um indivíduo no arquivo de dominação população.
%   Se o indivíduo não for descartado, refaz os cálculos 
%   das fronteiras de Pareto e do hyperbox.
%
%   Parâmetros de entrada:
%     - c: indivíduo a inserir;
%     - arqnd: população (arquivo de dominação);
%     - inf: limites inferiores do grid para cada dimensão;
%     - sup: limites superiores do grid para cada dimensão;
%     - ndbox: boxes do arquivo;
%     - ndref: XXX
%     - max: número máximo de elementos em no arquivo;
%     - boxsize: XXX
%     - npontos: XXX
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - arqnd: população atualizada.
%     - ndbox: XXX
%     - ndref: XXX
%     - extrap: XXX

extrap = 0;

% Verifica se indivíduo não existe na população
idx = indice(c,arqnd,L);
if idx > 0
    %display('c já existe: descartado');
    return;
end

% Calcula os indivíduos de pop dominados e que dominam c.
dom = dominacao(c,arqnd,L);

extrap = extrapola(c,inf,sup,L);
descartado = 0;

if dom.m > 0
    % c é dominado por algum indivíduo do arquivo e é descartado
    descartado = 1;
elseif dom.n > 0
    % c domina indivíduos do arquivo: remove os dominados e insere c
    arqnd = remove(arqnd,dom.idn,dom.n,L);
    arqnd = [arqnd; c];
elseif size(arqnd,1) < maxpontos
    % c não domina nem é dominado e há espaço no arquivo: insere c
    arqnd = [arqnd; c];
else
    % c não domina nem é dominado; não há espaço no arquivo
    % Verifica se c estende o grid do arquivo
    if extrap
        % Substitui o pior elemento do grid por c
        pr = pior(arqnd,1:size(arqnd,1),@compara,L);
        arqnd(pr,:) = c;
    else
        % Encontra o nicho e a distância do elemento
        [iref,dist] = indref(c,ndref,L);
        c(L.COLPT) = 1;
        c(L.COLRF) = iref;
        c(L.COLDT) = dist;

        % Procura por elementos do mesmo nicho
        rfinds = find(arqnd(:,L.COLRF) == iref);
        c(L.COLNH) = length(rfinds) + 1;

        % Se houver, encontra o pior e tenta substituir
        inserido = 0;
        if ~isempty(rfinds)
            rfinds = transpose(rfinds);
            pr = pior(arqnd,rfinds,@comprefs,L);
            cpr = arqnd(pr,:);
            % Compara o pior com c
            [m,~] = comprefs(1,2,[c;cpr],L);
            if m == 1
                % c é melhor: substitui
                idx = indice(cpr,arqnd,L);
                arqnd(idx,:) = c;
                inserido = 1;
            end
        end
        
        % Se não foi substituído, tenta contra o pior do arquivo
        if ~inserido
            % Encontra o box do indivíduo
            hy = indbox(c,inf,ndbox,boxsize,L);
            % Indivíduos no mesmo hyperbox (índices)
            hyinds = find(arqnd(:,L.COLHY) == hy);
            sq = size(hyinds,1) + 1;
            c(L.COLHY) = hy;
            c(L.COLSQ) = sq;
            % Encontra o pior da população
            pr = pior(arqnd,1:size(arqnd,1),@compara,L);
            cpr = arqnd(pr,:);
            if cpr(L.COLHY) == hy
                cpr(L.COLSQ) = cpr(L.COLSQ) + 1;
            end
            % Compara com o indivíduo
            [m,~] = compara(1,2,[c;cpr],L);
            if m == 1
                % c é melhor: substitui
                arqnd(pr,:) = c;
            else
                % c é pior: descarta
                descartado = 1;
            end
        end
        
    end
end

if ~descartado
    % c foi inserido na população: refaz os cálculos de dominação
    % e distribuição.
    if ~extrap
        % Observação: este código pode ser otimizado para realizar as
        % atualizações estritamente necessárias, tal como foi feito
        % no caso da população geral.
        no = max(L.COLF) - min(L.COLF) + 1;
        arqnd = pareto(arqnd,L);
        [arqnd,ndbox] = hyperbox(arqnd,inf,boxsize,L);
        [arqnd,ndref] = calcRefPontos(inf,sup,npontos,no,arqnd,L);
    end
else
    extrap = 0; % Não considera extrapolação para ponto descartado
end
    
end

function [arqsq,sqbox,sqref,extrap] = atualizadis(c,arqsq,inf,sup,sqbox,sqref, ...
                                           maxpontos,boxsize,npontos,L)
%ATUALIZADIS Atualiza o arquivo de dominação com um indivíduo.
%   Tenta inserir um indivíduo no arquivo de distribuição.
%   Se o indivíduo não for descartado, refaz os cálculos 
%   das fronteiras de Pareto e do hyperbox.
%
%   Parâmetros de entrada:
%     - c: indivíduo a inserir;
%     - arqsq: população (arquivo de distribuição);
%     - inf: limites inferiores do grid para cada dimensão;
%     - sup: limites superiores do grid para cada dimensão;
%     - sqbox: boxes do arquivo;
%     - sqref: XXX
%     - maxpontos: número máximo de elementos no arquivo;
%     - boxsize: XXX
%     - npontos: número de pontos de referência;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - arqsq: população atualizada.
%     - sqbox: boxes atualizados;
%     - sqref: XXX
%     - extrap: XXX

extrap = 0;

% Verifica se indivíduo não existe na população
idx = indice(c,arqsq,L);
if idx > 0
    %display('c já existe: descartado');
    return;
end

extrap = extrapola(c,inf,sup,L);
descartado = 0;

if extrap
    % c estende grid do arquivo: deve ser inserido.
    if size(arqsq,1) < maxpontos
        % Arquivo tem espaço: insere c.
        arqsq = [arqsq; c];
    else
        % Arquivo está cheio: substitui o pior.
        pr = pior(arqsq,1:size(arqsq,1),@compdist,L);
        arqsq(pr,:) = c;
    end
else
    % C não estende grid: verifica de deve ser inserido
    % Calcula o nicho e a distância de c.
    [iref,dist] = indref(c,sqref,L);
    c(L.COLRF) = iref;
    c(L.COLDT) = dist;
    
    % Procura por elementos do mesmo nicho
    rfinds = find(arqsq(:,L.COLRF) == iref);
    nh = length(rfinds) + 1;
    c(L.COLNH) = 1; % Se inserido, terá valor 1
    
    % Calcula a fronteira de pareto do indivíduo.
    c(L.COLPT) = indpar(c,arqsq,L);

    if nh == 1
        % c em referência isolada: deve ser inserido
        if size(arqsq,1) < maxpontos
            % Arquivo tem espaço: insere c.
            arqsq = [arqsq; c];
        else
            % Arquivo está cheio: substitui o pior.
            % Encontra o box do indivíduo
            hy = indbox(c,inf,sqbox,boxsize,L);
            % Indivíduos no mesmo hyperbox (índices)
            hyinds = find(arqsq(:,L.COLHY) == hy);
            sq = size(hyinds,1) + 1;
            c(L.COLHY) = hy;
            c(L.COLSQ) = sq;
            %Encontra o pior
            pr = pior(arqsq,1:size(arqsq,1),@compdist,L);
            cpr = arqsq(pr,:);
            if cpr(L.COLHY) == hy
                cpr(L.COLSQ) = cpr(L.COLSQ) + 1;
            end
            % Compara com o indivíduo
            [m,~] = compdist(1,2,[c;cpr],L);
            if m == 1
                % c é melhor: substitui
                arqsq(pr,:) = c;
            else
                % c é pior: descarta
                descartado = 1;
            end
        end
    else
        % c em nicho ocupado
        % Obtém o outro indivíduo do nicho.
        inref = arqsq(rfinds,:);
        % Compara as distâncias
        if c(L.COLDT) < inref(L.COLDT)
            % c é melhor: substitui
            arqsq(rfinds,:) = c;
        else
            % c é pior: descarta
            descartado = 1;
        end
    end
    
end

if ~descartado
    % c foi inserido na população: refaz os cálculos de dominação
    % e distribuição.
    if ~extrap
        % Observação: este código pode ser otimizado para realizar as
        % atualizações estritamente necessárias, tal como foi feito
        % no caso da população geral.
        no = max(L.COLF) - min(L.COLF) + 1;
        arqsq = pareto(arqsq,L);
        [arqsq,sqbox] = hyperbox(arqsq,inf,boxsize,L);
        [arqsq,sqref] = calcRefPontos(inf,sup,npontos,no,arqsq,L);
        arqsq = selecionaDistribuicao(arqsq,L);
    end
else
    extrap = 0; % Não considera extrapolação para ponto descartado
end
    
end

function m = indice(c,pop,L)
%INDICE Verifica se um indivíduo está em uma população.
%   Verifica se um indivíduo participa de uma população.
%   Considera apenas os campos das variáveis.
%
%   Parâmetros de entrada:
%     - c: indivíduo a considerar;
%     - pop: população;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - m: 0 (não está na população), índice (caso contrário).

m = 0;
for i = 1:size(pop,1)
    if c(L.COLX) == pop(i,L.COLX)
        m = i;
        break;
    end
end

end

function dom = dominacao(c,pop,L)
%DOMINACAO Calcula elementos dominados por e que dominam um indivíduo.
%   Encontra os indivíduos de uma população que dominam um determinado
%   indivíduo e que são dominados por ele.
%
%   Parâmetros de entrada:
%     - c: indivíduo a comparar;
%     - pop: população;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - dom: estrutura de dominação, com os campos:
%         n: número de indivíduos dominados por c;
%         idn: vetor com índices dos indivíduos dominados por c;
%         m: número de indivíduos que dominam c;
%         idm: vetor com índices dos indivíduos que dominam c.

no = max(L.COLF) - min(L.COLF) + 1;
npop = size(pop,1);
elemento = struct('n',0, 'idn',zeros(1,npop), 'm',0, 'idm',zeros(1,npop));
dom = repmat(elemento,1,1);

for i = 1:npop
    % Compara com o i-ésimo indivíduo da população
    difs = c(L.COLF) - pop(i,L.COLF);
    menores = sum(difs < 0);
    iguais = sum(difs == 0);
    maiores = sum(difs > 0);
    
    if maiores == 0 && iguais ~= no
        % c domina i: nenhum maior, eplo menos um menor
        dom.n = dom.n + 1;
        dom.idn(dom.n) = i;
    elseif menores == 0 && iguais ~= no
        % i domina c: nenhum menor, pelo menos um diferente
        dom.m = dom.m + 1;
        dom.idm(dom.m) = i;
    end
end

end

function ext = extrapola(c,inf,sup,L)
%EXTRAPOLA Verifica se um indivíduo extrapola os limites do grid.
%   Verifica se um indivíduo c extrapola os limites do grid de uma
%   determinada população (e, portanto, aumenta sua diversidade).
%
%   Parâmetros de entrada:
%     - c: indivíduo a verificar;
%     - inf: limites inferiores do grid (para todas as dimensões);
%     - sup: limites superiores do grid (para todas as dimensões);
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - ext: 0 (não viola) ou 1 (viola).

  infdif = (c(L.COLF) - inf);
  supdif = (sup - c(L.COLF));
  
  extinf = sum(infdif <= 0);
  extsup = sum(supdif <= 0);
  
  if extinf > 0 || extsup > 0
    ext = 1;
  else
    ext = 0;
  end
  
end

function pt = indpar(c,pop,L)
%INDPAR Calcula a fronteira de Pareto de um indivíduo.
%   Encontra a qual fronteira de Pareto um indivíduo pertence
%   em relação a uma população.
%
%   Parâmetros de entrada:
%     - c: indivíduo a verificar;
%     - pop: população;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - pt: fronteira de Pareto do indivíduo.

dom = dominacao(c,pop,L);

if dom.m == 0
    % Não dominado: primeira fronteira
    pt = 1;
else
    % Dominado: verifica a última fronteira dos dominadores e soma 1.
    pt = 0;
    for i = 1:dom.m
        front = pop(dom.idm(i),L.COLPT);
        if front > pt
            pt = front;
        end
    end
    pt = pt + 1;
end

end

function [rf,dist] = indref(c,refs,L)
%INDREF Calcula o ponto de referência mais próximo de um indivíduo.
%   XXX
%
%   Parâmetros de entrada:
%     - c: indivíduo a verificar;
%     - refs: pontos de referência da população
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - rf: ponto de referência do indivíduo
  
npontos = length(refs);

prox = 1;
distprox = norm(c(L.COLF) - refs(1).ponto);
    
for j = 2:npontos
    dist = norm(c(L.COLF) - refs(j).ponto);
    if dist < distprox
        prox = j;
        distprox = dist;
    end
end
    
rf = prox;
dist = distprox;

end

function hb = indbox(c,inf,boxes,boxsize,L)
%INDBOX Calcula o hyperbox de um indivíduo.
%   Encontra o hyperbox ocupado por um indivíduo em uma determinada
%   população.
%
%   Parâmetros de entrada:
%     - c: indivíduo a verificar;
%     - inf: limites inferiores do grid (para todas as dimensões);
%     - sup: limites superiores do grid (para todas as dimensões);
%     - boxes: boxes da população
%     - resolucao: tamanho do grid (o mesmo para todas as dimensões);
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - hb: hyperbox do indivíduo
  
  hbox = floor((c(L.COLF) - inf) ./ boxsize);
  
  ind = 0;
  for k = 1:length(boxes)
      if boxes(k).hbox == hbox
          ind = k;
          break
      end
  end
  
  hb = ind;
  
end

function [pop,ref] = insereref(pop,ind,ref,L)
%INSEREBOX Insere um indivíduo em um ponto de referência.
%   Insere um indivíduo em um ponto de referência e atualiza os dados
%   do ponto e dos indivíduos ligados a ele.
%
%   Parâmetros de entrada:
%     - pop: população;
%     - ind: índice do indivíduo na população;
%     - ref: pontos de referência da população
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - pop: XXX
%     - ref: pontos de referência atualizados.

c = pop(ind,:);

[iref,dist] = indref(c,ref,L);

ref(iref).n = ref(iref).n + 1;
ref(iref).ids(ref(iref).n) = ind;

for i = 1:ref(iref).n
    pop(ref(iref).ids(i), L.COLNH) = ref(iref).n; 
end

pop(ind,L.COLDT) = dist;

end

function [pop,ref] = removeref(pop,ind,ref,L)
%REMOVEBOX Remove um indivíduo de um ponto de referência.
%   Remove um indivíduo de um ponto de referência e atualiza os dados
%   do ponto.
%
%   Parâmetros de entrada:
%     - pop: população;
%     - ind: índice do indivíduo na população;
%     - ref: pontos de referência da população;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - pop: XXX
%     - ref: pontos de referência atualizados.

iref = (pop(ind,L.COLRF));
for i = 1:ref(iref).n
    if ref(iref).ids(i) == ind
        break;
    end
end

for j = (i+1):ref(iref).n
    ref(iref).ids(j-1) = ref(iref).ids(j);
end

ref(iref).n = ref(iref).n - 1;

for i = 1:ref(iref).n
    pop(ref(iref).ids(i),L.COLNH) = ref(iref).n;
end

end

function [pop,box] = inserebox(pop,ind,inf,box,boxsize,L)
%INSEREBOX Insere um indivíduo em um box.
%   Insere um indivíduo em um box e atualiza os dados do box.
%
%   Parâmetros de entrada:
%     - pop: população;
%     - ind: índice do indivíduo na população;
%     - inf: XXX
%     - box: boxes da população
%     - boxsize: XXX
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - pop: XXX
%     - box: boxes atualizados.

c = pop(ind,:);

ibox = indbox(c,inf,box,boxsize,L);

box(ibox).n = box(ibox).n + 1;
box(ibox).ids(box(ibox).n) = ind;

for i = 1:box(ibox).n
    pop(box(ibox).ids(i), L.COLSQ) = box(ibox).n; 
end

end

function [pop,box] = removebox(pop,ind,box,L)
%REMOVEBOX Remove um indivíduo de uma box.
%   Remove um indivíduo de um box e atualiza os dados do box.
%
%   Parâmetros de entrada:
%     - pop: população;
%     - ind: índice do indivíduo na população;
%     - box: boxes da população
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - pop: XXX
%     - box: boxes atualizados.

indbox = (pop(ind,L.COLHY));
for i = 1:box(indbox).n
    if box(indbox).ids(i) == ind
        break;
    end
end

for j = (i+1):box(indbox).n
    box(indbox).ids(j-1) = box(indbox).ids(j);
end

box(indbox).n = box(indbox).n - 1;

for i = 1:box(indbox).n
    pop(box(indbox).ids(i),L.COLSQ) = box(indbox).n;
end

end

function res = calcResolucao(pop,precisao,L)
%CALCRESOLUCAO Calcula a resolucao do grid de hyperboxes.
%   Dada uma precisao, encontra a resolução do grid de
%   hyperboxes para uma determinada população.
%
%   Parâmetros de entrada:
%     - pop: população;
%     - precisao: precisão desejada;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - res: vetor com a resolução para cada dimensão.

% Encontra os valores máximos e mínimos das funções objetivo.
zmax = max(pop(:,L.COLF),[],1);
zmin = min(pop(:,L.COLF),[],1);

res = floor((zmax - zmin) / precisao) + 1;

end

function ind = pior(pop,inds,comp,L)
%PIOR Encontra o pior indivíduo de uma (sub)população.
%   Encontra o pior indivíduo dentre um conjunto de indivíduos
%   (cujos índices são fornecidos) de uma população.
%   Utiliza uma função fornecida para fazer as comparações.
%
%   Parâmetros de entrada:
%     - pop: população;
%     - inds: índices dos indivíduos a considerar;
%     - comp: função de comparação.
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - ind: índice do pior indivíduo.

n = size(inds,2);

ind = inds(1);
for k = 2:n
    [~,ind] = comp(ind,inds(k),pop,L);
end

end

function pop = remove(pop,inds,n,L)
%REMOVE remove indivíduos de uma população.
%   Remove de uma população os indivíduos cujos índices são fornecidos.
%   Observação: pressupõe que os índices estão em ordem crescente.
%
%   Parâmetros de entrada:
%     - pop: população;
%     - inds: índices dos indivíduos a remover;
%     - n: quantidade de elementos a remover;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - pop: população atualizada.

% Aloca espaço para a nova população
npop = size(pop,1);
rpop = zeros(npop-n,L.NC);

% Efetua as remoções
inicio = 1;
c = 1;
for i = 1:n
    % Copia elementos anteriores ao índice
    for k = inicio:inds(i)-1
        rpop(c,:) = pop(k,:);
        c = c+1;
    end
    
    % Novo início é após o elemento removido
    inicio = inds(i) + 1;
end

% Copia os elementos restantes
for k = inicio : npop
    rpop(c,:) = pop(k,:);
    c = c + 1;
end

pop = rpop;

end

function f = cruzamento(ind1,ind2,pc,L)
%CRUZAMENTO Realiza o cruzamento entre dois indivíduos.
%   Realiza o cruzamento entre dois indivíduos e retorna o filho
%   resultante.
%
%   Parâmetros de entrada:
%     - ind1: primeiro indivíduo;
%     - ind2: segundo indivíduo;
%     - pc: probabilidade de cruzamento;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - f: filho gerado pelo cruzamento.

% Determina o melhor e o pior pai.
pais = [ind1; ind2];
[m,p] = compara(1,2,pais,L);
mp = pais(m,:);
pp = pais(p,:);

% Verifica se o cruzamento deve ser realizado;
% caso contrário, apenas copia o melhor pai.
if rand > pc
    f = mp;
else
    f1 = zeros(1,L.NC);
    %f2 = zeros(1,L.NC);
    
    % Determina os coeficientes
    nvar = max(L.COLX) - min(L.COLX) + 1;
    kcross = randi(nvar);
    apol = 0.5 * rand + 0.5;
    %a = 1.2 * rand - 0.1;
    dir = randi(2) - 1;
    
    % Realiza o cruzamento
    if dir == 0
        f1(1:kcross) = apol * mp(1:kcross) + (1-apol) * pp(1:kcross);
        f1((kcross+1):nvar) = mp((kcross+1):nvar);
        %f2(1:kcross) = (1-a) * mp(1:kcross) + a * pp(1:kcross);
        %f2((kcross+1):nvar) = pp((kcross+1):nvar);
    else
        f1(1:(kcross-1)) = mp(1:(kcross-1));
        f1(kcross:nvar) = apol*mp(kcross:nvar) + (1-apol)*pp(kcross:nvar);
        %f2(1:(kcross-1)) = pp(1:(kcross-1));
        %f2(kcross:nvar) = (1-a) * mp(kcross:nvar) + a * pp(kcross:nvar);
    end
    
    f = f1;
end

end

function gamma = perturbacao(c,pm,xmin,xmax,L)
%PERTURBACAO Calcula o vetor de perturbações.
%   O vetor de perturbações corresponde aos valores a serem
%   adicionados às variáveis do indivíduo para a implementação
%   da mutação.
%
%   Parâmetros de entrada:
%     - c: indivíduo a perturbar;
%     - pm: probabilidade de mutação;
%     - xmin: valor mínimo de uma variável;
%     - xmax: valor máximo de uma variável;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - gamma: vetor de perturbações.

% Inicializa o vetor de perturbações.
nvar = max(L.COLX) - min(L.COLX) + 1;
gamma = zeros(1,nvar);

if rand < pm
    kmut = randi(nvar);
    range = xmax - xmin;
    dir = randi(2) - 1;
    if dir == 0
        faixa = 1:kmut;
    else
        faixa = kmut:nvar;
    end
    
    for i = faixa
        beta = (2 * rand - 1);
        gamma(i) = 0.05 * beta * range;
    end
end

end

function [pc, pm] = calcProbs(pop, pc, pm, L)
%CALCPROBS Calcula as probabilidades de cruzamento e mutação.
%   Calcula as probabilidades de cruzamento e mutação da presente geração.
%   Essas probabilidades são influenciadas pelo mdg.
%
%   Parâmetros de entrada:
%     - pop: array com os indivíduos da população;
%     - pc: probabilidade de cruzamento atual;
%     - pm: probabilidade de mutação atual;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - pc: nova probabilidade de cruzamento;
%     - pm: nova probabilidade de mutação.


vinf = 0.2;
vmax = 0.7;
kc = 1.2;
km = 1.2;

mdg = calcMdg(pop,L);

if mdg < vinf
  pc = pc / kc;
  pm = pm * km;
elseif mdg > vmax
  pc = pc * kc;
  pm = pm / km;
end

end

function [mdg] = calcMdg(pop,L)
%CALCMDG Calcula a medida da diversidade genética.
%   A medida da diversidade genética da população influencia as
%   probabilidades de cruzamento e mutação. O cálculo da medida
%   baseia-se no valor agregado das funções objetivo.
%
%   Parâmetros de entrada:
%     - pop: array com os indivíduos da população;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%       - mdg: medida de diversidade genética da população.

npop = size(pop,1);
fit = pop(:,L.COLAG);
fmed = sum(fit) / npop;
fmax = max(fit);

mdg = (fmax - fmed) / fmax;

end

function pop = calcFos(pop,objFunc,L)
%CALCFOS Calcula as funções objetivo de uma população.
%   Utiliza a função fornecida como parâmetro para calcular
%   os valores para cada indivíduo.
%
%   Parâmetros de entrada:
%     - pop: array com os indivíduos da população;
%     - objfunc: função objetivo a utilizar;
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - pop: população atualizada.

npop = size(pop,1);
%nvar = max(L.COLX) - min(L.COLX) + 1;
no = max(L.COLF) - min(L.COLF) + 1;

for i = 1:npop
    pop(i,L.COLF) = objFunc(pop(i,L.COLX),no);
end

end

%---------------------

function f = dtlz1(x,m)

% x = indivíduo
% m = número de funções objetivo

x = x(:);

n = length(x);
k = n - m + 1;

s = 0;
for i = m:n
    s = s + (x(i)-0.5)^2 - cos(20*pi*(x(i)-0.5));
end
g = 100*(k+s);

f = zeros(1,m);

f(1) = 0.5 * prod(x(1:m-1)) * (1+g);
for i = 2:m-1
    f(i) = 0.5 * prod(x(1:m-i)) * (1-x(m-i+1)) * (1+g);
end
f(m) = 0.5 * (1-x(1)) * (1+g);

end

%---------------------


function f = dtlz2(x,m)

% x = indivíduo
% m = número de funções objetivo

x = x(:);

n = length(x);
k = n - m + 1;

%display(['n = ', num2str(n)]);
%display(['m = ', num2str(m)]);
%display(['k = ', num2str(k)]);

s = 0;
for i = m:n
    s = s + (x(i)-0.5)^2;
end
g = s;

cosx = cos(x*pi/2);
sinx = sin(x*pi/2);

f = zeros(1,m);

f(1) =  (1+g) * prod(cosx(1:m-1));
for i = 2:m-1
    f(i) = (1+g) * prod(cosx(1:m-i)) * sinx(m-i+1);
end
f(m) = (1+g) * sinx(1);

%F = f(:);

end
