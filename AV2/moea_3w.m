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
npop = 10; % XXX Cálculo?

% Estabelecimento da polulação inicial.
pop = popinit(npop,xmin,xmax,L);

% Avaliação dos elementos da população inicial.
%fos = @dtlz1; % XXX Parâmetro?
fos = @dtlz2;
pop = fos(pop,nvar,no,L);
ncal = ncal - npop;

pop = agregacao(pop,L); % cálculo do valor agregado
pop = pareto(pop,L);    % cálculo da fronteira de pareto
resolucao = 5; % XXX Parâmetro?
[pop, pinf, psup] = hyperbox(pop,resolucao,L); % cálculo do hyperbox

% Cópia dos elementos não dominados da população para o arquivo arqnd
arqnd = naodominado(pop,L);
ndinf = pinf;
ndsup = psup;

% Cópia dos elementos não agregados da população para o arquivo arqsq
arqsq = naoaglomerado(pop,L);
sqinf = pinf;
sqsup = psup;

display(pop);

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
        c = fos(c,nvar,no,L);
        c = agregacao(c,L);
        
        % Teste
        display(c);
        [pop,pinf,psup] = atualizapop(c,pop,pinf,psup,resolucao,L);
        display(pop);
    end
    
    % Foram avaliadas as funções objetivo de três indivíduos
    ncal = ncal - 3;
end

%for i = 1:npop
%    display(decodifica(pop(i,L.COLHY),no,resolucao));
%end

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

function [pop, inf, sup] = hyperbox(pop,resolucao,L)
%HYPERBOX Calcula o hyperbox de cada indivíduo.
%   Para cada indivíduo da população, calcula o hyperbox ao qual
%   pertence e o fator de aglomeração. Retorna também os limites
%   inferior e superior do grid para cada dimensão.
%
%   Parâmetros de entrada:
%     - pop: array contendo a população inicial;
%     - resolucao: tamanho do grid (o mesmo para todas as dimensões).
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - pop: população inicial hyperbox e aglomeração atualizados;
%     - inf: limites inferiores do grid para cada dimensão;
%     - sup: limites superiores do grid para cada dimensão.

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
%   pareto > fator de agregação > valor agregado.
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
%   fator de agregação > pareto > valor agregado.
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

if ind1(L.COLSQ) < ind2(L.COLSQ)
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

function [pop,pinf,psup] = atualizapop(c,pop,pinf,psup,resolucao,L)
%ATUALIZAPOP Atualiza a população com um indivíduo.
%   Tenta inserir um indivíduo em uma população. Se o indivíduo
%   não for descartado, refaz os cálculos das fronteiras de Pareto
%   e do hyperbox.
%
%   Parâmetros de entrada:
%     - c: indivíduo a inserir;
%     - pop: população;
%     - inf: limites inferiores do grid para cada dimensão;
%     - sup: limites superiores do grid para cada dimensão;
%     - resolucao: tamanho do grid.
%     - L: layout de um indivíduo.
%
%   Parâmetros de saída:
%     - pop: população atualizada.
%     - inf: limites inferiores atualizados;
%     - sup: limites superiores atualizados.

% Calcula os indivíduos de pop dominados e que dominam c.
dom = dominacao(c,pop,L);

descartado = 0;

if dom.n > 0
    % c domina um grupo de indivíduos:
    % substitui o pior deles por c.
    inds = dom.idn(1:dom.n);
    pr = pior(pop,inds,@compara,L);
    pop(pr,:) = c;
elseif dom.m > 0
    % c é dominado e descartado.
    descartado = 1;
else
    % c não domina nem é dominado:
    % substitui o pior indivíduo da população.
    inds = 1:size(pop,1);
    pr = pior(pop,inds,@compara,L);
    pop(pr,:) = c;
end

if ~descartado
    % c foi inserido na população: refaz os cálculos de dominação
    % e distribuição.
    pop = pareto(pop,L);
    [pop, pinf, psup] = hyperbox(pop,resolucao,L); % cálculo do hyperbox
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
    %display(['menor = ', num2str(menores), ' maior = ', num2str(maiores), ' igual = ', num2str(iguais)]);
    
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
%       - pop: array com os indivíduos da população;
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
