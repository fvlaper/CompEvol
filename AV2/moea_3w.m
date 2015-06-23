function [ ps ] = moea_3w( ncal, npop, nvar, no, precisao, objfunc )
%MOEA_3W Multi Objective Evolutionary Algorithm - 3 Way
%   Esta fun��o implementa um algoritmo para otimiza��o de fun��es
%   Multi-objetivo. XXX Explicar o algoritmo.
%
%   Par�metros de entrada:
%     - ncal: n�mero m�ximo de c�lculos das fun��es objetivo;
%     - npop: tamanho da popula��o;
%     - nvar: n�mero de vari�veis;
%     - no: n�mero de fun��es objetivo;
%     - precisao: XXX
%     - objfunc: handle da fun��o objetivo
%
%   Par�metros de sa�da:
%     - ps: conjunto de Pareto.

% Layout de array de popula��o: array contendo uma linha por
% indiv�duo. Para cada linha:
%   - Vari�veis de decis�o: nvar colunas.
%   - Fun��es objetivo: no colunas.
%   - Valor agregado: uma coluna.
%   - Fronteira de Pareto: uma coluna;
%   - Hyperbox: uma coluna;
%   - Fator de aglomera��o (squeeze factor): uma coluna;
%   - Ponto de refer�ncia: uma coluna;
%   - Fator de nicho: uma coluna;
%   - Dist�ncia ao ponto de refer�ncia: uma coluna.
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

% Dados iniciais da popula��o: faixa das vari�veis e 
% quantidade de indiv�duos.
xmin = 0; xmax = 1;

% Estabelecimento da popula��o inicial.
pop = popinit(npop,xmin,xmax,L);

% Avalia��o dos elementos da popula��o inicial.
pop = calcFos(pop,objfunc,L);
ncal = ncal - npop;

pop = agregacao(pop,L); % c�lculo do valor agregado
pop = pareto(pop,L);    % c�lculo da fronteira de pareto
[inf,sup,boxsize] = calParEspaco(pop,precisao,L); % XXX
[pop, pbox] = hyperbox(pop,inf,boxsize,L); % c�lculo do hyperbox
[pop, pref] = calcRefPontos(inf,sup,npop,no,pop,L); % c�lculo dos pontos de refer�ncia

% C�pia dos elementos n�o dominados da popula��o para o arquivo arqnd
arqnd = naodominado(pop,L);
[arqnd,ndbox] = hyperbox(arqnd,inf,boxsize,L);
[arqnd,ndref] = calcRefPontos(inf,sup,npop,no,arqnd,L); % c�lculo dos pontos de refer�ncia

% C�pia dos elementos n�o agregados da popula��o para o arquivo arqsq
arqsq = selecionaDistribuicao(pop,L);
[arqsq,sqbox] = hyperbox(arqsq,inf,boxsize,L);
[arqsq,sqref] = calcRefPontos(inf,sup,npop,no,arqsq,L); % c�lculo dos pontos de refer�ncia

% Probabilidades iniciais de cruzamento e muta��o
pc = 0.6;
pm = 0.05;

% La�o principal
while ncal >= 3
    % Seleciona os pais para os cruzamentos (um de cada popula�ao/arquivo)
    p1 = seleciona(pop,@compara,L);
    p2 = seleciona(arqnd,@compara,L);
    p3 = seleciona(arqsq,@compdist,L);
    
    % Calcula as probabilidades de cruzamento e muta��o
    [pc,pm] = calcProbs(pop,pc,pm,L);

    % Realiza os cruzamentos para gerar tr�s filhos.
    c1 = cruzamento(p1,p2,pc,L);
    c2 = cruzamento(p1,p3,pc,L);
    c3 = cruzamento(p2,p3,pc,L);
    filhos = [c1;c2;c3];
    
    % Tratamentos dos filhos
    for i = 1:size(filhos,1)
        
        c = filhos(i,:);
        
        % Realiza muta��es
        gamma = perturbacao(c,pm,xmin,xmax,L);
        c(L.COLX) = c(L.COLX) + gamma;

        % Garante que as condi��es de contorno sejam respeitadas
        c(:,L.COLX) = min(xmax, max(xmin,c(:,L.COLX)));
        
        % Calcula as fun��es objetivo e o valor agregado
        c = calcFos(c,objfunc,L);
        c = agregacao(c,L);
        
        % Tenta inserir c na popula��o e nos arquivos
        [pop,pbox,pref,pextrapola] = atualizapop(c,pop,inf,sup,pbox,pref,boxsize,npop,L);

        [arqnd,ndbox,ndref,ndextrapola] = atualizadom(c,arqnd,inf,sup,ndbox,ndref,...
                                          npop,boxsize,npop,L);
        [arqsq,sqbox,sqref,sqextrap] = atualizadis(c,arqsq,inf,sup,sqbox,sqref,...
                                          npop,boxsize,npop,L);
        
        % Algum novo ponto extrapolou o espa�o; recalcula
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
  
    % Foram avaliadas as fun��es objetivo de tr�s indiv�duos
    ncal = ncal - 3;
    display(ncal);
end

% Retorno do resultado (elementos na primeira fronteira de Pareto)
ps = arqnd;
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

function pop = agregacao(pop,L)
%AGREGACAO Calcula o valor agregado das fun��es objetivo.
%   Calcula o valor agregado das fun��es objetivo utilizando
%   o m�todo da soma. Considera o mesmo pelo para todas as fun��es.
%
%   Par�metros de entrada:
%     - pop: array contendo a popula��o inicial;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - pop: popula��o inicial com valor agregado atualizado.

no = max(L.COLF) - min(L.COLF) + 1; % n�mero de fun��es objetivo
pop(:,L.COLAG) = sum(pop(:,L.COLF),2) ./ no;

end

function pop = pareto(pop,L)
%PARETO Calcula a fronteira de Pareto.
%   Para cada indiv�duo da popula��o, calcula a que fronteira de
%   Pareto ele pertence.
%   Inspirado em XXX.
%
%   Par�metros de entrada:
%     - pop: array contendo a popula��o inicial;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - pop: popula��o inicial com fronteira de Pareto atualizada.

npop = size(pop,1); % n�mero de indiv�duos
no = max(L.COLF) - min(L.COLF) + 1; % n�mero de fun��es objetivo
f = 1;  % f-�sima fronteira de Pareto

% Fronteira: vetor com uma estrutura por fronteira de Pareto:
%   n: n�mero de indiv�duos na fronteira;
%   ids: �ndices dos indiv�duos na fronteira.
indfront = struct('n',0, 'ids',zeros(1,npop));
front = repmat(indfront,1,npop);

% Domina��o: vetor com uma estrutura por indiv�duo i:
%   n: n�mero de indiv�duos que dominam i;
%   nd: n�mero de indiv�duos dominados por i;
%   d: �ndices dos indiv�duos dominados por i.
elemento = struct('n',0, 'nd',0, 'd',zeros(1,npop));
dominacao = repmat(elemento,1,npop);

% Compara cada indiv�duo com todos os demais para determinar
% a domina��o.
for i = 1:npop
    for j = 1:npop
        if i ~= j
            % Compara indiv�duos i e j
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
    
    % Trata os indiv�duos na primeira fronteira.
    if dominacao(i).n == 0
        front(f).n = front(f).n + 1;
        front(f).ids(front(f).n) = i;
        pop(i,L.COLPT) = 1;
    end
end

% Tratamento das demais fronteiras
while front(f).n ~= 0
    % Percorre os elementos da f-�sima fronteira, removendo-os
    % e fazendo os c�lculos para a (f+1)-�sima fronteira.
    for i = 1:front(f).n
        individuo = front(f).ids(i);
        % Trata os indiv�duos dominados pelo indiv�duo atual.
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
%CALCPARESPACO Calcula par�metros do espa�o de pontos.
%   XXX
%
%   Par�metros de entrada:
%     - pop: array contendo a popula��o;
%     - precisao: XXX
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - inf: limites inferiores do grid para cada dimens�o;
%     - sup: limites superiores do grid para cada dimens�o;
%     - boxsize: XXX

% Encontra os valores m�ximos e m�nimos das fun��es objetivo.
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
%HYPERBOX Calcula o hyperbox de cada indiv�duo.
%   Para cada indiv�duo da popula��o, calcula o hyperbox ao qual
%   pertence e o fator de aglomera��o.
%
%   Par�metros de entrada:
%     - pop: array contendo a popula��o;
%     - inf: XXX;
%     - boxsize: XXX;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - pop: popula��o atualizados;
%     - pbox: vetor de estruturas com as informa��es dos boxes:
%             hbox: coordenadas do box;
%             n: n�mero de elementos do box;
%             ids: �ndices dos indiv�duos do box.

% Calcula o hyperbox (para cada dimens�o) e o fator de aglomeracao
npop = size(pop,1); % n�mero de indiv�duos
no = max(L.COLF) - min(L.COLF) + 1; % n�mero de fun��es objetivo
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

% Registra o valor de aglomera��o para cada indiv�duo.
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
%   Par�metros de entrada:
%     - inf: XXX
%     - sup: XXX
%     - npontos: XXX
%     - no: XXX
%     - pop: XXX
%     - L: XXX
%
%   Par�metros de sa�da:
%     - pop: XXX
%     - refpontos: XXX

% Observa��o: o c�lculo do n�mero de subdivis�es para cada
% dimens�o abaixo pode gerar um n�mero muito grande de pontos
% de refer�ncia para problemas com muitos objetivos. Deve ser
% estudada uma forma de estabelecer um limite superior para esse
% n�mero para n�o comprometer a efici�ncia do algoritmo.
n = max(2,round(nthroot(npontos,no))); % nro de subdivisoes por dimens�o

npop = size(pop,1);

delta = (sup - inf) / n; % espa�amento dos pontos (por dimens�o)
npontos = (n+1) ^ no; % nro total de pontos de refer�ncia

strucponto = struct ('ponto',zeros(1,no), 'n', 0, 'ids',zeros(1,npop));
refpontos = repmat(strucponto,1,npontos);

% Preenche os pontos de refer�ncia
ponto = inf; % Come�a pelo "canto inferior esquerdo"
for k = 1:npontos
    refpontos(k).ponto = ponto; % grava o ponto
    
    % Faz incremento em uma das dimens�es para o pr�ximo ponto
    tol = 1e-6;
    d = 1; % dimens�o a incrementar
    incrementado = 0; % incremento bem sucedido?
    while(~incrementado)
        ponto(d) = ponto(d) + delta(d); % incrementa
        if ponto(d) <= sup(d) + tol % se superou m�ximo, muda dimens�o
            incrementado = 1;
        else
            ponto(d) = inf(d);
            d = d+1;
        end
        
        if d > no  % seguran�a
            break;
        end
    end
end

% Encontra o ponto de refer�ncia mais pr�ximo de cada indiv�duo
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
    
    % XXX Gravar aqui informa��es no indiv�duo
    pop(i,L.COLDT) = distprox;
    refpontos(prox).n = refpontos(prox).n + 1;
    refpontos(prox).ids(refpontos(prox).n) = i;
end

% Grava as informa��es nos indiv�duos
for k = 1:npontos
    for i = 1:refpontos(k).n
        pop(refpontos(k).ids(i), L.COLRF) = k;
        pop(refpontos(k).ids(i), L.COLNH) = refpontos(k).n;        
    end
end

end

function p = seleciona(pop,comp,L)
%SELECIONA Seleciona um indiv�duo da popula��o.
%   Seleciona um indiv�duo para participar dos cruzamentos.
%   Utiliza torneio bin�rio.
%
%   Par�metros de entrada:
%     - pop: popula��o;
%     - comp: handle da fun��o de compara��o;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - p: indiv�duo selecionado.

npop = size(pop,1); % n�mero de indiv�duos

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
%COMPARA Compara dois indiv�duos e retorna-os ordenados.
%   Utiliza crit�rios normais de compara��o:
%   pareto > fator de nicho> fator de agrega��o > valor agregado.
%
%   Par�metros de entrada:
%     - i1: �ndice do indiv�duo 1;
%     - i2: �ndice do indiv�duo 2;
%     - pop: popula��o;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - melhor: �ndice do melhor indiv�duo;
%     - pior: �ndice do pior indiv�duo.

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
%COMPDIST Compara dois indiv�duos e retorna-os ordenados.
%   Utiliza crit�rio alternativo (prioriza distribui��o):
%   fator de nicho > fator de agrega��o > pareto > valor agregado.
%
%   Par�metros de entrada:
%     - i1: �ndice do indiv�duo 1;
%     - i2: �ndice do indiv�duo 2;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - melhor: melhor indiv�duo;
%     - pior: pior indiv�duo.

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
%COMPDIST Compara dois indiv�duos e retorna-os ordenados.
%   Utiliza dist�ncia a ponto de refer�ncia para compara��o.
%
%   Par�metros de entrada:
%     - i1: �ndice do indiv�duo 1;
%     - i2: �ndice do indiv�duo 2;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - melhor: melhor indiv�duo;
%     - pior: pior indiv�duo.

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
%NAODOMINADO Encontra os elementos n�o dominados de uma populacao.
%   Encontra os elementos n�o dominados de uma popula��o e retorna-os
%   em uma subpopula��o separada.
%
%   Par�metros de entrada:
%     - pop: array contendo a popula��o inicial;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - arqnd: subpopula��o de elementos n�o dominados.

idx = (pop(:,L.COLPT) == 1);
arqnd = pop(idx,:);

end

function arqsq = selecionaDistribuicao(pop,L)
%SelecionaDistribuicao Seleciona os elementos do arquivo de distribui��o.
%   Examina os elementos de uma popula��o e seleciona os melhores de
%   cada nicho para compor o arquivo de distribui��o.
%
%   Par�metros de entrada:
%     - pop: array contendo a popula��o inicial;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - arqsq: subpopula��o de elementos n�o aglomerados.

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
%ATUALIZAPOP Atualiza a popula��o com um indiv�duo.
%   Tenta inserir um indiv�duo em uma popula��o. XXX
%
%   Par�metros de entrada:
%     - c: indiv�duo a inserir;
%     - pop: popula��o;
%     - inf: limites inferiores do grid para cada dimens�o;
%     - sup: limites superiores do grid para cada dimens�o;
%     - pbox: box da popula��o;
%     - pref: XXX
%     - boxsize: XXX
%     - npontos: XXX
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - pop: popula��o atualizada.
%     - pbox: XXX
%     - pbof: XXX
%     - extrapola: XXX

% Verifica se indiv�duo n�o existe na popula��o
%idx = indice(c,pop,L);
%if idx > 0
%    display('c j� existe: descartado');
%    return;
%end

% Calcula os indiv�duos de pop dominados e que dominam c.
dom = dominacao(c,pop,L);

extrap = extrapola(c,inf,sup,L);
descartado = 0;

if dom.n > 0
    % c domina um grupo de indiv�duos:
    % substitui o pior deles por c.
    inds = dom.idn(1:dom.n);
    pr = pior(pop,inds,@compara,L);
    % XXX atualizar pbox e pref
    [pop,pbox] = removebox(pop,pr,pbox,L);
    [pop,pref] = removeref(pop,pr,pref,L);
    pop(pr,:) = c;
elseif dom.m > 0
    % c � dominado e descartado.
    descartado = 1;
else
    % c n�o domina nem � dominado:
    % substitui o pior indiv�duo da popula��o.
    inds = 1:size(pop,1);
    pr = pior(pop,inds,@compara,L);
    % XXX atualizar pbox e pref
    [pop,pbox] = removebox(pop,pr,pbox,L);
    [pop,pref] = removeref(pop,pr,pref,L);
    pop(pr,:) = c;
    % XXX atualiza pbox e pref
end

if ~descartado
    % c foi inserido na popula��o: refaz os c�lculos de domina��o
    % e distribui��o.
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
            % Caiu em novo ponto de refer�ncia; recalcula
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
    extrap = 0; % N�o considera extrapola��o para ponto descartado
end

end

function [arqnd,ndbox,ndref,extrap] = atualizadom(c,arqnd,inf,sup,ndbox,ndref, ...
                                           maxpontos,boxsize,npontos,L)
%ATUALIZAPOP Atualiza o arquivo de domina��o com um indiv�duo.
%   Tenta inserir um indiv�duo no arquivo de domina��o popula��o.
%   Se o indiv�duo n�o for descartado, refaz os c�lculos 
%   das fronteiras de Pareto e do hyperbox.
%
%   Par�metros de entrada:
%     - c: indiv�duo a inserir;
%     - arqnd: popula��o (arquivo de domina��o);
%     - inf: limites inferiores do grid para cada dimens�o;
%     - sup: limites superiores do grid para cada dimens�o;
%     - ndbox: boxes do arquivo;
%     - ndref: XXX
%     - max: n�mero m�ximo de elementos em no arquivo;
%     - boxsize: XXX
%     - npontos: XXX
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - arqnd: popula��o atualizada.
%     - ndbox: XXX
%     - ndref: XXX
%     - extrap: XXX

extrap = 0;

% Verifica se indiv�duo n�o existe na popula��o
idx = indice(c,arqnd,L);
if idx > 0
    %display('c j� existe: descartado');
    return;
end

% Calcula os indiv�duos de pop dominados e que dominam c.
dom = dominacao(c,arqnd,L);

extrap = extrapola(c,inf,sup,L);
descartado = 0;

if dom.m > 0
    % c � dominado por algum indiv�duo do arquivo e � descartado
    descartado = 1;
elseif dom.n > 0
    % c domina indiv�duos do arquivo: remove os dominados e insere c
    arqnd = remove(arqnd,dom.idn,dom.n,L);
    arqnd = [arqnd; c];
elseif size(arqnd,1) < maxpontos
    % c n�o domina nem � dominado e h� espa�o no arquivo: insere c
    arqnd = [arqnd; c];
else
    % c n�o domina nem � dominado; n�o h� espa�o no arquivo
    % Verifica se c estende o grid do arquivo
    if extrap
        % Substitui o pior elemento do grid por c
        pr = pior(arqnd,1:size(arqnd,1),@compara,L);
        arqnd(pr,:) = c;
    else
        % Encontra o nicho e a dist�ncia do elemento
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
                % c � melhor: substitui
                idx = indice(cpr,arqnd,L);
                arqnd(idx,:) = c;
                inserido = 1;
            end
        end
        
        % Se n�o foi substitu�do, tenta contra o pior do arquivo
        if ~inserido
            % Encontra o box do indiv�duo
            hy = indbox(c,inf,ndbox,boxsize,L);
            % Indiv�duos no mesmo hyperbox (�ndices)
            hyinds = find(arqnd(:,L.COLHY) == hy);
            sq = size(hyinds,1) + 1;
            c(L.COLHY) = hy;
            c(L.COLSQ) = sq;
            % Encontra o pior da popula��o
            pr = pior(arqnd,1:size(arqnd,1),@compara,L);
            cpr = arqnd(pr,:);
            if cpr(L.COLHY) == hy
                cpr(L.COLSQ) = cpr(L.COLSQ) + 1;
            end
            % Compara com o indiv�duo
            [m,~] = compara(1,2,[c;cpr],L);
            if m == 1
                % c � melhor: substitui
                arqnd(pr,:) = c;
            else
                % c � pior: descarta
                descartado = 1;
            end
        end
        
    end
end

if ~descartado
    % c foi inserido na popula��o: refaz os c�lculos de domina��o
    % e distribui��o.
    if ~extrap
        % Observa��o: este c�digo pode ser otimizado para realizar as
        % atualiza��es estritamente necess�rias, tal como foi feito
        % no caso da popula��o geral.
        no = max(L.COLF) - min(L.COLF) + 1;
        arqnd = pareto(arqnd,L);
        [arqnd,ndbox] = hyperbox(arqnd,inf,boxsize,L);
        [arqnd,ndref] = calcRefPontos(inf,sup,npontos,no,arqnd,L);
    end
else
    extrap = 0; % N�o considera extrapola��o para ponto descartado
end
    
end

function [arqsq,sqbox,sqref,extrap] = atualizadis(c,arqsq,inf,sup,sqbox,sqref, ...
                                           maxpontos,boxsize,npontos,L)
%ATUALIZADIS Atualiza o arquivo de domina��o com um indiv�duo.
%   Tenta inserir um indiv�duo no arquivo de distribui��o.
%   Se o indiv�duo n�o for descartado, refaz os c�lculos 
%   das fronteiras de Pareto e do hyperbox.
%
%   Par�metros de entrada:
%     - c: indiv�duo a inserir;
%     - arqsq: popula��o (arquivo de distribui��o);
%     - inf: limites inferiores do grid para cada dimens�o;
%     - sup: limites superiores do grid para cada dimens�o;
%     - sqbox: boxes do arquivo;
%     - sqref: XXX
%     - maxpontos: n�mero m�ximo de elementos no arquivo;
%     - boxsize: XXX
%     - npontos: n�mero de pontos de refer�ncia;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - arqsq: popula��o atualizada.
%     - sqbox: boxes atualizados;
%     - sqref: XXX
%     - extrap: XXX

extrap = 0;

% Verifica se indiv�duo n�o existe na popula��o
idx = indice(c,arqsq,L);
if idx > 0
    %display('c j� existe: descartado');
    return;
end

extrap = extrapola(c,inf,sup,L);
descartado = 0;

if extrap
    % c estende grid do arquivo: deve ser inserido.
    if size(arqsq,1) < maxpontos
        % Arquivo tem espa�o: insere c.
        arqsq = [arqsq; c];
    else
        % Arquivo est� cheio: substitui o pior.
        pr = pior(arqsq,1:size(arqsq,1),@compdist,L);
        arqsq(pr,:) = c;
    end
else
    % C n�o estende grid: verifica de deve ser inserido
    % Calcula o nicho e a dist�ncia de c.
    [iref,dist] = indref(c,sqref,L);
    c(L.COLRF) = iref;
    c(L.COLDT) = dist;
    
    % Procura por elementos do mesmo nicho
    rfinds = find(arqsq(:,L.COLRF) == iref);
    nh = length(rfinds) + 1;
    c(L.COLNH) = 1; % Se inserido, ter� valor 1
    
    % Calcula a fronteira de pareto do indiv�duo.
    c(L.COLPT) = indpar(c,arqsq,L);

    if nh == 1
        % c em refer�ncia isolada: deve ser inserido
        if size(arqsq,1) < maxpontos
            % Arquivo tem espa�o: insere c.
            arqsq = [arqsq; c];
        else
            % Arquivo est� cheio: substitui o pior.
            % Encontra o box do indiv�duo
            hy = indbox(c,inf,sqbox,boxsize,L);
            % Indiv�duos no mesmo hyperbox (�ndices)
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
            % Compara com o indiv�duo
            [m,~] = compdist(1,2,[c;cpr],L);
            if m == 1
                % c � melhor: substitui
                arqsq(pr,:) = c;
            else
                % c � pior: descarta
                descartado = 1;
            end
        end
    else
        % c em nicho ocupado
        % Obt�m o outro indiv�duo do nicho.
        inref = arqsq(rfinds,:);
        % Compara as dist�ncias
        if c(L.COLDT) < inref(L.COLDT)
            % c � melhor: substitui
            arqsq(rfinds,:) = c;
        else
            % c � pior: descarta
            descartado = 1;
        end
    end
    
end

if ~descartado
    % c foi inserido na popula��o: refaz os c�lculos de domina��o
    % e distribui��o.
    if ~extrap
        % Observa��o: este c�digo pode ser otimizado para realizar as
        % atualiza��es estritamente necess�rias, tal como foi feito
        % no caso da popula��o geral.
        no = max(L.COLF) - min(L.COLF) + 1;
        arqsq = pareto(arqsq,L);
        [arqsq,sqbox] = hyperbox(arqsq,inf,boxsize,L);
        [arqsq,sqref] = calcRefPontos(inf,sup,npontos,no,arqsq,L);
        arqsq = selecionaDistribuicao(arqsq,L);
    end
else
    extrap = 0; % N�o considera extrapola��o para ponto descartado
end
    
end

function m = indice(c,pop,L)
%INDICE Verifica se um indiv�duo est� em uma popula��o.
%   Verifica se um indiv�duo participa de uma popula��o.
%   Considera apenas os campos das vari�veis.
%
%   Par�metros de entrada:
%     - c: indiv�duo a considerar;
%     - pop: popula��o;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - m: 0 (n�o est� na popula��o), �ndice (caso contr�rio).

m = 0;
for i = 1:size(pop,1)
    if c(L.COLX) == pop(i,L.COLX)
        m = i;
        break;
    end
end

end

function dom = dominacao(c,pop,L)
%DOMINACAO Calcula elementos dominados por e que dominam um indiv�duo.
%   Encontra os indiv�duos de uma popula��o que dominam um determinado
%   indiv�duo e que s�o dominados por ele.
%
%   Par�metros de entrada:
%     - c: indiv�duo a comparar;
%     - pop: popula��o;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - dom: estrutura de domina��o, com os campos:
%         n: n�mero de indiv�duos dominados por c;
%         idn: vetor com �ndices dos indiv�duos dominados por c;
%         m: n�mero de indiv�duos que dominam c;
%         idm: vetor com �ndices dos indiv�duos que dominam c.

no = max(L.COLF) - min(L.COLF) + 1;
npop = size(pop,1);
elemento = struct('n',0, 'idn',zeros(1,npop), 'm',0, 'idm',zeros(1,npop));
dom = repmat(elemento,1,1);

for i = 1:npop
    % Compara com o i-�simo indiv�duo da popula��o
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
%EXTRAPOLA Verifica se um indiv�duo extrapola os limites do grid.
%   Verifica se um indiv�duo c extrapola os limites do grid de uma
%   determinada popula��o (e, portanto, aumenta sua diversidade).
%
%   Par�metros de entrada:
%     - c: indiv�duo a verificar;
%     - inf: limites inferiores do grid (para todas as dimens�es);
%     - sup: limites superiores do grid (para todas as dimens�es);
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - ext: 0 (n�o viola) ou 1 (viola).

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
%INDPAR Calcula a fronteira de Pareto de um indiv�duo.
%   Encontra a qual fronteira de Pareto um indiv�duo pertence
%   em rela��o a uma popula��o.
%
%   Par�metros de entrada:
%     - c: indiv�duo a verificar;
%     - pop: popula��o;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - pt: fronteira de Pareto do indiv�duo.

dom = dominacao(c,pop,L);

if dom.m == 0
    % N�o dominado: primeira fronteira
    pt = 1;
else
    % Dominado: verifica a �ltima fronteira dos dominadores e soma 1.
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
%INDREF Calcula o ponto de refer�ncia mais pr�ximo de um indiv�duo.
%   XXX
%
%   Par�metros de entrada:
%     - c: indiv�duo a verificar;
%     - refs: pontos de refer�ncia da popula��o
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - rf: ponto de refer�ncia do indiv�duo
  
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
%INDBOX Calcula o hyperbox de um indiv�duo.
%   Encontra o hyperbox ocupado por um indiv�duo em uma determinada
%   popula��o.
%
%   Par�metros de entrada:
%     - c: indiv�duo a verificar;
%     - inf: limites inferiores do grid (para todas as dimens�es);
%     - sup: limites superiores do grid (para todas as dimens�es);
%     - boxes: boxes da popula��o
%     - resolucao: tamanho do grid (o mesmo para todas as dimens�es);
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - hb: hyperbox do indiv�duo
  
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
%INSEREBOX Insere um indiv�duo em um ponto de refer�ncia.
%   Insere um indiv�duo em um ponto de refer�ncia e atualiza os dados
%   do ponto e dos indiv�duos ligados a ele.
%
%   Par�metros de entrada:
%     - pop: popula��o;
%     - ind: �ndice do indiv�duo na popula��o;
%     - ref: pontos de refer�ncia da popula��o
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - pop: XXX
%     - ref: pontos de refer�ncia atualizados.

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
%REMOVEBOX Remove um indiv�duo de um ponto de refer�ncia.
%   Remove um indiv�duo de um ponto de refer�ncia e atualiza os dados
%   do ponto.
%
%   Par�metros de entrada:
%     - pop: popula��o;
%     - ind: �ndice do indiv�duo na popula��o;
%     - ref: pontos de refer�ncia da popula��o;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - pop: XXX
%     - ref: pontos de refer�ncia atualizados.

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
%INSEREBOX Insere um indiv�duo em um box.
%   Insere um indiv�duo em um box e atualiza os dados do box.
%
%   Par�metros de entrada:
%     - pop: popula��o;
%     - ind: �ndice do indiv�duo na popula��o;
%     - inf: XXX
%     - box: boxes da popula��o
%     - boxsize: XXX
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
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
%REMOVEBOX Remove um indiv�duo de uma box.
%   Remove um indiv�duo de um box e atualiza os dados do box.
%
%   Par�metros de entrada:
%     - pop: popula��o;
%     - ind: �ndice do indiv�duo na popula��o;
%     - box: boxes da popula��o
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
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
%   Dada uma precisao, encontra a resolu��o do grid de
%   hyperboxes para uma determinada popula��o.
%
%   Par�metros de entrada:
%     - pop: popula��o;
%     - precisao: precis�o desejada;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - res: vetor com a resolu��o para cada dimens�o.

% Encontra os valores m�ximos e m�nimos das fun��es objetivo.
zmax = max(pop(:,L.COLF),[],1);
zmin = min(pop(:,L.COLF),[],1);

res = floor((zmax - zmin) / precisao) + 1;

end

function ind = pior(pop,inds,comp,L)
%PIOR Encontra o pior indiv�duo de uma (sub)popula��o.
%   Encontra o pior indiv�duo dentre um conjunto de indiv�duos
%   (cujos �ndices s�o fornecidos) de uma popula��o.
%   Utiliza uma fun��o fornecida para fazer as compara��es.
%
%   Par�metros de entrada:
%     - pop: popula��o;
%     - inds: �ndices dos indiv�duos a considerar;
%     - comp: fun��o de compara��o.
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - ind: �ndice do pior indiv�duo.

n = size(inds,2);

ind = inds(1);
for k = 2:n
    [~,ind] = comp(ind,inds(k),pop,L);
end

end

function pop = remove(pop,inds,n,L)
%REMOVE remove indiv�duos de uma popula��o.
%   Remove de uma popula��o os indiv�duos cujos �ndices s�o fornecidos.
%   Observa��o: pressup�e que os �ndices est�o em ordem crescente.
%
%   Par�metros de entrada:
%     - pop: popula��o;
%     - inds: �ndices dos indiv�duos a remover;
%     - n: quantidade de elementos a remover;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - pop: popula��o atualizada.

% Aloca espa�o para a nova popula��o
npop = size(pop,1);
rpop = zeros(npop-n,L.NC);

% Efetua as remo��es
inicio = 1;
c = 1;
for i = 1:n
    % Copia elementos anteriores ao �ndice
    for k = inicio:inds(i)-1
        rpop(c,:) = pop(k,:);
        c = c+1;
    end
    
    % Novo in�cio � ap�s o elemento removido
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
%CRUZAMENTO Realiza o cruzamento entre dois indiv�duos.
%   Realiza o cruzamento entre dois indiv�duos e retorna o filho
%   resultante.
%
%   Par�metros de entrada:
%     - ind1: primeiro indiv�duo;
%     - ind2: segundo indiv�duo;
%     - pc: probabilidade de cruzamento;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - f: filho gerado pelo cruzamento.

% Determina o melhor e o pior pai.
pais = [ind1; ind2];
[m,p] = compara(1,2,pais,L);
mp = pais(m,:);
pp = pais(p,:);

% Verifica se o cruzamento deve ser realizado;
% caso contr�rio, apenas copia o melhor pai.
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
%PERTURBACAO Calcula o vetor de perturba��es.
%   O vetor de perturba��es corresponde aos valores a serem
%   adicionados �s vari�veis do indiv�duo para a implementa��o
%   da muta��o.
%
%   Par�metros de entrada:
%     - c: indiv�duo a perturbar;
%     - pm: probabilidade de muta��o;
%     - xmin: valor m�nimo de uma vari�vel;
%     - xmax: valor m�ximo de uma vari�vel;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - gamma: vetor de perturba��es.

% Inicializa o vetor de perturba��es.
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
%CALCPROBS Calcula as probabilidades de cruzamento e muta��o.
%   Calcula as probabilidades de cruzamento e muta��o da presente gera��o.
%   Essas probabilidades s�o influenciadas pelo mdg.
%
%   Par�metros de entrada:
%     - pop: array com os indiv�duos da popula��o;
%     - pc: probabilidade de cruzamento atual;
%     - pm: probabilidade de muta��o atual;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - pc: nova probabilidade de cruzamento;
%     - pm: nova probabilidade de muta��o.


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
%CALCMDG Calcula a medida da diversidade gen�tica.
%   A medida da diversidade gen�tica da popula��o influencia as
%   probabilidades de cruzamento e muta��o. O c�lculo da medida
%   baseia-se no valor agregado das fun��es objetivo.
%
%   Par�metros de entrada:
%     - pop: array com os indiv�duos da popula��o;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%       - mdg: medida de diversidade gen�tica da popula��o.

npop = size(pop,1);
fit = pop(:,L.COLAG);
fmed = sum(fit) / npop;
fmax = max(fit);

mdg = (fmax - fmed) / fmax;

end

function pop = calcFos(pop,objFunc,L)
%CALCFOS Calcula as fun��es objetivo de uma popula��o.
%   Utiliza a fun��o fornecida como par�metro para calcular
%   os valores para cada indiv�duo.
%
%   Par�metros de entrada:
%     - pop: array com os indiv�duos da popula��o;
%     - objfunc: fun��o objetivo a utilizar;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - pop: popula��o atualizada.

npop = size(pop,1);
%nvar = max(L.COLX) - min(L.COLX) + 1;
no = max(L.COLF) - min(L.COLF) + 1;

for i = 1:npop
    pop(i,L.COLF) = objFunc(pop(i,L.COLX),no);
end

end

%---------------------

function f = dtlz1(x,m)

% x = indiv�duo
% m = n�mero de fun��es objetivo

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

% x = indiv�duo
% m = n�mero de fun��es objetivo

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
