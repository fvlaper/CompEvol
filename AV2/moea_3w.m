function [ ps ] = moea_3w( ncal, nvar, no )
%MOEA_3W Multi Objective Evolutionary Algorithm - 3 Way
%   Esta fun��o implementa um algoritmo para otimiza��o de fun��es
%   Multi-objetivo. XXX Explicar o algoritmo.
%
%   Par�metros de entrada:
%     - ncal: n�mero m�ximo de c�lculos das fun��es objetivo;
%     - nvar: n�mero de vari�veis;
%     - no: n�mero de fun��es objetivo.
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
%   - Fator de aglomera��o (squeeze factor): uma coluna.
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

% Dados iniciais da popula��o: faixa das vari�veis e 
% quantidade de indiv�duos.
xmin = 0; xmax = 1;
npop = 10; % XXX C�lculo?

% Estabelecimento da polula��o inicial.
pop = popinit(npop,xmin,xmax,L);

% Avalia��o dos elementos da popula��o inicial.
%fos = @dtlz1; % XXX Par�metro?
fos = @dtlz2;
pop = fos(pop,nvar,no,L);
ncal = ncal - npop;

pop = agregacao(pop,L); % c�lculo do valor agregado
pop = pareto(pop,L);    % c�lculo da fronteira de pareto
resolucao = 5; % XXX Par�metro?
[pop, pinf, psup] = hyperbox(pop,resolucao,L); % c�lculo do hyperbox

% C�pia dos elementos n�o dominados da popula��o para o arquivo arqnd
arqnd = naodominado(pop,L);
ndinf = pinf;
ndsup = psup;

% C�pia dos elementos n�o agregados da popula��o para o arquivo arqsq
arqsq = naoaglomerado(pop,L);
sqinf = pinf;
sqsup = psup;

display(pop);

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
        c = fos(c,nvar,no,L);
        c = agregacao(c,L);
        
        % Teste
        display(c);
        [pop,pinf,psup] = atualizapop(c,pop,pinf,psup,resolucao,L);
        display(pop);
    end
    
    % Foram avaliadas as fun��es objetivo de tr�s indiv�duos
    ncal = ncal - 3;
end

%for i = 1:npop
%    display(decodifica(pop(i,L.COLHY),no,resolucao));
%end

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
%   Para cada indiv�duo da polula��o, calcula a que fronteira de
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
    
    % Trata os indiv�duos na primeira fronteira.
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
    % Percorre os elementos da f-�sima fronteira, removendo-os
    % e fazendo os c�lculos para a (f+1)-�sima fronteira.
    for i = 1:front(f).n
        individuo = front(f).ids(i);
%        display(['Tratando indiv�duo ', num2str(individuo)]);
        % Trata os indiv�duos dominados pelo indiv�duo atual.
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
%HYPERBOX Calcula o hyperbox de cada indiv�duo.
%   Para cada indiv�duo da popula��o, calcula o hyperbox ao qual
%   pertence e o fator de aglomera��o. Retorna tamb�m os limites
%   inferior e superior do grid para cada dimens�o.
%
%   Par�metros de entrada:
%     - pop: array contendo a popula��o inicial;
%     - resolucao: tamanho do grid (o mesmo para todas as dimens�es).
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - pop: popula��o inicial hyperbox e aglomera��o atualizados;
%     - inf: limites inferiores do grid para cada dimens�o;
%     - sup: limites superiores do grid para cada dimens�o.

% Encontra os valores m�ximos e m�nimos das fun��es objetivo.
zmax = max(pop(:,L.COLF),[],1);
zmin = min(pop(:,L.COLF),[],1);

% Encontra os limites superiores e inferiores e o tamanho da caixa
delta = (zmax - zmin) / (2 * resolucao);
inf = zmin - delta;
sup = zmax + delta;
boxsize = (sup - inf) / resolucao;

% Calcula o hyperbox (para cada dimens�o) e o fator de aglomeracao
npop = size(pop,1); % n�mero de indiv�duos
no = max(L.COLF) - min(L.COLF) + 1; % n�mero de fun��es objetivo
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

% Registra o valor de aglomera��o para cada indiv�duo.
for k = 1:p
%    display(pbox(k));
    for j = 1:pbox(k).n
        pop(pbox(k).ids(j),L.COLSQ) = pbox(k).n;
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
%   pareto > fator de agrega��o > valor agregado.
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
%   fator de agrega��o > pareto > valor agregado.
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
%CODIFICA Codifica um hyperbox como um valor num�rico.
%   Codifica as coordenadas do hyperbox em um valor num�rico �nico.
%
%   Par�metros de entrada:
%     - hbox: hyperbox a codificar;
%     - resolucao: resolucao utilizada para o c�lculo do hyperbox.
%
%   Par�metros de sa�da:
%     - box: hyperbox codificado como inteiro.

no = size(hbox,2);
box = 0;
for i = 1:no
    box = (box * resolucao) + hbox(i);
end

end

function hbox = decodifica(box,no,resolucao)
%DECODIFICA Decodifica um hyperbox como um vetor de dimens�es.
%   Obt�m o vetor de dimens�es de um hyperbox a aprtir do valor
%   codificado.
%
%   Par�metros de entrada:
%     - box: hyperbox a decodificar;
%     - no: n�mero de fun��es objetivo;
%     - resolucao: resolucao utilizada para o c�lculo do hyperbox.
%
%   Par�metros de sa�da:
%     - box: hyperbox decodificado.

hbox = zeros(1,no);
for i = 0:no-1
    hbox(no-i) = rem(box,resolucao);
    box = fix(box/resolucao);
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

function arqsq = naoaglomerado(pop,L)
%NAOAGLOMERADO Encontra os elementos n�o aglomerados de uma populacao.
%   Encontra os elementos n�o aglomerados de uma popula��o (ou seja,
%   aqueles sozinhos em um hyperbox) e retorna-os
%   em uma subpopula��o separada.
%
%   Par�metros de entrada:
%     - pop: array contendo a popula��o inicial;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - arqnd: subpopula��o de elementos n�o aglomerados.

idx = (pop(:,L.COLSQ) == 1);
arqsq = pop(idx,:);

end

function [pop,pinf,psup] = atualizapop(c,pop,pinf,psup,resolucao,L)
%ATUALIZAPOP Atualiza a popula��o com um indiv�duo.
%   Tenta inserir um indiv�duo em uma popula��o. Se o indiv�duo
%   n�o for descartado, refaz os c�lculos das fronteiras de Pareto
%   e do hyperbox.
%
%   Par�metros de entrada:
%     - c: indiv�duo a inserir;
%     - pop: popula��o;
%     - inf: limites inferiores do grid para cada dimens�o;
%     - sup: limites superiores do grid para cada dimens�o;
%     - resolucao: tamanho do grid.
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - pop: popula��o atualizada.
%     - inf: limites inferiores atualizados;
%     - sup: limites superiores atualizados.

% Calcula os indiv�duos de pop dominados e que dominam c.
dom = dominacao(c,pop,L);

descartado = 0;

if dom.n > 0
    % c domina um grupo de indiv�duos:
    % substitui o pior deles por c.
    inds = dom.idn(1:dom.n);
    pr = pior(pop,inds,@compara,L);
    pop(pr,:) = c;
elseif dom.m > 0
    % c � dominado e descartado.
    descartado = 1;
else
    % c n�o domina nem � dominado:
    % substitui o pior indiv�duo da popula��o.
    inds = 1:size(pop,1);
    pr = pior(pop,inds,@compara,L);
    pop(pr,:) = c;
end

if ~descartado
    % c foi inserido na popula��o: refaz os c�lculos de domina��o
    % e distribui��o.
    pop = pareto(pop,L);
    [pop, pinf, psup] = hyperbox(pop,resolucao,L); % c�lculo do hyperbox
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
%       - pop: array com os indiv�duos da popula��o;
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

function pop = dtlz1 (pop,nvar,no,L)
%DTLZ1 Fun��o dtlz1.
%   Calcula as fun��es objetivo DTLZ1 para todos os indiv�duos
%   da popula��o.
%
%   Par�metros de entrada:
%     - pop: popula��o;
%     - nvar: n�mero de vari�veis;
%     - no: n�mero de fun��es objetivo;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - pop: popula��o com os valores das fun��es objetivo.

if no > nvar
    ME = MException('dtlz1:nvarError','Invalid number of variables');
    throw(ME);
end

% Colunas das vari�veis.
xi = min(L.COLX);
fi = min(L.COLF);
ff = max(L.COLF);
k = nvar - no + 1;

xg = pop(:,no:nvar);

% C�lculo
gx = 100 * (k + sum(((xg-0.5) .^ 2) - cos(20*pi*(xg-0.5)),2));

pop(:,fi) = 0.5 * (prod(pop(:,xi:(no-1)),2) .* (1+gx));

for i = 2:(no-1)
  pop(:,(fi+i-1)) = 0.5 * (prod((pop(:,xi:(no-i))),2) .* ...
                           (1 - pop(:,(no-i+1))) .* (1+gx));
end

pop(:,ff) = 0.5 * ((1 - pop(:,xi)) .* (1+gx));

end

function pop = dtlz2 (pop,nvar,no,L)
%DTLZ1 Fun��o dtlz1.
%   Calcula as fun��es objetivo DTLZ2 para todos os indiv�duos
%   da popula��o.
%
%   Par�metros de entrada:
%     - pop: popula��o;
%     - nvar: n�mero de vari�veis;
%     - no: n�emro de fun��es objetivo;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - pop: popula��o com os valores das fun��es objetivo.

if no > nvar
    ME = MException('dtlz2:nvarError','Invalid number of variables');
    throw(ME);
end

% Colunas das vari�veis.
xi = min(L.COLX);
fi = min(L.COLF);
ff = max(L.COLF);

xg = pop(:,no:nvar);

% C�lculo
gx = sum(((xg-0.5) .^ 2),2);

pop(:,fi) = (1+gx) .* prod(cos(0.5*pi*pop(:,xi:(no-1))),2);

for i = 2:(no-1)
    pop(:,(fi+i-1)) = (1+gx) .* prod(cos(0.5*pi*pop(:,xi:(no-i))),2) .* ...
                                sin(0.5*pi*pop(:,(no-i+1)));
end

pop(:,ff) = (1+gx) .* sin(0.5*pi*pop(:,xi));

end
