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
%   - Fator de nicho: uma coluna.
xi = 1; xf = xi+nvar-1;
fi = xf+1; ff = fi+no-1;
ag = ff+1;
pt = ag+1;
hy = pt+1;
sq = hy+1;
rf = sq+1;
nh = rf+1;
nc = nh;

L.COLX  = xi:xf;
L.COLF  = fi:ff;
L.COLAG = ag;
L.COLPT = pt;
L.COLHY = hy;
L.COLSQ = sq;
L.COLRF = rf;
L.COLNH = nh;
L.NC = nc;

% Dados iniciais da popula��o: faixa das vari�veis e 
% quantidade de indiv�duos.
xmin = 0; xmax = 1;
%npop = 10; % XXX C�lculo?

% Estabelecimento da popula��o inicial.
pop = popinit(npop,xmin,xmax,L);

% Avalia��o dos elementos da popula��o inicial.
%fos = @dtlz1; % XXX Par�metro?
%fos = @dtlz2;
%pop = fos(pop,nvar,no,L);
pop = calcFos(pop,objfunc,L);
ncal = ncal - npop;

pop = agregacao(pop,L); % c�lculo do valor agregado
pop = pareto(pop,L);    % c�lculo da fronteira de pareto
%resolucao = 5; % XXX Par�metro?
%[pop, pinf, psup] = hyperbox(pop,resolucao,L); % c�lculo do hyperbox
%[pop, pinf, psup, pbox] = hyperbox2(pop,resolucao,L); % c�lculo do hyperbox
[pop, pinf, psup, pbox,pres] = hyperbox3(pop,precisao,L); % c�lculo do hyperbox
%for k = 1:npop
%    display(pbox(k));
%end
[pop, pref] = calcRefPontos(pinf,psup,npop,no,pop,L); % c�lculo dos pontos de refer�ncia
%for mi = 1:length(refpontos)
%    display(refpontos(mi));
%end

% C�pia dos elementos n�o dominados da popula��o para o arquivo arqnd
arqnd = naodominado(pop,L);
ndbox = pbox;
ndinf = pinf;
ndsup = psup;
ndres = pres;
%ndref = pref;
%[arqnd, ndref] = calcRefPontos(ndinf,ndsup,npop,no,arqnd,L);

% C�pia dos elementos n�o agregados da popula��o para o arquivo arqsq
%arqsq = naoaglomerado(pop,L);
%arqsq = selecionaDistribuicao(pop,L);
arqsq = selecionaDistribuicao2(pop,L);
sqbox = pbox;
sqinf = pinf;
sqsup = psup;
sqres = pres;
%sqref = pref;
[arqsq, sqref] = calcRefPontos(sqinf,sqsup,npop,no,arqsq,L);
%for k = 1:npop
%    display(sqbox(k).hbox);
%end

% Probabilidades iniciais de cruzamento e muta��o
pc = 0.6;
pm = 0.05;

%display(arqsq);

% La�o principal
while ncal >= 3
    
    % Seleciona os pais para os cruzamentos (um de cada popula�ao/arquivo)
    %p1 = seleciona(pop,@compara,L);
    %p2 = seleciona(arqnd,@compara,L);
    p1 = seleciona(pop,@compara2,L);
    p2 = seleciona(arqnd,@compara2,L);
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
        display(c);
%{
        display(pop);
        %[pop,pinf,psup,] = atualizapop(c,pop,pinf,psup,resolucao,L);
        [pop,pinf,psup,pbox,pres] = atualizapop3(c,pop,pinf,psup,pbox,pres,precisao,L);
        display(pop)
        display(arqnd);
        [arqnd,ndinf,ndsup,ndbox,ndres] = atualizadom3(c,arqnd,ndinf,ndsup,ndbox,ndres,...
                                          npop,precisao,L);
        display(arqnd);
        
        display(arqsq);
        %[arqsq,sqinf,sqsup] = atualizadis(c,arqsq,sqinf,sqsup,...
        %                                  npop,resolucao,L);
        %[arqsq,sqinf,sqsup,sqbox] = atualizadis2(c,arqsq,sqinf,sqsup,sqbox,...
        %                                  npop,resolucao,L);
        [arqsq,sqinf,sqsup,sqbox,sqres] = atualizadis3(c,arqsq,sqinf,sqsup,sqbox,sqres,...
                                          npop,precisao,L);
        display(arqsq);                                              
%}
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
    %display(pbox(k));
    for j = 1:pbox(k).n
        pop(pbox(k).ids(j),L.COLSQ) = pbox(k).n;
    end
end

end

function [pop, inf, sup, pbox] = hyperbox2(pop,resolucao,L)
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
%     - sup: limites superiores do grid para cada dimens�o;
%     - pbox: vetor de estruturas com as informa��es dos boxes:
%             hbox: coordenadas do box;
%             n: n�mero de elementos do box;
%             ids: �ndices dos indiv�duos do box.

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
%display(hbox);

% Registra o valor de aglomera��o para cada indiv�duo.
for k = 1:p
    %display(pbox(k));
    for j = 1:pbox(k).n
        pop(pbox(k).ids(j),L.COLSQ) = pbox(k).n;
    end
end

end

function [pop, inf, sup, pbox, resolucao] = hyperbox3(pop,precisao,L)
%HYPERBOX Calcula o hyperbox de cada indiv�duo.
%   Para cada indiv�duo da popula��o, calcula o hyperbox ao qual
%   pertence e o fator de aglomera��o. Retorna tamb�m os limites
%   inferior e superior do grid para cada dimens�o.
%
%   Par�metros de entrada:
%     - pop: array contendo a popula��o inicial;
%     - precisao: XXX
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - pop: popula��o inicial hyperbox e aglomera��o atualizados;
%     - inf: limites inferiores do grid para cada dimens�o;
%     - sup: limites superiores do grid para cada dimens�o;
%     - pbox: vetor de estruturas com as informa��es dos boxes:
%             hbox: coordenadas do box;
%             n: n�mero de elementos do box;
%             ids: �ndices dos indiv�duos do box;
%     - resolucao: XXX

% Encontra os valores m�ximos e m�nimos das fun��es objetivo.
zmax = max(pop(:,L.COLF),[],1);
zmin = min(pop(:,L.COLF),[],1);

resolucao = calcResolucao(pop,precisao,L);

% Encontra os limites superiores e inferiores e o tamanho da caixa
delta = (zmax - zmin) ./ (2 * resolucao);
inf = zmin - delta;
sup = zmax + delta;
boxsize = (sup - inf) ./ resolucao;

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
%display(hbox);

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
%refpontos = zeros(npontos,no);

% Preenche os pontos de refer�ncia
ponto = inf; % Come�a pelo "canto inferior esquerdo"
for k = 1:npontos
    %refpontos(k,:) = ponto; % grava o ponto
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

function [melhor, pior] = compara2(i1, i2, pop, L)
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

function [melhor, pior] = compdist2(i1, i2, pop, L)
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

function arqsq = selecionaDistribuicao(pop,L)
%SelecionaDistribuicao Seleciona os elementos do arquivo de distribui��o.
%   Examina os elementos de uma popula��o e seleciona os melhores de
%   cada box para compor o arquivo de destribui��o.
%
%   Par�metros de entrada:
%     - pop: array contendo a popula��o inicial;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - arqsq: subpopula��o de elementos n�o aglomerados.

popsort = sortrows(pop,[L.COLHY L.COLPT L.COLAG]);

numboxes = length(unique(popsort(:,L.COLHY)));
arqsq = zeros(numboxes,L.NC);

lastbox = -1;
ind = 1;

for i = 1:size(popsort,1)
    if popsort(i,L.COLHY) ~= lastbox
        arqsq(ind,:) = popsort(i,:);
        arqsq(ind,L.COLSQ) = 1;
        ind = ind + 1;
        lastbox = popsort(i,L.COLHY);
    end
end

end

function arqsq = selecionaDistribuicao2(pop,L)
%SelecionaDistribuicao Seleciona os elementos do arquivo de distribui��o.
%   Examina os elementos de uma popula��o e seleciona os melhores de
%   cada nicho para compor o arquivo de destribui��o.
%
%   Par�metros de entrada:
%     - pop: array contendo a popula��o inicial;
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - arqsq: subpopula��o de elementos n�o aglomerados.

popsort = sortrows(pop,[L.COLRF L.COLHY L.COLPT L.COLAG]);

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

display('Introdu��o na popula��o');

% Verifica se indiv�duo n�o existe na popula��o
%idx = indice(c,pop,L);
%if idx > 0
%    display('c j� existe: descartado');
%    return;
%end

% Calcula os indiv�duos de pop dominados e que dominam c.
dom = dominacao(c,pop,L);

descartado = 0;

if dom.n > 0
    % c domina um grupo de indiv�duos:
    % substitui o pior deles por c.
    inds = dom.idn(1:dom.n);
    pr = pior(pop,inds,@compara,L);
    pop(pr,:) = c;
    display(['c domina: substitui ', num2str(pr)]);
elseif dom.m > 0
    % c � dominado e descartado.
    descartado = 1;
    display('c dominado: descartado');
else
    % c n�o domina nem � dominado:
    % substitui o pior indiv�duo da popula��o.
    inds = 1:size(pop,1);
    pr = pior(pop,inds,@compara,L);
    pop(pr,:) = c;
    display(['c n�o domina nem � dominado: substitui ', num2str(pr)]);
end

if ~descartado
    % c foi inserido na popula��o: refaz os c�lculos de domina��o
    % e distribui��o.
    display('c inserido (recalcula)');
    pop = pareto(pop,L);
    [pop, pinf, psup] = hyperbox(pop,resolucao,L); % c�lculo do hyperbox
end

end

function [pop,pinf,psup,pbox] = atualizapop2(c,pop,pinf,psup,pbox,resolucao,L)
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
%     - pbox: box da popula��o;
%     - resolucao: tamanho do grid.
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - pop: popula��o atualizada.
%     - inf: limites inferiores atualizados;
%     - sup: limites superiores atualizados.

display('Introdu��o na popula��o');

% Verifica se indiv�duo n�o existe na popula��o
%idx = indice(c,pop,L);
%if idx > 0
%    display('c j� existe: descartado');
%    return;
%end

% Calcula os indiv�duos de pop dominados e que dominam c.
dom = dominacao(c,pop,L);

descartado = 0;

if dom.n > 0
    % c domina um grupo de indiv�duos:
    % substitui o pior deles por c.
    inds = dom.idn(1:dom.n);
    pr = pior(pop,inds,@compara,L);
    pop(pr,:) = c;
    display(['c domina: substitui ', num2str(pr)]);
elseif dom.m > 0
    % c � dominado e descartado.
    descartado = 1;
    display('c dominado: descartado');
else
    % c n�o domina nem � dominado:
    % substitui o pior indiv�duo da popula��o.
    inds = 1:size(pop,1);
    pr = pior(pop,inds,@compara,L);
    pop(pr,:) = c;
    display(['c n�o domina nem � dominado: substitui ', num2str(pr)]);
end

if ~descartado
    % c foi inserido na popula��o: refaz os c�lculos de domina��o
    % e distribui��o.
    display('c inserido (recalcula)');
    pop = pareto(pop,L);
    [pop, pinf, psup, pbox] = hyperbox2(pop,resolucao,L); % c�lculo do hyperbox
end

end

function [pop,pinf,psup,pbox,pres] = atualizapop3(c,pop,pinf,psup,pbox,pres,precisao,L)
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
%     - pbox: box da popula��o;
%     - pres: XXX
%     - precisao: XXX
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - pop: popula��o atualizada.
%     - inf: limites inferiores atualizados;
%     - sup: limites superiores atualizados.
%     - pbox: XXX
%     - pres: XXX

display('Introdu��o na popula��o');

% Verifica se indiv�duo n�o existe na popula��o
%idx = indice(c,pop,L);
%if idx > 0
%    display('c j� existe: descartado');
%    return;
%end

% Calcula os indiv�duos de pop dominados e que dominam c.
dom = dominacao(c,pop,L);

descartado = 0;

if dom.n > 0
    % c domina um grupo de indiv�duos:
    % substitui o pior deles por c.
    inds = dom.idn(1:dom.n);
    pr = pior(pop,inds,@compara,L);
    pop(pr,:) = c;
    display(['c domina: substitui ', num2str(pr)]);
elseif dom.m > 0
    % c � dominado e descartado.
    descartado = 1;
    display('c dominado: descartado');
else
    % c n�o domina nem � dominado:
    % substitui o pior indiv�duo da popula��o.
    inds = 1:size(pop,1);
    pr = pior(pop,inds,@compara,L);
    pop(pr,:) = c;
    display(['c n�o domina nem � dominado: substitui ', num2str(pr)]);
end

if ~descartado
    % c foi inserido na popula��o: refaz os c�lculos de domina��o
    % e distribui��o.
    display('c inserido (recalcula)');
    pop = pareto(pop,L);
    [pop, pinf, psup, pbox, pres] = hyperbox3(pop,precisao,L); % c�lculo do hyperbox
end

end

function [arqnd,ndinf,ndsup] = atualizadom(c,arqnd,ndinf,ndsup, ...
                                           max,resolucao,L)
%ATUALIZAPOP Atualiza o arquivo de domina��o com um indiv�duo.
%   Tenta inserir um indiv�duo no arquivo de domina��o popula��o.
%   Se o indiv�duo n�o for descartado, refaz os c�lculos 
%   das fronteiras de Pareto e do hyperbox.
%
%   Par�metros de entrada:
%     - c: indiv�duo a inserir;
%     - arqnd: popula��o (arquivo de domina��o);
%     - ndinf: limites inferiores do grid para cada dimens�o;
%     - ndsup: limites superiores do grid para cada dimens�o;
%     - max: n�mero m�ximo de elementos em no arquivo;
%     - resolucao: tamanho do grid.
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - arqnd: popula��o atualizada.
%     - ndinf: limites inferiores atualizados;
%     - ndsup: limites superiores atualizados.

display('Introdu��o no arquivo de domina��o');

% Verifica se indiv�duo n�o existe na popula��o
idx = indice(c,arqnd,L);
if idx > 0
    display('c j� existe: descartado');
    return;
end

% Calcula os indiv�duos de pop dominados e que dominam c.
dom = dominacao(c,arqnd,L);

descartado = 0;

if dom.m > 0
    % c � dominado por algum indiv�duo do arquivo e � descartado
    display('c dominado: descartado');
    descartado = 1;
elseif dom.n > 0
    % c domina indiv�duos do arquivo: remove os dominados e insere c
    display('c domina: substitui dominados');
    display(dom.idn);
    arqnd = remove(arqnd,dom.idn,dom.n,L);
    arqnd = [arqnd; c];
elseif size(arqnd,1) < max
    % c n�o domina nem � dominado e h� espa�o no arquivo: insere c
    display('n�o domina nem � dominado; h� espa�o: insere');
    arqnd = [arqnd; c];
else
    % c n�o domina nem � dominado; n�o h� espa�o no arquivo
    display('n�o domina nem � dominado; n�o h� espa�o');
    % Verifica se c estende o grid do arquivo
    if extrapola(c,ndinf,ndsup,L)
        % Substitui o pior elemento do grid por c
        display('c extrapola o grid');
        pr = pior(arqnd,1:size(arqnd,1),@compara,L);
        display(['substitui ', num2str(pr)]);
        arqnd(pr,:) = c;
    else
        % Calcula o hyperbox e o squeeze factor de c
        display('c n�o extrapola o grid');
        hy = indbox(c,ndinf,ndsup,resolucao,L);
        % Indiv�duos no mesmo hyperbox (�ndices)
        hyinds = find(arqnd(:,L.COLHY) == hy);
        sq = size(hyinds,1) + 1;
        c(L.COLPT) = 1;
        c(L.COLHY) = hy;
        c(L.COLSQ) = sq;
        % Encontra os elementos com squeeze factor igual ou superior
        sqinds = transpose(find(arqnd(:,L.COLSQ) >= sq));
        if isempty(sqinds)
            sqsups = [];
        else
            sqsups = arqnd(sqinds,:);
        end
        % Encontra indiv�duos no mesmo hyperbox e atualiza squeeze factor
        if isempty(hyinds)
            hysups = [];
        else
            hysups = arqnd(hyinds,:);
            hysups(:,L.COLSQ) = sq;
        end
        % Une os dois conjuntos e determina o pior indiv�duo.
        sups = [sqsups; hysups];
        pr = sups(pior(sups,1:size(sups,1),@compara,L),:);
        % Compara o pior com c
        [m,~] = compara(1,2,[c;pr],L);
        if m == 1
            % c � melhor: substitui
            idx = indice(pr,arqnd,L);
            display(['substitui ', num2str(idx)]);
            arqnd(idx,:) = c;
        else
            % c � pior: descarta
            display('n�o substitui (descartado)');
            descartado = 1;
        end
    end
end

if ~descartado
    display('c inserido (recalcula)');
    % c foi inserido na popula��o: refaz os c�lculos de domina��o
    % e distribui��o.
    arqnd = pareto(arqnd,L);
    [arqnd, ndinf, ndsup] = hyperbox(arqnd,resolucao,L); % c�lculo do hyperbox
end
    
end

function [arqnd,ndinf,ndsup,ndbox] = atualizadom2(c,arqnd,ndinf,ndsup,ndbox, ...
                                           max,resolucao,L)
%ATUALIZAPOP Atualiza o arquivo de domina��o com um indiv�duo.
%   Tenta inserir um indiv�duo no arquivo de domina��o popula��o.
%   Se o indiv�duo n�o for descartado, refaz os c�lculos 
%   das fronteiras de Pareto e do hyperbox.
%
%   Par�metros de entrada:
%     - c: indiv�duo a inserir;
%     - arqnd: popula��o (arquivo de domina��o);
%     - ndinf: limites inferiores do grid para cada dimens�o;
%     - ndsup: limites superiores do grid para cada dimens�o;
%     - nxbox: boxes do arquivo;
%     - max: n�mero m�ximo de elementos em no arquivo;
%     - resolucao: tamanho do grid.
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - arqnd: popula��o atualizada.
%     - ndinf: limites inferiores atualizados;
%     - ndsup: limites superiores atualizados;
%     - ndbox: boxes atualizado.

display('Introdu��o no arquivo de domina��o');

% Verifica se indiv�duo n�o existe na popula��o
idx = indice(c,arqnd,L);
if idx > 0
    display('c j� existe: descartado');
    return;
end

% Calcula os indiv�duos de pop dominados e que dominam c.
dom = dominacao(c,arqnd,L);

descartado = 0;

if dom.m > 0
    % c � dominado por algum indiv�duo do arquivo e � descartado
    display('c dominado: descartado');
    descartado = 1;
elseif dom.n > 0
    % c domina indiv�duos do arquivo: remove os dominados e insere c
    display('c domina: substitui dominados');
    display(dom.idn);
    arqnd = remove(arqnd,dom.idn,dom.n,L);
    arqnd = [arqnd; c];
elseif size(arqnd,1) < max
    % c n�o domina nem � dominado e h� espa�o no arquivo: insere c
    display('n�o domina nem � dominado; h� espa�o: insere');
    arqnd = [arqnd; c];
else
    % c n�o domina nem � dominado; n�o h� espa�o no arquivo
    display('n�o domina nem � dominado; n�o h� espa�o');
    % Verifica se c estende o grid do arquivo
    if extrapola(c,ndinf,ndsup,L)
        % Substitui o pior elemento do grid por c
        display('c extrapola o grid');
        pr = pior(arqnd,1:size(arqnd,1),@compara,L);
        display(['substitui ', num2str(pr)]);
        arqnd(pr,:) = c;
    else
        % Calcula o hyperbox e o squeeze factor de c
        display('c n�o extrapola o grid');
        %hy = indbox(c,ndinf,ndsup,resolucao,L);
        hy = indbox2(c,ndinf,ndsup,ndbox,resolucao,L);
        % Indiv�duos no mesmo hyperbox (�ndices)
        hyinds = find(arqnd(:,L.COLHY) == hy);
        sq = size(hyinds,1) + 1;
        c(L.COLPT) = 1;
        c(L.COLHY) = hy;
        c(L.COLSQ) = sq;
        % Encontra os elementos com squeeze factor igual ou superior
        sqinds = transpose(find(arqnd(:,L.COLSQ) >= sq));
        if isempty(sqinds)
            sqsups = [];
        else
            sqsups = arqnd(sqinds,:);
        end
        % Encontra indiv�duos no mesmo hyperbox e atualiza squeeze factor
        if isempty(hyinds)
            hysups = [];
        else
            hysups = arqnd(hyinds,:);
            hysups(:,L.COLSQ) = sq;
        end
        % Une os dois conjuntos e determina o pior indiv�duo.
        sups = [sqsups; hysups];
        pr = sups(pior(sups,1:size(sups,1),@compara,L),:);
        % Compara o pior com c
        [m,~] = compara(1,2,[c;pr],L);
        if m == 1
            % c � melhor: substitui
            idx = indice(pr,arqnd,L);
            display(['substitui ', num2str(idx)]);
            arqnd(idx,:) = c;
        else
            % c � pior: descarta
            display('n�o substitui (descartado)');
            descartado = 1;
        end
    end
end

if ~descartado
    display('c inserido (recalcula)');
    % c foi inserido na popula��o: refaz os c�lculos de domina��o
    % e distribui��o.
    arqnd = pareto(arqnd,L);
    %[arqnd, ndinf, ndsup] = hyperbox(arqnd,resolucao,L); % c�lculo do hyperbox
    [arqnd, ndinf, ndsup, ndbox] = hyperbox2(arqnd,resolucao,L); % c�lculo do hyperbox
end
    
end

function [arqnd,ndinf,ndsup,ndbox,ndres] = atualizadom3(c,arqnd,ndinf,ndsup,ndbox,ndres, ...
                                           max,precisao,L)
%ATUALIZAPOP Atualiza o arquivo de domina��o com um indiv�duo.
%   Tenta inserir um indiv�duo no arquivo de domina��o popula��o.
%   Se o indiv�duo n�o for descartado, refaz os c�lculos 
%   das fronteiras de Pareto e do hyperbox.
%
%   Par�metros de entrada:
%     - c: indiv�duo a inserir;
%     - arqnd: popula��o (arquivo de domina��o);
%     - ndinf: limites inferiores do grid para cada dimens�o;
%     - ndsup: limites superiores do grid para cada dimens�o;
%     - nxbox: boxes do arquivo;
%     - ndres: XXX
%     - max: n�mero m�ximo de elementos em no arquivo;
%     - recisao: XXX
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - arqnd: popula��o atualizada.
%     - ndinf: limites inferiores atualizados;
%     - ndsup: limites superiores atualizados;
%     - ndbox: boxes atualizado;
%     - ndres: XXX

display('Introdu��o no arquivo de domina��o');

% Verifica se indiv�duo n�o existe na popula��o
idx = indice(c,arqnd,L);
if idx > 0
    display('c j� existe: descartado');
    return;
end

% Calcula os indiv�duos de pop dominados e que dominam c.
dom = dominacao(c,arqnd,L);

descartado = 0;

if dom.m > 0
    % c � dominado por algum indiv�duo do arquivo e � descartado
    display('c dominado: descartado');
    descartado = 1;
elseif dom.n > 0
    % c domina indiv�duos do arquivo: remove os dominados e insere c
    display('c domina: substitui dominados');
    display(dom.idn);
    arqnd = remove(arqnd,dom.idn,dom.n,L);
    arqnd = [arqnd; c];
elseif size(arqnd,1) < max
    % c n�o domina nem � dominado e h� espa�o no arquivo: insere c
    display('n�o domina nem � dominado; h� espa�o: insere');
    arqnd = [arqnd; c];
else
    % c n�o domina nem � dominado; n�o h� espa�o no arquivo
    display('n�o domina nem � dominado; n�o h� espa�o');
    % Verifica se c estende o grid do arquivo
    if extrapola(c,ndinf,ndsup,L)
        % Substitui o pior elemento do grid por c
        display('c extrapola o grid');
        pr = pior(arqnd,1:size(arqnd,1),@compara,L);
        display(['substitui ', num2str(pr)]);
        arqnd(pr,:) = c;
    else
        % Calcula o hyperbox e o squeeze factor de c
        display('c n�o extrapola o grid');
        %hy = indbox(c,ndinf,ndsup,resolucao,L);
        %hy = indbox2(c,ndinf,ndsup,ndbox,resolucao,L);
        hy = indbox3(c,ndinf,ndsup,ndbox,ndres,L);
        % Indiv�duos no mesmo hyperbox (�ndices)
        hyinds = find(arqnd(:,L.COLHY) == hy);
        sq = size(hyinds,1) + 1;
        c(L.COLPT) = 1;
        c(L.COLHY) = hy;
        c(L.COLSQ) = sq;
        % Encontra os elementos com squeeze factor igual ou superior
        sqinds = transpose(find(arqnd(:,L.COLSQ) >= sq));
        if isempty(sqinds)
            sqsups = [];
        else
            sqsups = arqnd(sqinds,:);
        end
        % Encontra indiv�duos no mesmo hyperbox e atualiza squeeze factor
        if isempty(hyinds)
            hysups = [];
        else
            hysups = arqnd(hyinds,:);
            hysups(:,L.COLSQ) = sq;
        end
        % Une os dois conjuntos e determina o pior indiv�duo.
        sups = [sqsups; hysups];
        pr = sups(pior(sups,1:size(sups,1),@compara,L),:);
        % Compara o pior com c
        [m,~] = compara(1,2,[c;pr],L);
        if m == 1
            % c � melhor: substitui
            idx = indice(pr,arqnd,L);
            display(['substitui ', num2str(idx)]);
            arqnd(idx,:) = c;
        else
            % c � pior: descarta
            display('n�o substitui (descartado)');
            descartado = 1;
        end
    end
end

if ~descartado
    display('c inserido (recalcula)');
    % c foi inserido na popula��o: refaz os c�lculos de domina��o
    % e distribui��o.
    arqnd = pareto(arqnd,L);
    %[arqnd, ndinf, ndsup] = hyperbox(arqnd,resolucao,L); % c�lculo do hyperbox
    %[arqnd, ndinf, ndsup, ndbox] = hyperbox2(arqnd,resolucao,L); % c�lculo do hyperbox
    [arqnd, ndinf, ndsup, ndbox,ndres] = hyperbox3(arqnd,precisao,L); % c�lculo do hyperbox
end
    
end

function [arqsq,sqinf,sqsup] = atualizadis(c,arqsq,sqinf,sqsup, ...
                                           max,resolucao,L)
%ATUALIZADIS Atualiza o arquivo de domina��o com um indiv�duo.
%   Tenta inserir um indiv�duo no arquivo de distribui��o.
%   Se o indiv�duo n�o for descartado, refaz os c�lculos 
%   das fronteiras de Pareto e do hyperbox.
%
%   Par�metros de entrada:
%     - c: indiv�duo a inserir;
%     - arqsq: popula��o (arquivo de distribui��o);
%     - sqinf: limites inferiores do grid para cada dimens�o;
%     - sqsup: limites superiores do grid para cada dimens�o;
%     - max: n�mero m�ximo de elementos no arquivo;
%     - resolucao: tamanho do grid.
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - arqsq: popula��o atualizada.
%     - sqinf: limites inferiores atualizados;
%     - sqsup: limites superiores atualizados.

display('Introdu��o no arquivo de distribui��o');

% Verifica se indiv�duo n�o existe na popula��o
idx = indice(c,arqsq,L);
if idx > 0
    display('c j� existe: descartado');
    return;
end

descartado = 0;

if extrapola(c,sqinf,sqsup,L)
    % c estende grid do arquivo: deve ser inserido.
    display('c extrapola o grid');
    if size(arqsq,1) < max
        % Arquivo tem espa�o: insere c.
        arqsq = [arqsq; c];
        display('arquivo tem espa�o: inserido');
    else
        % Arquivo est� cheio: substitui o pior.
        pr = pior(arqsq,1:size(arqsq,1),@compdist,L);
        arqsq(pr,:) = c;
        display(['arquivo cheio: substui ', num2str(pr)]);
    end
else
    % C n�o estende grid: verifica de deve ser inserido
    display('c n�o extrapola grid');
    
    % Calcula o hyperbox e o squeeze factor de c.
    hy = indbox(c,sqinf,sqsup,resolucao,L);
    hyinds = find(arqsq(:,L.COLHY) == hy);
    sq = size(hyinds,1) + 1;
    c(L.COLHY) = hy;
    c(L.COLSQ) = sq;
    
    % Calcula a fronteira de pareto do indiv�duo.
    c(L.COLPT) = indpar(c,arqsq,L);
    display(c);
    
    if sq == 1
        % c em box individual: deve ser inserido
        if size(arqsq,1) < max
            % Arquivo tem espa�o: insere c.
            arqsq = [arqsq; c];
            display('arquivo tem espa�o: inserido');
        else
            % Arquivo est� cheio: substitui o pior.
            pr = pior(arqsq,1:size(arqsq,1),@compdist,L);
            arqsq(pr,:) = c;
            display(['arquivo cheio: substui ', num2str(pr)]);
        end
    else
        % c em box ocupado
        display('c em box ocupado');
        % Obt�m o outro indiv�duo do box.
        inbox = arqsq(hyinds,:);
        inbox(L.COLSQ) = sq;
        display(inbox);
        % Compara com c
        [m,~] = compdist(1,2,[c;inbox],L);
        if m == 1
            % c � melhor: substitui
            idx = indice(inbox,arqsq,L);
            display('substitui ');
            display(idx);
            arqsq(idx,:) = c;
        else
            % c � pior: descarta
            display('n�o substitui (descartado)');
            descartado = 1;
        end
    end
end

if ~descartado
    display('c inserido (recalcula)');
    % c foi inserido na popula��o: refaz os c�lculos de domina��o
    % e distribui��o.
    arqsq = pareto(arqsq,L);
    [arqsq, sqinf, sqsup] = hyperbox(arqsq,resolucao,L);
    %arqsq = naoaglomerado(arqsq,L);
    arqsq = selecionaDistribuicao(arqsq,L);
end
    
end

function [arqsq,sqinf,sqsup,sqbox] = atualizadis2(c,arqsq,sqinf,sqsup,sqbox, ...
                                           max,resolucao,L)
%ATUALIZADIS Atualiza o arquivo de domina��o com um indiv�duo.
%   Tenta inserir um indiv�duo no arquivo de distribui��o.
%   Se o indiv�duo n�o for descartado, refaz os c�lculos 
%   das fronteiras de Pareto e do hyperbox.
%
%   Par�metros de entrada:
%     - c: indiv�duo a inserir;
%     - arqsq: popula��o (arquivo de distribui��o);
%     - sqinf: limites inferiores do grid para cada dimens�o;
%     - sqsup: limites superiores do grid para cada dimens�o;
%     - sqbox: boxes do arquivo;
%     - max: n�mero m�ximo de elementos no arquivo;
%     - resolucao: tamanho do grid.
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - arqsq: popula��o atualizada.
%     - sqinf: limites inferiores atualizados;
%     - sqsup: limites superiores atualizados;
%     - sqbox: boxes atualizados.

display('Introdu��o no arquivo de distribui��o');

% Verifica se indiv�duo n�o existe na popula��o
idx = indice(c,arqsq,L);
if idx > 0
    display('c j� existe: descartado');
    return;
end

descartado = 0;

if extrapola(c,sqinf,sqsup,L)
    % c estende grid do arquivo: deve ser inserido.
    display('c extrapola o grid');
    if size(arqsq,1) < max
        % Arquivo tem espa�o: insere c.
        arqsq = [arqsq; c];
        display('arquivo tem espa�o: inserido');
    else
        % Arquivo est� cheio: substitui o pior.
        pr = pior(arqsq,1:size(arqsq,1),@compdist,L);
        arqsq(pr,:) = c;
        display(['arquivo cheio: substui ', num2str(pr)]);
    end
else
    % C n�o estende grid: verifica de deve ser inserido
    display('c n�o extrapola grid');
    
    % Calcula o hyperbox e o squeeze factor de c.
    %hy = indbox(c,sqinf,sqsup,resolucao,L);
    hy = indbox2(c,sqinf,sqsup,sqbox,resolucao,L);
    hyinds = find(arqsq(:,L.COLHY) == hy);
    sq = size(hyinds,1) + 1;
    c(L.COLHY) = hy;
    c(L.COLSQ) = sq;
    
    % Calcula a fronteira de pareto do indiv�duo.
    c(L.COLPT) = indpar(c,arqsq,L);
    %display(c);
    
    if sq == 1
        % c em box individual: deve ser inserido
        if size(arqsq,1) < max
            % Arquivo tem espa�o: insere c.
            arqsq = [arqsq; c];
            display('arquivo tem espa�o: inserido');
        else
            % Arquivo est� cheio: substitui o pior.
            pr = pior(arqsq,1:size(arqsq,1),@compdist,L);
            arqsq(pr,:) = c;
            display(['arquivo cheio: substui ', num2str(pr)]);
        end
    else
        % c em box ocupado
        display('c em box ocupado');
        % Obt�m o outro indiv�duo do box.
        inbox = arqsq(hyinds,:);
        inbox(L.COLSQ) = sq;
        display(inbox);
        % Compara com c
        [m,~] = compdist(1,2,[c;inbox],L);
        if m == 1
            % c � melhor: substitui
            idx = indice(inbox,arqsq,L);
            display('substitui ');
            display(idx);
            arqsq(idx,:) = c;
        else
            % c � pior: descarta
            display('n�o substitui (descartado)');
            descartado = 1;
        end
    end
end

if ~descartado
    display('c inserido (recalcula)');
    % c foi inserido na popula��o: refaz os c�lculos de domina��o
    % e distribui��o.
    arqsq = pareto(arqsq,L);
    %[arqsq, sqinf, sqsup] = hyperbox(arqsq,resolucao,L);
    [arqsq, sqinf, sqsup, sqbox] = hyperbox2(arqsq,resolucao,L); % c�lculo do hyperbox
    %arqsq = naoaglomerado(arqsq,L);
    arqsq = selecionaDistribuicao(arqsq,L);
end
    
end

function [arqsq,sqinf,sqsup,sqbox,sqres] = atualizadis3(c,arqsq,sqinf,sqsup,sqbox,sqres, ...
                                           max,precisao,L)
%ATUALIZADIS Atualiza o arquivo de domina��o com um indiv�duo.
%   Tenta inserir um indiv�duo no arquivo de distribui��o.
%   Se o indiv�duo n�o for descartado, refaz os c�lculos 
%   das fronteiras de Pareto e do hyperbox.
%
%   Par�metros de entrada:
%     - c: indiv�duo a inserir;
%     - arqsq: popula��o (arquivo de distribui��o);
%     - sqinf: limites inferiores do grid para cada dimens�o;
%     - sqsup: limites superiores do grid para cada dimens�o;
%     - sqbox: boxes do arquivo;
%     - sqres: XXX
%     - max: n�mero m�ximo de elementos no arquivo;
%     - precisao: XXX
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - arqsq: popula��o atualizada.
%     - sqinf: limites inferiores atualizados;
%     - sqsup: limites superiores atualizados;
%     - sqbox: boxes atualizados;
%     - sqres: XXX

display('Introdu��o no arquivo de distribui��o');

% Verifica se indiv�duo n�o existe na popula��o
idx = indice(c,arqsq,L);
if idx > 0
    display('c j� existe: descartado');
    return;
end

descartado = 0;

if extrapola(c,sqinf,sqsup,L)
    % c estende grid do arquivo: deve ser inserido.
    display('c extrapola o grid');
    if size(arqsq,1) < max
        % Arquivo tem espa�o: insere c.
        arqsq = [arqsq; c];
        display('arquivo tem espa�o: inserido');
    else
        % Arquivo est� cheio: substitui o pior.
        pr = pior(arqsq,1:size(arqsq,1),@compdist,L);
        arqsq(pr,:) = c;
        display(['arquivo cheio: substui ', num2str(pr)]);
    end
else
    % C n�o estende grid: verifica de deve ser inserido
    display('c n�o extrapola grid');
    
    % Calcula o hyperbox e o squeeze factor de c.
    %hy = indbox(c,sqinf,sqsup,resolucao,L);
    %hy = indbox2(c,sqinf,sqsup,sqbox,resolucao,L);
    hy = indbox3(c,sqinf,sqsup,sqbox,sqres,L);
    hyinds = find(arqsq(:,L.COLHY) == hy);
    sq = size(hyinds,1) + 1;
    c(L.COLHY) = hy;
    c(L.COLSQ) = sq;
    
    % Calcula a fronteira de pareto do indiv�duo.
    c(L.COLPT) = indpar(c,arqsq,L);
    %display(c);
    
    if sq == 1
        % c em box individual: deve ser inserido
        if size(arqsq,1) < max
            % Arquivo tem espa�o: insere c.
            arqsq = [arqsq; c];
            display('arquivo tem espa�o: inserido');
        else
            % Arquivo est� cheio: substitui o pior.
            pr = pior(arqsq,1:size(arqsq,1),@compdist,L);
            arqsq(pr,:) = c;
            display(['arquivo cheio: substui ', num2str(pr)]);
        end
    else
        % c em box ocupado
        display('c em box ocupado');
        % Obt�m o outro indiv�duo do box.
        inbox = arqsq(hyinds,:);
        inbox(L.COLSQ) = sq;
        display(inbox);
        % Compara com c
        [m,~] = compdist(1,2,[c;inbox],L);
        if m == 1
            % c � melhor: substitui
            idx = indice(inbox,arqsq,L);
            display('substitui ');
            display(idx);
            arqsq(idx,:) = c;
        else
            % c � pior: descarta
            display('n�o substitui (descartado)');
            descartado = 1;
        end
    end
end

if ~descartado
    display('c inserido (recalcula)');
    % c foi inserido na popula��o: refaz os c�lculos de domina��o
    % e distribui��o.
    arqsq = pareto(arqsq,L);
    %[arqsq, sqinf, sqsup] = hyperbox(arqsq,resolucao,L);
    %[arqsq, sqinf, sqsup, sqbox] = hyperbox2(arqsq,resolucao,L); % c�lculo do hyperbox
    [arqsq, sqinf, sqsup, sqbox,sqres] = hyperbox3(arqsq,precisao,L); % c�lculo do hyperbox
    %arqsq = naoaglomerado(arqsq,L);
    arqsq = selecionaDistribuicao(arqsq,L);
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

function hb = indbox(c,inf,sup,resolucao,L)
%INDBOX Calcula o hyperbox de um indiv�duo.
%   Encontra o hyperbox ocupado por um indiv�duo em uma determinada
%   popula��o.
%
%   Par�metros de entrada:
%     - c: indiv�duo a verificar;
%     - inf: limites inferiores do grid (para todas as dimens�es);
%     - sup: limites superiores do grid (para todas as dimens�es);
%     - resolucao: tamanho do grid (o mesmo para todas as dimens�es);
%     - L: layout de um indiv�duo.
%
%   Par�metros de sa�da:
%     - hb: hyperbox do indiv�duo
  
  boxsize = (sup - inf) / resolucao;
  
  hbox = floor((c(L.COLF) - inf) ./ boxsize);
  
  hb = codifica(hbox,resolucao);
  
end

function hb = indbox2(c,inf,sup,boxes,resolucao,L)
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
  
  boxsize = (sup - inf) / resolucao;
  
  hbox = floor((c(L.COLF) - inf) ./ boxsize);
  
  ind = 0;
  for k = 1:length(boxes)
      if boxes(k).hbox == hbox
          ind = k;
          break
      end
  end
  
  %hb = codifica(hbox,resolucao);
  hb = ind;
  
end

function hb = indbox3(c,inf,sup,boxes,resolucao,L)
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
  
  boxsize = (sup - inf) ./ resolucao;
  
  hbox = floor((c(L.COLF) - inf) ./ boxsize);
  
  ind = 0;
  for k = 1:length(boxes)
      if boxes(k).hbox == hbox
          ind = k;
          break
      end
  end
  
  %hb = codifica(hbox,resolucao);
  hb = ind;
  
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

%{
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

%}
