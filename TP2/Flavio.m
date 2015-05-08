function [ x, f, g, h ] = Flavio( ncal, nvar )
%FLAVIO Trabalho 2 de Computa��o Evolucion�ria - Fl�vio Laper.
%   Esta fun��o implementa um algoritmo gen�tico para minimizar a
%   fun��o de Rastrigin para nvar vari�veis com restri��es.
%   
%   Par�metros de entrada:
%       - ncal: n�mero m�ximo de chamadas da fun��o de c�lculo da
%         fitness do problema;
%       - nvar: n�mero de vari�veis.
%
%   Par�metros de sa�da:
%       - x: vetor das vari�veis de decis�o do melhor indiv�duo;
%       - f: melhor fun��o objetivo;
%       - g: vetor restri��o de desigualdade avaliado no ponto x;
%       - h: vetor restri��o de igualdade avaliado no ponto x.

% Estabelecimento dos parametros iniciais
npop = 10;                % n�mero de indiv�duos na popula��o.
ngen = floor((ncal / npop) - 1); % n�mero de gera��es. Calculado a partir
                                 % do n�mero m�ximo de c�lculos da fun��o
                                 % de fitness, considerando npop c�lculos 
                                 % por gera��o, mais npop c�lculos para 
                                 % a popula��o inicial.
                          
f   = @rastrigin;         % Handle da fun��o objetivo.
gle = @rastrigin_le;      % Handle da fun��o de restri��o de desigualdade.
geq = @rastrigin_eq;      % Handle da fun��o de restri��o de igualdade.
xmin = -5.12;             % Valores de contorno para as vari�veis
xmax = 5.12;              % de decis�o.
                          
% Estabelecimento da popula��o inicial.
pop = popinit(npop,nvar,xmin,xmax);

% Execu��o do algoritmo gen�tico para otimiza��o.
[x, f, g, h] = ga(pop, ngen, f, gle, geq);

end

function [pop] = popinit ( npop, nvar, xmin, xmax )
%POPINIT Gera��o da popula��o inicial.
%   Gera a popula��o com npop indiv�duos e nvar vari�veis. Cada
%   vari�vel tem valores limitados � faixa [xmin,xmax].
%
%   Par�metros de entrada:
%       - npop: n�mero de indiv�duos da popula��o;
%       - nvar: n�mero de vari�veis por indiv�duo;
%       - xmin: valor m�nimo de uma vari�vel;
%       - xmax: valor m�ximo de uma vari�vel.
%
%   Par�metros de sa�da:
%       - pop: array contendo a popula��o (npop linhas, nvar colunas). 

pop = xmin * ones(npop,nvar) + (xmax - xmin) * rand(npop,nvar);

end

function [ x, f, g, h ] = ga(pop,ngen,f,gle,geq)
%GA Fun��o principal do algoritmo gen�tico.
%   Executa o algoritmo gen�tico (ga) para otimiza��o da fun��o objetivo.
%
%   Par�metros de entrada:
%       - pop: array contendo a popula��o inicial;
%       - ngen: n�mero de gera��es a processar;
%       - f: handle da fun��o de otimiza��o original;
%       - gle: handle da fun��o de restri��o de desigualdade;
%       - geq: handle da fun��o de restri��o de igualdade.
%
%   Par�metros de sa�da:
%       - x: vetor das vari�veis de decis�o do melhor indiv�duo;
%       - f: melhor fun��o objetivo;
%       - g: vetor restri��o de desigualdade avaliado no ponto x;
%       - h: vetor restri��o de igualdade avaliado no ponto x.

% Popula��o de pais para a pr�xima gera��o. Array contendo uma linha
% por indiv�duo da popula��o. Para cada linha:
%   - colunas 1 - nvar: vari�veis de decis�o;
%   - colunas (nvar+1) - (2*nvar): restri��es de desigualdade;
%   - colunas (2*nvar+1) - (3*nvar): restri��es de igualdade;
%   - coluna (3*nvar+1): fun��o objetivo;
%   - coluna (3*nvar+2): fitness;
%   - coluna (3*nvar+3): n�mero de restri��es violadas.
nvar = size(pop,2);
columnsX  = 1:nvar;
columnsG  = nvar+1 : 2*nvar;
columnsH  = 2*nvar+1 : 3*nvar;
columnF   = 3*nvar+1;
%columnFit = 3*nvar+2;
%columnV   = 3*nvar+3;
pais = [pop zeros(size(pop)) zeros(size(pop)) zeros(size(pop,1),3)];

%Par�metros do algoritmo
ft = @fitness_desc; % handle da fun��o de fitness.
fp = @fitpar_desc;  % handle da fun��o de par�metros da fitness.

nelite = 1;         % n�mero de elementos da elite (preservados para a
                    % pr�xima gera��o).

% Par�metros de penalidade para a fun��o objetivo modificada para
% problema sem restri��es. A penalidade pelo desrespeito �s restri��es
% aumenta a cada gera��o.
r = 1 : ((10^8 - 1) / ngen) : 10^8;
s = 1 : ((10^8 - 1) / ngen) : 10^8;

% Par�metros da fun��o de fitness
fitpar = fp(ngen);

% C�lculos iniciais para a primeira popula��o de pais
[pais] = fitness(ft, f, gle, geq, r(1), s(1), fitpar(1), pais, nvar);
display(pais);
pais = popSort(pais,nvar);
display(pais);

% Inicializa��o da popula��o de filhos
filhos = zeros(size(pais,1)-nelite, size(pais,2));

% Loop de evolu��o
%for g = 1 : ngen
for g = 1 : 1
    % Guarda os indiv�duos da elite
    elite = pais(1 : nelite, :);
    
    % Selec�o de pais
    selecionados = selecao(pais,size(filhos,1));
    display(selecionados);
    
    % Realiza os cruzamentos
    filhos = cruzamento(filhos,pais,selecionados(1),selecionados(2),nvar,1,0.5);
    display(filhos);
end

x = pais(1,columnsX);
f = pais(1,columnF);
g = pais(1,columnsG);
h = pais(1,columnsH);

end

function [filhos] = cruzamento(filhos,pais,pai1,pai2,nvar,ind,pc)
%CRUZAMENTO Faz o cruzamento entre dois indiv�duos.
%   Realiza o cruzamento estre os dois indiv�duos cujos �ndices s�o
%   fornecidos e coloca os filhos sequencialmente no vetor de filhos.
%   Se o cruzamento n�o for realizado, apenas copia os pais para o
%   vetor de filhos.
%
%   Par�metros de entrada:
%       - filhos: array dos filhos;
%       - pais: array contendo os pais;
%       - pai1: �ndice do primenro pai;
%       - pai2: �ndice do segundo pai;
%       - nvar: n�mero de vari�veis;
%       - ind: �ndice a primeira linha vaga no array de filhos;
%       - pc: probabilidade de realiza��o do cruzamento.
%
%   Par�metros de sa�da:
%       - filhos: array de filhos (preenchido com novos indiv�duos).


% Verifica se o cruzamento deve ser realizado; caso contr�rio, copia.
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
%SELECTION Seleciona os pais que participar�o dos pr�ximos processos.
%   Utiliza torneio bin�rio.
%
%   Par�metros de entrada:
%       - pais: array contendo os pais;
%       - nfilhos: n�mero de filhos.
%
%   Par�metros de sa�da:
%       - filhos: vetor com �ndices dos pais selecionados.

npais = size(pais,1);   % quantidade de pais.
ultimo = -1;            % �ndice do �ltimo pai selecionado.
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
    
    % Evita o cruzamento de un indiv�duo consigo mesmo
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
%POPSORT Ordena o array de popula�ao.
%   Ordena o array de popula��o do melhor para o pior indiv�duo.
%
%   Par�metros de entrada:
%       - pop: array de popula��o;
%       - nvar: n�mero de vari�veis.
%
%   Par�metros de sa�da:
%       - pop: array de popula��o ordenado.

columnFit = 3*nvar+2; % Coluna do valor de fitness.
columnV = 3*nvar+3;   % Coluna da quantidade de viola��es de restri��o.

%Ordena pela fun��o de fitness apenas (decrescente).
%pop = sortrows(pop, -columnFit);

%Ordena pela quantidade de viola��es e pela fitness (decrescente).
pop = sortrows(pop,[columnV -columnFit]);
end

function [v] = violacoes(pop, nvar)
%VIOLACOES Calcula o n�mero restri��es violadas por um indiv�duo.
%   Retorna, para cada indiv�duo, o n�mero de restri��es (igualdade
%   e desigualdade) violadas por suas vari�veis).
%
%   Par�metros de entrada:
%       - pop: array de popula��o.
%       - nvar: n�mero de vari�veis.
%
%   Par�metros de sa�da:
%       - v: n�mero de restri��es violadas por indiv�duo(vetor coluna).

columnsG  = nvar+1 : 2*nvar;
columnsH  = 2*nvar+1 : 3*nvar;

% Varia��o admiss�vel no teste de viola��o da restri��o de igualdade
% (testes de igualdade com valores em ponto flutuante n�o s�o precisos).
delta = 0.0000001;

%v = sum((pop(:,columnsG) > 0),2) + sum((pop(:,columnsH) ~= 0),2);
v = sum((pop(:,columnsG) > 0),2) + sum(abs(pop(:,columnsH)) > delta,2);
end

function [pop] = fitness(ft, f, gle, geq, r, s, n, pop, nvar)
%FITNESS Fun��o de fitness.
%   Calcula a fitness dos indiv�duos da popula��o. Esta fun��o �
%   apenas um "stub" que recebe a fun��o de fitness espec�fica a
%   utilizar no c�lculo.
%
%   Par�metros de entrada:
%       - ft: handle da fun��o de fitness;
%       - f: handle da fun��o de otimiza��o original;
%       - gle: handle da fun��o de restri��o de desigualdade;
%       - geq: handle da fun��o de restri��o de igualdade;
%       - r: par�metro de penalidade para a restri��o de desigualdade;
%       - s: par�metro de penalidade para a restri��o de igualdade;
%       - n: par�metro de escalonamento da fun��o de fitness;
%       - pop: array contendo a popula��o (um indiv�duo por linha);
%       - nvar: n�mero de vari�veis.
%
%   Par�metros de sa�da:
%       - pop: array de popula��o com valores calculados.

columnFit = 3*nvar+2;
columnV = 3*nvar+3;

pop = ft(f, gle, geq, r, s, n, pop, nvar);
pop(:,columnFit) = escalonamento(pop,nvar);
pop(:,columnV) = violacoes(pop,nvar);
end

function [e] = escalonamento(pop,nvar)
%ESCALONAMENTO Faz o escalonamento da fitness dos indiv�duos.
%   Escalona a fitness dos indiv�duos para evitar que um "superindiv�duo"
%   domine a popula��o. Utiliza escalonamento linear.
%
%   Par�metros de entrada:
%       - pop: array de popula��o.
%       - nvar: n�mero de vari�veis.
%
%   Par�metros de sa�da:
%       - e: fitness escalonada (vetor coluna).

columnFit = 3*nvar+2; % Coluna da fitness
Cmax = 2; % n�mero m�ximo de c�pias do melhor indiv�duo.

fit = pop(:,columnFit);
fmed = sum(fit) / size(pop,1);
fmax = max(fit);
a = (Cmax-1)*fmed / (fmax-fmed);
b = (1-a) / fmed;
e = max((a * fit + b), zeros(size(fit)));
end

function [pop] = fitness_inv(f, gle, geq, r, s, n, pop, nvar)
%FITNESS_INV Fun��o de fitness pelo m�todo de invers�o.
%   Calcula a fitness dos indiv�duos da popula��o pelo m�todo de
%   invers�o (transformando um problema de minimiza��o em maximiza��o).
%
%   Par�metros de entrada:
%       - f: handle da fun��o de otimiza��o original;
%       - gle: handle da fun��o de restri��o de desigualdade;
%       - geq: handle da fun��o de restri��o de igualdade;
%       - r: par�metro de penalidade para a restri��o de desigualdade;
%       - s: par�metro de penalidade para a restri��o de igualdade;
%       - n: par�metro de escalonamento para a invers�o;
%       - pop: array contendo a polula��o (um indiv�duo por linha);
%       - nvar: n�mero de vari�veis.
%
%   Par�metros de sa�da:
%       - pop: array de popula��o com valores calculados.

columnFit = 3*nvar+2;

[h,pop] = fmod(f, gle, geq, r, s, pop, nvar);
d = 1 ./ (h - (min(h) - (10 ^ -n)));
pop(:,columnFit) = d;
end

function [par] = fitpar_inv(ngen)
%FITPAR_DESC Par�metro para c�lculo da fitness pelo m�todo de invers�o.
%   Retorna, para cada gera��o, o valor do par�metro n para o c�lculo
%   da fitness pelo m�todo de deslocamento.
%
%   Par�metros de entrada:
%       - ngen: n�mero de gera��es.
%
%   Par�metros de sa�da:
%       - vetor com (ngen+1) valores do par�metro (um para cada gera��o,
%         mais a popula��o inicial).

par = 0.1 : ((1.0 - 0.1) / ngen) : 1.0;
end

function [pop] = fitness_desc(f, gle, geq, r, s, k, pop, nvar)
%FITNESS_DESC Fun��o de fitness pelo m�todo de deslocamento.
%   Calcula a fitness dos indiv�duos da popula��o pelo m�todo de
%   deslocamento (transformando um problema de minimiza��o em maximiza��o).
%
%   Par�metros de entrada:
%       - f: handle da fun��o de otimiza��o original;
%       - gle: handle da fun��o de restri��o de desigualdade;
%       - geq: handle da fun��o de restri��o de igualdade;
%       - r: par�metro de penalidade para a restri��o de desigualdade;
%       - s: par�metro de penalidade para a restri��o de igualdade;
%       - k: par�metro de escalonamento para o deslocamento;
%       - pop: array contendo a polula��o (um indiv�duo por linha);
%       - nvar: n�mero de vari�veis.
%
%   Par�metros de sa�da:
%       - pop: array de popula��o com valores calculados.

columnFit = 3*nvar+2;

[h,pop] = fmod(f, gle, geq, r, s, pop, nvar);
cmax = k * sum(h) / size(h,1);
d = max(zeros(size(h)), (cmax-h));
pop(:,columnFit) = d;
end

function [par] = fitpar_desc(ngen)
%FITPAR_DESC Par�metro para c�lculo da fitness pelo m�todo de deslocamento.
%   Retorna, para cada gera��o, o valor do par�metro K para o c�lculo
%   da fitness pelo m�todo de deslocamento. Ser� utilizado um valor
%   fixo para todas as gera��es.
%
%   Par�metros de entrada:
%       - ngen: n�mero de gera��es.
%
%   Par�metros de sa�da:
%       - vetor com (ngen+1) valores do par�metro (um para cada gera��o,
%         mais a popula��o inicial).

k = 1.2; % Valor fixo
par = k * ones(1,ngen+1);
end

function [h,pop] = fmod (f, gle, geq, r, s, pop, nvar)
%FMOD Fun��o de otimiza��o modificada.
%   Calcula o valor da fun��o de otimiza��o (por indiv�duo) modificada
%   para um problema irrestrito.
%
%   Par�metros de entrada:
%       - f: handle da fun��o de otimiza��o original;
%       - gle: handle da fun��o de restri��o de desigualdade;
%       - geq: handle da fun��o de restri��o de igualdade;
%       - r: par�metro de penalidade para a restri��o de desigualdade;
%       - s: par�metro de penalidade para a restri��o de igualdade;
%       - pop: array contendo a polula��o (um indiv�duo por linha);
%       - nvar: n�mero de vari�veis.
%
%   Par�metros de sa�da:
%       - h: valor da fun��o de otimiza��o modificada
%         para cada indiv�duo (vetor coluna);
%       - pop: array de popula��o com valores calculados.

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
%ROLETA Sele��o pelo m�todo da roleta.
%   Seleciona os pais para a pr�xima gera��o pelo m�todo da roleta.
%
%   Par�metros de entrada:
%       - pop: array contendo a popula��o (sup�e-se que a �ltima
%         coluna cont�m a fitness do indiv�duo).
%
%   Par�metros de sa�da:
%       - pais: popula��o de pais selecionados para os procedimentos
%         gen�ticos seguintes.

%npop = size(pop,1);                    % tamanho da popula��o
%fitcol = size(pop,2);                  % coluna com o valor da fitness
%totfitness = sum(pop(:,fitcol));       % fitness total
%pais = zeros(size(pop));               % pais selecionados

% Loop de sele��o
%for n = 1:npop
%    r = totfitness * rand(); % valor do lance    
%    i = 1;                   % �ndice do indiv�duo selecionado
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
% Fun��o de Rastrigin e restri��es
%------------------------------------------------------------------------%

function [f] = rastrigin (pop)
%RASTRIGIN Fun��o de Rastrigin.
%   Calcula a fun��o de Rastrigin para todos os indiv�duos da popula��o.
%
%   Par�metros de entrada:
%       - pop: array contendo a popula��o (unm indiv�duo por linha).
%
%   Par�metros de sa�da:
%       - f: valor da fun��o de Rastrigin para cada indiv�duo (vetor
%         coluna com npop linhas).

f = 10 * size(pop,2) + sum(pop .^ 2 - 10 * cos(2 * pi * pop), 2);

end

function [g] = rastrigin_le (pop)
%RASTRIGIN_LE Fun��o de restri��o de desigualdade - Rastrigin.
%   Calcula o valor da restri��o de desigualdade para cada vari�vel
%   de todos os indiv�duos da popula��o.
%
%   Par�metros de entrada:
%       - pop: array contendo a popula��o (um indiv�duo por linha).
%
%   Par�metros de sa�da:
%       - g: array contendo a fun��o de restri��o (npop linhas,
%         nvar colunas). 

g = sin(2 * pi * pop) + 0.5;

end

function [h] = rastrigin_eq (pop)
%RASTRIGIN_EQ Fun��o de restri��o de igualdade - Rastrigin.
%   Calcula o valor da restri��o de igualdade para cada vari�vel
%   de todos os indiv�duos da popula��o.
%
%   Par�metros de entrada:
%       - pop: array contendo a popula��o (um indiv�duo por linha).
%
%   Par�metros de sa�da:
%       - g: array contendo a fun��o de restri��o (npop linhas,
%         nvar colunas). 

h = cos(2 * pi * pop) + 0.5;

end
