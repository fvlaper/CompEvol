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

function [d,fo,gleo,geqo] = fitness(ft, f, gle, geq, r, s, n, pop)
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
%       - pop: array contendo a polula��o (um indiv�duo por linha).
%
%   Par�metros de sa�da:
%       - d: valor da fitness para cada indiv�duo (vetor coluna);
%       - fo: valor da fun��o de otimiza��o original
%         para cada indiv�duo (vetor coluna);
%       - gleo: valores da fun��o de restri��o de desigualdade
%         (npop linhas, nvar colunas);
%       - geqo: valores da fun��o de restri��o de igualdade
%         (npop linhas, nvar colunas). 

[d,fo,gleo,geqo] = ft(f, gle, geq, r, s, n, pop);
end

function [d,fo,gleo,geqo] = fitness_inv(f, gle, geq, r, s, n, pop)
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
%       - pop: array contendo a polula��o (um indiv�duo por linha).
%
%   Par�metros de sa�da:
%       - d: valor da fitness para cada indiv�duo (vetor coluna);
%       - fo: valor da fun��o de otimiza��o original
%         para cada indiv�duo (vetor coluna);
%       - gleo: valores da fun��o de restri��o de desigualdade
%         (npop linhas, nvar colunas);
%       - geqo: valores da fun��o de restri��o de igualdade
%         (npop linhas, nvar colunas). 

[h,fo,gleo,geqo] = fmod(f, gle, geq, r, s, pop);
d = 1 ./ (h - (min(h) - (10 ^ -n)));
end

function [d,fo,gleo,geqo] = fitness_desc(f, gle, geq, r, s, k, pop)
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
%       - pop: array contendo a polula��o (um indiv�duo por linha).
%
%   Par�metros de sa�da:
%       - d: valor da fitness para cada indiv�duo (vetor coluna);
%       - fo: valor da fun��o de otimiza��o original
%         para cada indiv�duo (vetor coluna);
%       - gleo: valores da fun��o de restri��o de desigualdade
%         (npop linhas, nvar colunas);
%       - geqo: valores da fun��o de restri��o de desigualdade
%         (npop linhas, nvar colunas). 

[h,fo,gleo,geqo] = fmod(f, gle, geq, r, s, pop);
cmax = k * sum(h) / size(h,1);
d = max(zeros(size(h)), (cmax-h));
end

function [h,fo,gleo,geqo] = fmod (f, gle, geq, r, s, pop)
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
%       - pop: array contendo a polula��o (um indiv�duo por linha).
%
%   Par�metros de sa�da:
%       - h: valor da fun��o de otimiza��o modificada
%         para cada indiv�duo (vetor coluna);
%       - fo: valor da fun��o de otimiza��o original
%         para cada indiv�duo (vetor coluna).
%       - gleo: valores da fun��o de restri��o de desigualdade
%         (npop linhas, nvar colunas);
%       - geqo: valores da fun��o de restri��o de igualdade
%         (npop linhas, nvar colunas). 

fo = f(pop);
gleo = gle(pop);
geqo = geq(pop);

h = fo ...                                             % func. original
  + r * sum((max(zeros(size(pop)), gleo) .^ 2) ,2) ... % desigualdades
  + s * sum(geqo .^ 2, 2);                             % igualdades
end

function [ x, f, g, h ] = ga(pop,ngen,nelite,ft,f,gle,geq,r,s,n)
%GA Fun��o principal do algoritmo gen�tico.
%   Executa o algoritmo gen�tico (ga) para otimiza��o da fun��o objetivo.
%
%   Par�metros de entrada:
%       - pop: array contendo a popula��o inicial;
%       - ngen: n�mero de gera��es a processar;
%       - nelite: n�mero de elementos da elite (preservados para a
%         pr�xima gera��o);
%       - ft: handle da fun��o de fitness;
%       - f: handle da fun��o de otimiza��o original;
%       - gle: handle da fun��o de restri��o de desigualdade;
%       - geq: handle da fun��o de restri��o de igualdade;
%       - r: par�metro de penalidade para a restri��o de desigualdade;
%       - s: par�metro de penalidade para a restri��o de igualdade;
%       - n: par�metro de escalonamento da fun��o de fitness.
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
%   - coluna (3*nvar+2): fitness.
columnsX = 1:size(pop,2);
columnsG = size(pop,2)+1 : 2*size(pop,2);
columnsH = 2*size(pop,2)+1 : 3*size(pop,2);
columnF = 3*size(pop,2)+1;
columnFit = 3*size(pop,2)+2;
pais = [pop zeros(size(pop)) zeros(size(pop)) zeros(size(pop,1),2)];

% C�lculos iniciais para a primeira popula��o de pais
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
%       - pop: array contendo a popula��o (unm indiv�duo por linha)
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
%       - pop: array contendo a popula��o (unm indiv�duo por linha)
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
%       - pop: array contendo a popula��o (unm indiv�duo por linha)
%
%   Par�metros de sa�da:
%       - g: array contendo a fun��o de restri��o (npop linhas,
%         nvar colunas). 

h = cos(2 * pi * pop) + 0.5;

end
