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

pop = popinit(3,3,-5.12,5.12);

x = [1,1,1];
f = 2;
g = pop;
%h = rastrigin(pop);
%h = restr_desig(pop);
h = restr_ig(pop);

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

function [g] = restr_desig (pop)
%RESTR_DESIG Fun��o de restri��o de desigualdade.
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

function [h] = restr_ig (pop)
%RESTR_IG Fun��o de restri��o de igualdade.
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