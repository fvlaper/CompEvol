function [ x, f, p ] = Flavio( ncal, nvar )
%FLAVIO Avalia��o 1 de Computa��o Evolucion�ria - Fl�vio Laper.
%   Esta fun��o implementa um algoritmo PSO para minimizar a
%   fun��o de Rastrigin para nvar vari�veis com restri��es.
%   
%   Par�metros de entrada:
%       - ncal: n�mero m�ximo de chamadas da fun��o objetivo.
%       - nvar: n�mero de vari�veis.
%
%   Par�metros de sa�da:
%       - x: vetor das vari�veis da melhor part�cula;
%       - f: melhor fun��o objetivo;
%       - p: valor da fun��o de penalidade avaliada em x.

% N�mero de part�culas do exame
% Calculado a partir
% do n�mero m�ximo de c�lculos da fun��o
% objetivo, considerando npop c�lculos 
% por iteracao, mais npop c�lculos para 
% o enxame inicial.
nexm = 10;
%nexm = 200;
nit = floor((ncal / nexm) - 1);
%if nit < 50
%    nexm = 100;
%    nit = floor((ncal / nexm) - 1);
%    if nit < 50
%        nexm = 50;
%        nit = floor((ncal / nexm) - 1);
%    end
%end

fc   = @rastrigin;        % Handle da fun��o objetivo.
gle = @rastrigin_le;      % Handle da fun��o de restri��o de desigualdade.
geq = @rastrigin_eq;      % Handle da fun��o de restri��o de igualdade.
xmin = -5.12;             % Valores de contorno para as vari�veis
xmax = 5.12;              % de decis�o.

% Layout do array de enxame: array contendo uma linha
% por part�cula. Para cada linha:
%   - Vari�veis de decis�o: nvar colunas;
%   - Restri��es de desigualdade: nvar colunas;
%   - Restri��es de igualdade: nvar colunas;
%   - Viola��es de desigualdade: 1 coluna;
%   - Viola��es de igualdade: 1 coluna;
%   - Fun��o objetivo: 1 coluna;
%   - Vari�veis do melhor resultado: nvar colunas;
%   - Restri��es de desigualdade (melhor resultado): nvar colunas;
%   - Restri��es de igualdade (melhor resultado): nvar colunas;
%   - Viola��es de desigualdade (melhor resultado): 1 coluna;
%   - Viola��es de igualdade (melhor resultado): 1 coluna;
%   - Fun��o objetivo (melhor resultado): 1 coluna;
%   - Velocidade (gradiente): nvar colunas;
%   - Label de vizinhan�a: 1 coluna.
xi  = 1;     xf  = xi + nvar - 1;
gi  = xf+1;  gf  = gi + nvar -1;
hi  = gf+1;  hf  = hi + nvar -1;
vg  = hf+1;
vh  = vg+1;
fo  = vh+1;
xbi = fo+1;  xbf = xbi + nvar -1;
gbi = xbf+1; gbf = gbi + nvar -1;
hbi = gbf+1; hbf = hbi + nvar -1;
vgb = hbf+1;
vhb = vgb+1;
fbo = vhb+1;
vi  = fbo+1; vf  = vi + nvar - 1;
vz  = vf+1;
nc  = vz;    % Total de colunas;

colX   = xi:xf;
colG   = gi:gf;
colH   = hi:hf;
colVG  = vg;
colVH  = vh;
colF   = fo;
colXB  = xbi:xbf;
colGB  = gbi:gbf;
colHB  = hbi:hbf;
colVGB = vgb;
colVHB = vhb;
colFB  = fbo;
colV   = vi:vf;
colVZ  = vz;

% Determina��o da velocidade m�xima das part�culas
vmax = abs((xmax - xmin)) * 0.2;

% Estabelecimento do enxame inicial.
exm = exminit(nexm,nc,xmin,xmax,vmax,colX,colV,colVZ);

% Execu��o do algoritmo gen�tico para otimiza��o.
[x, f, p] = pso(exm,nit,fc,gle,geq,xmin,xmax, ...
                colX,colG,colH,colVG,colVH,colF,colXB,colGB,colHB,colVGB,colVHB,colFB,colV,colVZ);

end

function [exm] = exminit (nexm, nc, xmin, xmax, vmax, colX, colV, colVZ)
%EXMINIT Gera��o do enxame inicial.
%   Gera o enxame com nexm part�culas Cada
%   vari�vel tem valores limitados � faixa [xmin,xmax].
%
%   Par�metros de entrada:
%       - nexm: n�mero de part�culas do enxame;
%       - nc: n�mero de colunas;
%       - xmin: valor m�nimo de uma vari�vel;
%       - xmax: valor m�ximo de uma vari�vel;
%       - vmax: velocidade m�xima das part�culas;
%       - colX: colunas das vari�veis;
%       - colV: coluna das velocidades;
%       - colVZ: coluna da vizinhan�a.
%
%   Par�metros de sa�da:
%       - exm: array contendo o enxame.

% Aloca o array.
exm = zeros(nexm,nc);

% Inicializa as vari�veis.
exm(:,colX) = xmin * ones(nexm,max(colX))  + ...
            (xmax - xmin) * rand(nexm,max(colX));
        
% Inicializa as velocidades.
exm(:,colV) = 2 * vmax * rand(nexm,max(colX)) - vmax;
        
% Inicializa a vizinhan�a. Cada part�cula interage (em m�dia)
% com 10*p% do enxame;
p = 0.1;
exm(:,colVZ) = randi(round(nexm * p),nexm,1);

end

function [ x, f, p ] = pso(exm,nit,fc,gle,geq,xmin,xmax, ...
           colX,colG,colH,colVG,colVH,colF,colXB,colGB,colHB,colVGB,colVHB,colFB,colV,colVZ)
%PSO Fun��o principal do algoritmo.
%   Executa o algoritmo particle swarm (pso) para otimiza��o
%   da fun��o objetivo.
%
%   Par�metros de entrada:
%       - exm: enxame inicial;
%       - nit: n�mero de itera��es;
%       - fc: handle da fun��o de otimiza��o original;
%       - gle: handle da fun��o de restri��o de desigualdade;
%       - geq: handle da fun��o de restri��o de igualdade;
%       - xmin: valor m�nimo de uma vari�vel;
%       - xmax: valor m�ximo de uma vari�vel;
%       - colX: colunas das vari�veis;
%       - colG: colunas das restri��es de desigualdade;
%       - colH: colunas das restri��es de igualdade;
%       - colVG: coluna de viola��es de desigualdade;
%       - colVH: coluna de viola��es de igualdade;
%       - colF: coluna da fun��o objetivo;
%       - colXB: colunas das vari�veis (melhor resultado);
%       - colGB: colunas das restri��es de desigualdade (melhor resultado);
%       - colHB: colunas das restri��es de igualdade (melhor resultado);
%       - colVGB: coluna de viola��es de desigualdade (melhor resultado);
%       - colVHB: coluna de viola��es de igualdade (melhor resultado);
%       - colFB: coluna da fun��o objetivo (melhor resultado);
%       - colV: colunas das velocidades;
%       - colVZ: coluna da vizinhan�a.
%
%   Par�metros de sa�da:
%       - x: vetor das vari�veis da melhor part�cula;
%       - f: melhor fun��o objetivo;
%       - p: valor da fun��o de penalidade avaliada em x.

% C�lculos iniciais para a primeira itera��o
exm = avaliacao(exm,fc,gle,geq,colX,colG,colH,colVG,colVH,colF);

% O melhor resultado (local) � o �nico resultado
exm(:,colXB)  = exm(:,colX);
exm(:,colGB)  = exm(:,colG);
exm(:,colHB)  = exm(:,colH);
exm(:,colVGB) = exm(:,colVG);
exm(:,colVHB) = exm(:,colVH);
exm(:,colFB)  = exm(:,colF);
display(exm);

% Retorna o melhor resultado
best = 1;
x = exm(best,colX);
f = exm(best,colF);
p = sum(max(0,exm(best,colG))) + sum(abs(exm(best,colH)));

end

function [exm] = avaliacao(exm,fc,gle,geq,colX,colG,colH,colVG,colVH,colF)
%AVALIACAO Avalia os par�metros de uma part�cula.
%   Calcula os valores das diversas colunas das part�culas de um enxame
%   a partir do valor das vari�veis.
%
%   Par�metros de entrada:
%       - exm: enxame;
%       - fc: handle da fun��o de otimiza��o original;
%       - gle: handle da fun��o de restri��o de desigualdade;
%       - geq: handle da fun��o de restri��o de igualdade;
%       - colX: colunas das vari�veis;
%       - colG: colunas das restri��es de desigualdade;
%       - colH: colunas das restri��es de igualdade;
%       - colVG: coluna de viola��es de desigualdade;
%       - colVH: coluna de viola��es de igualdade;
%       - colF: coluna da fun��o objetivo.
%
%   Par�metros de sa�da:
%       - exm: enxame atualizado.

exm(:,colF) = fc(exm(:,colX));  % fun��o objetivo
exm(:,colG) = gle(exm(:,colX)); % restri��es de desigualdade
exm(:,colH) = geq(exm(:,colX)); % restri��es de igualdade

%exm(:,colVG) = sum((exm(:,colG) > 0),2);  % viola��es de desigualdade
exm(:,colVG) = sum(max(0,exm(:,colG)),2); % viola��es de desigualdade

%exm(:,colVH) = sum(pot10(exm(:,colH)),2); % viola��es de igualdade
%exm(:,colVH) = sum((exm(:,colH) ~= 0),2); % viola��es de igualdade
exm(:,colVH) = sum(abs(exm(:,colH)),2);   % viola��es de igualdade

end

function p = pot10(n)
%POT10 Pr�xima pot�ncia de 10.
%   Calcula p tal que 10^p >= abs(n).
%
%   Par�metros de entrada:
%       - n: array de n�meros a examinar.
%
%   Par�metros de sa�da:
%       - p: array com as pr�ximas pot�ncias de 10.

i = find(n == 0);
if ~isempty(i)  % Evita logaritmo de 0
    n(i) = realmin;
end

f = log10(abs(n));
p = floor(f) + 1;
z = find((floor(f)-f) == 0); % � pot�ncia de 10 exata?
if ~isempty(z)
   p(z) = f(z);
end

end
%------------------------------------------------------------------------%
% Fun��o de Rastrigin e restri��es
%------------------------------------------------------------------------%

function [f] = rastrigin (exm)
%RASTRIGIN Fun��o de Rastrigin.
%   Calcula a fun��o de Rastrigin para todas as part�culas do enxame.
%
%   Par�metros de entrada:
%       - exm: array contendo o enxame (uma part�cula por linha).
%
%   Par�metros de sa�da:
%       - f: valor da fun��o de Rastrigin para cada part�cula.

f = 10 * size(exm,2) + sum(exm .^ 2 - 10 * cos(2 * pi * exm), 2);

end

function [g] = rastrigin_le (pop)
%RASTRIGIN_LE Fun��o de restri��o de desigualdade - Rastrigin.
%   Calcula o valor da restri��o de desigualdade para cada vari�vel
%   de todas as part�culas do enxame.
%
%   Par�metros de entrada:
%       - exm: array contendo o enxame (uma part�cula por linha).
%
%   Par�metros de sa�da:
%       - g: array contendo a fun��o de restri��o.

g = sin(2 * pi * pop) + 0.5;

end

function [h] = rastrigin_eq (exm)
%RASTRIGIN_EQ Fun��o de restri��o de igualdade - Rastrigin.
%   Calcula o valor da restri��o de igualdade para cada vari�vel
%   de todas as part�culas do enxame.
%
%   Par�metros de entrada:
%       - exm: array contendo o enxame (uma part�cula por linha).
%
%   Par�metros de sa�da:
%       - g: array contendo a fun��o de restri��o. 

h = cos(2 * pi * exm) + 0.5;

end
