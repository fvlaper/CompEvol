function [ x, f, p ] = Flavio( ncal, nvar )
%FLAVIO Avaliação 1 de Computação Evolucionária - Flávio Laper.
%   Esta função implementa um algoritmo PSO para minimizar a
%   função de Rastrigin para nvar variáveis com restrições.
%   
%   Parâmetros de entrada:
%       - ncal: número máximo de chamadas da função objetivo.
%       - nvar: número de variáveis.
%
%   Parâmetros de saída:
%       - x: vetor das variáveis da melhor partícula;
%       - f: melhor função objetivo;
%       - p: valor da função de penalidade avaliada em x.

% Número de partículas do exame
% Calculado a partir
% do número máximo de cálculos da função
% objetivo, considerando npop cálculos 
% por iteracao, mais npop cálculos para 
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

fc   = @rastrigin;        % Handle da função objetivo.
gle = @rastrigin_le;      % Handle da função de restrição de desigualdade.
geq = @rastrigin_eq;      % Handle da função de restrição de igualdade.
xmin = -5.12;             % Valores de contorno para as variáveis
xmax = 5.12;              % de decisão.

% Layout do array de enxame: array contendo uma linha
% por partícula. Para cada linha:
%   - Variáveis de decisão: nvar colunas;
%   - Restrições de desigualdade: nvar colunas;
%   - Restrições de igualdade: nvar colunas;
%   - Violações de desigualdade: 1 coluna;
%   - Violações de igualdade: 1 coluna;
%   - Função objetivo: 1 coluna;
%   - Variáveis do melhor resultado: nvar colunas;
%   - Restrições de desigualdade (melhor resultado): nvar colunas;
%   - Restrições de igualdade (melhor resultado): nvar colunas;
%   - Violações de desigualdade (melhor resultado): 1 coluna;
%   - Violações de igualdade (melhor resultado): 1 coluna;
%   - Função objetivo (melhor resultado): 1 coluna;
%   - Velocidade (gradiente): nvar colunas;
%   - Label de vizinhança: 1 coluna.
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

% Determinação da velocidade máxima das partículas
vmax = abs((xmax - xmin)) * 0.2;

% Estabelecimento do enxame inicial.
exm = exminit(nexm,nc,xmin,xmax,vmax,colX,colV,colVZ);

% Execução do algoritmo genético para otimização.
[x, f, p] = pso(exm,nit,fc,gle,geq,xmin,xmax, ...
                colX,colG,colH,colVG,colVH,colF,colXB,colGB,colHB,colVGB,colVHB,colFB,colV,colVZ);

end

function [exm] = exminit (nexm, nc, xmin, xmax, vmax, colX, colV, colVZ)
%EXMINIT Geração do enxame inicial.
%   Gera o enxame com nexm partículas Cada
%   variável tem valores limitados à faixa [xmin,xmax].
%
%   Parâmetros de entrada:
%       - nexm: número de partículas do enxame;
%       - nc: número de colunas;
%       - xmin: valor mínimo de uma variável;
%       - xmax: valor máximo de uma variável;
%       - vmax: velocidade máxima das partículas;
%       - colX: colunas das variáveis;
%       - colV: coluna das velocidades;
%       - colVZ: coluna da vizinhança.
%
%   Parâmetros de saída:
%       - exm: array contendo o enxame.

% Aloca o array.
exm = zeros(nexm,nc);

% Inicializa as variáveis.
exm(:,colX) = xmin * ones(nexm,max(colX))  + ...
            (xmax - xmin) * rand(nexm,max(colX));
        
% Inicializa as velocidades.
exm(:,colV) = 2 * vmax * rand(nexm,max(colX)) - vmax;
        
% Inicializa a vizinhança. Cada partícula interage (em média)
% com 10*p% do enxame;
p = 0.1;
exm(:,colVZ) = randi(round(nexm * p),nexm,1);

end

function [ x, f, p ] = pso(exm,nit,fc,gle,geq,xmin,xmax, ...
           colX,colG,colH,colVG,colVH,colF,colXB,colGB,colHB,colVGB,colVHB,colFB,colV,colVZ)
%PSO Função principal do algoritmo.
%   Executa o algoritmo particle swarm (pso) para otimização
%   da função objetivo.
%
%   Parâmetros de entrada:
%       - exm: enxame inicial;
%       - nit: número de iterações;
%       - fc: handle da função de otimização original;
%       - gle: handle da função de restrição de desigualdade;
%       - geq: handle da função de restrição de igualdade;
%       - xmin: valor mínimo de uma variável;
%       - xmax: valor máximo de uma variável;
%       - colX: colunas das variáveis;
%       - colG: colunas das restrições de desigualdade;
%       - colH: colunas das restrições de igualdade;
%       - colVG: coluna de violações de desigualdade;
%       - colVH: coluna de violações de igualdade;
%       - colF: coluna da função objetivo;
%       - colXB: colunas das variáveis (melhor resultado);
%       - colGB: colunas das restrições de desigualdade (melhor resultado);
%       - colHB: colunas das restrições de igualdade (melhor resultado);
%       - colVGB: coluna de violações de desigualdade (melhor resultado);
%       - colVHB: coluna de violações de igualdade (melhor resultado);
%       - colFB: coluna da função objetivo (melhor resultado);
%       - colV: colunas das velocidades;
%       - colVZ: coluna da vizinhança.
%
%   Parâmetros de saída:
%       - x: vetor das variáveis da melhor partícula;
%       - f: melhor função objetivo;
%       - p: valor da função de penalidade avaliada em x.

% Cálculos iniciais para a primeira iteração
exm = avaliacao(exm,fc,gle,geq,colX,colG,colH,colVG,colVH,colF);

% O melhor resultado (local) é o único resultado
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
%AVALIACAO Avalia os parâmetros de uma partícula.
%   Calcula os valores das diversas colunas das partículas de um enxame
%   a partir do valor das variáveis.
%
%   Parâmetros de entrada:
%       - exm: enxame;
%       - fc: handle da função de otimização original;
%       - gle: handle da função de restrição de desigualdade;
%       - geq: handle da função de restrição de igualdade;
%       - colX: colunas das variáveis;
%       - colG: colunas das restrições de desigualdade;
%       - colH: colunas das restrições de igualdade;
%       - colVG: coluna de violações de desigualdade;
%       - colVH: coluna de violações de igualdade;
%       - colF: coluna da função objetivo.
%
%   Parâmetros de saída:
%       - exm: enxame atualizado.

exm(:,colF) = fc(exm(:,colX));  % função objetivo
exm(:,colG) = gle(exm(:,colX)); % restrições de desigualdade
exm(:,colH) = geq(exm(:,colX)); % restrições de igualdade

%exm(:,colVG) = sum((exm(:,colG) > 0),2);  % violações de desigualdade
exm(:,colVG) = sum(max(0,exm(:,colG)),2); % violações de desigualdade

%exm(:,colVH) = sum(pot10(exm(:,colH)),2); % violações de igualdade
%exm(:,colVH) = sum((exm(:,colH) ~= 0),2); % violações de igualdade
exm(:,colVH) = sum(abs(exm(:,colH)),2);   % violações de igualdade

end

function p = pot10(n)
%POT10 Próxima potência de 10.
%   Calcula p tal que 10^p >= abs(n).
%
%   Parâmetros de entrada:
%       - n: array de números a examinar.
%
%   Parâmetros de saída:
%       - p: array com as próximas potências de 10.

i = find(n == 0);
if ~isempty(i)  % Evita logaritmo de 0
    n(i) = realmin;
end

f = log10(abs(n));
p = floor(f) + 1;
z = find((floor(f)-f) == 0); % É potência de 10 exata?
if ~isempty(z)
   p(z) = f(z);
end

end
%------------------------------------------------------------------------%
% Função de Rastrigin e restrições
%------------------------------------------------------------------------%

function [f] = rastrigin (exm)
%RASTRIGIN Função de Rastrigin.
%   Calcula a função de Rastrigin para todas as partículas do enxame.
%
%   Parâmetros de entrada:
%       - exm: array contendo o enxame (uma partícula por linha).
%
%   Parâmetros de saída:
%       - f: valor da função de Rastrigin para cada partícula.

f = 10 * size(exm,2) + sum(exm .^ 2 - 10 * cos(2 * pi * exm), 2);

end

function [g] = rastrigin_le (pop)
%RASTRIGIN_LE Função de restrição de desigualdade - Rastrigin.
%   Calcula o valor da restrição de desigualdade para cada variável
%   de todas as partículas do enxame.
%
%   Parâmetros de entrada:
%       - exm: array contendo o enxame (uma partícula por linha).
%
%   Parâmetros de saída:
%       - g: array contendo a função de restrição.

g = sin(2 * pi * pop) + 0.5;

end

function [h] = rastrigin_eq (exm)
%RASTRIGIN_EQ Função de restrição de igualdade - Rastrigin.
%   Calcula o valor da restrição de igualdade para cada variável
%   de todas as partículas do enxame.
%
%   Parâmetros de entrada:
%       - exm: array contendo o enxame (uma partícula por linha).
%
%   Parâmetros de saída:
%       - g: array contendo a função de restrição. 

h = cos(2 * pi * exm) + 0.5;

end
