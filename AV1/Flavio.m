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
%   - Função objetivo penalizada: 1 coluna;
%   - Função objetivo: 1 coluna;
%   - Variáveis do melhor resultado: nvar colunas;
%   - Restrições de desigualdade (melhor resultado): nvar colunas;
%   - Restrições de igualdade (melhor resultado): nvar colunas;
%   - Violações de desigualdade (melhor resultado): 1 coluna;
%   - Violações de igualdade (melhor resultado): 1 coluna;
%   - Função objetivo penalizada (melhor resultado): 1 coluna;
%   - Função objetivo (melhor resultado): 1 coluna;
%   - Variáveis do resultado lbest: nvar colunas;
%   - Restrições de desigualdade (resultado lbest): nvar colunas;
%   - Restrições de igualdade (resultado lbest): nvar colunas;
%   - Violações de desigualdade (resultado lbest): 1 coluna;
%   - Violações de igualdade (resultado lbest): 1 coluna;
%   - Função objetivo penalizada (resultado lbest): 1 coluna;
%   - Função objetivo (resultado lbest): 1 coluna;
%   - Velocidade (gradiente): nvar colunas;
%   - Label de vizinhança: 1 coluna.
xi  = 1;     xf  = xi + nvar - 1;
gi  = xf+1;  gf  = gi + nvar -1;
hi  = gf+1;  hf  = hi + nvar -1;
vg  = hf+1;
vh  = vg+1;
fp  = vh+1;
fo  = fp+1;
xbi = fo+1;  xbf = xbi + nvar -1;
gbi = xbf+1; gbf = gbi + nvar -1;
hbi = gbf+1; hbf = hbi + nvar -1;
vgb = hbf+1;
vhb = vgb+1;
fpb = vhb+1;
fbo = fpb+1;
xli = fbo+1;  xlf = xli + nvar -1;
gli = xlf+1; glf = gli + nvar -1;
hli = glf+1; hlf = hli + nvar -1;
vgl = hlf+1;
vhl = vgl+1;
flp = vhl+1;
flo = flp+1;
vi  = flo+1; vf  = vi + nvar - 1;
vz  = vf+1;
nc  = vz;    % Total de colunas;

colX   = xi:xf;
colG   = gi:gf;
colH   = hi:hf;
colVG  = vg;
colVH  = vh;
colP   = fp;
colF   = fo;
colXB  = xbi:xbf;
colGB  = gbi:gbf;
colHB  = hbi:hbf;
colVGB = vgb;
colVHB = vhb;
colPB  = fpb;
colFB  = fbo;
colXL  = xli:xlf;
colGL  = gli:glf;
colHL  = hli:hlf;
colVGL = vgl;
colVHL = vhl;
colPL  = flp;
colFL  = flo;
colV   = vi:vf;
colVZ  = vz;

% Determinação da velocidade máxima das partículas
vmax = abs((xmax - xmin)) * 0.3;

% Estabelecimento do enxame inicial.
exm = exminit(nexm,nc,xmin,xmax,vmax,colX,colV,colVZ);

% Execução do algoritmo genético para otimização.
[x, f, p] = pso(exm,nit,fc,gle,geq,xmin,xmax,vmax, ...
                colX,colG,colH,colVG,colVH,colP,colF,colXB,colGB,colHB,colVGB,colVHB,colPB,colFB,colXL,colGL,colHL,colVGL,colVHL,colPL,colFL,colV,colVZ);

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
%p = 0.1;
p = 0.5;
exm(:,colVZ) = randi(round(nexm * p),nexm,1);

end

function [ x, f, p ] = pso(exm,nit,fc,gle,geq,xmin,xmax,vmax, ...
           colX,colG,colH,colVG,colVH,colP,colF,colXB,colGB,colHB,colVGB,colVHB,colPB,colFB,colXL,colGL,colHL,colVGL,colVHL,colPL,colFL,colV,colVZ)
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
%       - vmax: valor máximo para a velocidade de uma partícula;
%       - colX: colunas das variáveis;
%       - colG: colunas das restrições de desigualdade;
%       - colH: colunas das restrições de igualdade;
%       - colVG: coluna de violações de desigualdade;
%       - colVH: coluna de violações de igualdade;
%       - colP: coluna da função objetivo penalizada;
%       - colF: coluna da função objetivo;
%       - colXB: colunas das variáveis (melhor resultado);
%       - colGB: colunas das restrições de desigualdade (melhor resultado);
%       - colHB: colunas das restrições de igualdade (melhor resultado);
%       - colVGB: coluna de violações de desigualdade (melhor resultado);
%       - colVHB: coluna de violações de igualdade (melhor resultado);
%       - colPB: coluna da função objetivo penalizada (melhor resultado);
%       - colFB: coluna da função objetivo (melhor resultado);
%       - colXL: colunas das variáveis (resultado lbest);
%       - colGL: colunas das restrições de desigualdade (resultado lbest);
%       - colHL: colunas das restrições de igualdade (resultado lbest);
%       - colVGL: coluna de violações de desigualdade (resultado lbest);
%       - colVHL: coluna de violações de igualdade (resultado lbest);
%       - colPL: coluna da função objetivo penalizada (resultado lbest);
%       - colFL: coluna da função objetivo (resultado lbest);
%       - colV: colunas das velocidades;
%       - colVZ: coluna da vizinhança.
%
%   Parâmetros de saída:
%       - x: vetor das variáveis da melhor partícula;
%       - f: melhor função objetivo;
%       - p: valor da função de penalidade avaliada em x.

% Faixas de colunas (para comparações entre partículas)
colD = min(colX):colF;    % Colunas com os dados da partícula
colB = min(colXB):colFB;  % Colunas com os melhores dados da partícula.
colL = min(colXL):colFL;  % Colunas com os dados lbest;

% Parâmetros de penalidade para a função objetivo modificada para
% problema sem restrições. A penalidade pelo desrespeito às restrições
% aumenta a cada geração.
r = 1 : ((10^8 - 1) / nit) : 10^8;
s = 1 : ((10^8 - 1) / nit) : 10^8;

% Parâmetros de peso de inércia
w = 0.9 : -((0.9 - 0.4) / nit) : 0.4;
k = 0.73; % ATENÇÃO: variar


% Cálculos iniciais para a primeira iteração
exm = avaliacao(exm,fc,gle,geq,r(1),s(1),colX,colG,colH,colVG,colVH,colP,colF);

% O melhor resultado (local) é o único resultado
exm(:,colXB)  = exm(:,colX);
exm(:,colGB)  = exm(:,colG);
exm(:,colHB)  = exm(:,colH);
exm(:,colVGB) = exm(:,colVG);
exm(:,colVHB) = exm(:,colVH);
exm(:,colPB)  = exm(:,colP);
exm(:,colFB)  = exm(:,colF);

% Quantidade de partículas
npart = size(exm,1);

% Preenche os dados lbest de todas as partículas
for p = 1:npart
    exm(p,:) = lbest(p,exm,colB,colL,colVG,colVH,colP,colVZ);
end

% Loop principal
c1 = 2.05; % parâmetro de aprendizagem cognitiva. (1.49)
c2 = 2.05; % parâmetro de aprendizagem social.    (1.49)
for it = 1:nit
%for it = 1:1
    
    % Calcula as novas velocidades
    v = k * (w(it) * exm(:,colV) + ...
     c1 * rand([npart, (max(colX))]) .* (exm(:,colX) - exm(:,colXB)) + ...
     c2 * rand([npart, (max(colX))]) .* (exm(:,colX) - exm(:,colXL)));
 
    % Mantém velocidades dentro dos limites
    v = min(vmax, max(-vmax, v));
    
    % Atualiza as posições e as velocidades
    display(exm(1:5,:));
    exm(:,colX) = exm(:,colX) + v;
    exm(:,colV) = v;
    
    % Avalia os novos dados
    exm = avaliacao(exm,fc,gle,geq,r(it+1),s(it+1),colX,colG,colH,colVG,colVH,colP,colF);
    exm = calcBest(exm,colD,colB,colVG,colVH,colF);

    % Preenche os dados lbest de todas as partículas
    for p = 1:npart
        exm(p,:) = lbest(p,exm,colB,colL,colVG,colVH,colP,colVZ);
    end

end

display(exm);

% Retorna o melhor resultado
best = 1;
x = exm(best,colX);
f = exm(best,colF);
p = sum(max(0,exm(best,colG))) + sum(abs(exm(best,colH)));

end

function part =  lbest(ind,exm,colB,colL,colVG,colVH,colF,colVZ)
%LBEST Calcula o lbest para uma partícula.
%   Faz o cálculo dos valores lbest para uma partícula. Estes valores
%   são diferentes para cada partícula, dependendo da topologia
%   de vizinhança adotada.
%
%   Parâmetros de entrada:
%       - ind: índice da partícula;
%       - exm: enxame;
%       - colB: colunas com os melhores valores das partículas;
%       - colL: colunas com os dados lbest;
%       - colVG: coluna de violações de desigualdade;
%       - colVH: coluna de violações de igualdade;
%       - colF: coluna da função objetivo;
%       - colVZ: coluna do label de vizinhança.
%
%   Parâmetros de saída:
%       - part: partícula com valores lbest atualizados.

% Determina a topologia a utilizar.
r = rand;

part = exm(ind,:);

% Determina a topologia da vizinhança.
if r < (1/3)
    lexm = vring(ind,exm);
elseif r < (2/3)
    lexm = vmesh(ind,exm);
else
    lexm = vsocial(part,exm,colVZ);
end

% Procura os melhores valores da vizinhança
melhor = lexm(1,colB);

for p = 2:size(lexm,1)
    melhor = compara(melhor,lexm(p,colB),colVG,colVH,colF);
end

part(colL) = melhor;

end

function lexm = vsocial(part,exm,colVZ)
%VSOCIAL Retorna a vizinhança social.
%   Retorna a vizinhança social da partícula, baseado no label de
%   vizinhança.
%
%   Parâmetros de entrada:
%       - part: partícula.
%       - exm: enxame.
%       - colVZ: coluna do label de vizinhança.
%
%   Parâmetros de saída:
%       - lexm: vizinhança da partícula.

vz = part(colVZ);
lexm = exm(exm(:,colVZ) == vz,:);

end

function lexm = vring(ind,exm)
%VRING Retorna a vizinhança (topologia de anel).
%   Retorna a vizinhança da partícula, baseado em uma topologia de anel.
%
%   Parâmetros de entrada:
%       - ind: índice da partícula;
%       - exm: enxame.
%
%   Parâmetros de saída:
%       - lexm: vizinhança da partícula.

k = 2;  % anel tem k partículas antes e k depois.
nexm = size(exm,1);
aux = [exm(nexm-k+1:nexm,:); exm; exm(1:k,:)];
lexm = aux(ind:ind+2*k,:);

end

function lexm = vmesh(ind,exm)
%VRING Retorna a vizinhança (topologia von Neumann).
%   Retorna a vizinhança da partícula, baseado em uma topologia
%   von Neumann.
%
%   Parâmetros de entrada:
%       - ind: índice da partícula;
%       - exm: enxame.
%
%   Parâmetros de saída:
%       - lexm: vizinhança da partícula.

nexm = size(exm,1);
k = floor(sqrt(nexm));  % enxame distribuído como um quadrado.
norte = ind-k; if norte < 1,    norte = norte + nexm; end;
sul   = ind+k; if sul   > nexm, sul   = sul   - nexm; end
leste = ind+1; if leste > nexm, leste = leste - nexm; end
oeste = ind-1; if oeste < 1,    oeste = oeste + nexm; end
lexm = [exm(norte,:); exm(oeste,:); exm(ind,:); exm(leste,:); exm(sul,:)];

end

function [exm] = calcBest(exm,colD,colB,colVG,colVH,colF)
%CALCBEST Calcula o melhor valor local.
%   Calcula o melhor valor já encontrado por uma partícula
%   durante as iterações.
%
%   Parâmetros de entrada:
%       - exm: enxame;
%       - colD: colunas com os valores das partículas;
%       - colB: colunas com os melhores valores das partículas;
%       - colVG: coluna de violações de desigualdade;
%       - colVH: coluna de violações de igualdade;
%       - colF: coluna da função objetivo.
%
%   Parâmetros de saída:
%       - exm: enxame atualizado.

for p = 1:size(exm,1)
    melhor = compara(exm(p,colD),exm(p,colB),colVG,colVH,colF);
    exm(p,colB) = melhor;
end

end

function [melhor] = compara(part1,part2,colVG,colVH,colF)
%COMPARA Compara duas partículas e retorna a melhor.
%   Compara os valores de função objetivo e violações de restrições
%   e retorna os da mais adequada.
%
%   Parâmetros de entrada:
%       - part1: valores da partícula 1;
%       - part2: valores da partícula 2;
%       - colVG: coluna de violações de desigualdade;
%       - colVH: coluna de violações de igualdade;
%       - colF: coluna da função objetivo.
%
%   Parâmetros de saída:
%       - melhor: valores da melhor partícula.

% Quantidade de violações das partículas.
%v1 = part1(colVG) + part1(colVH);
%v2 = part2(colVG) + part2(colVH);

% Seleciona a que viola menos restrições; em caso de empate,
% a de menor função objetivo.
%if v1 < v2
%    melhor = part1;
%elseif v2 < v1
%    melhor = part2;
%elseif part1(colF) <= part2(colF)
%    melhor = part1;
%else
%    melhor = part2;
%end

% Compara apenas a função objetivo
if part1(colF) <= part2(colF)
    melhor = part1;
else
    melhor = part2;
end

end

function [exm] = avaliacao(exm,fc,gle,geq,r,s,colX,colG,colH,colVG,colVH,colP,colF)
%AVALIACAO Avalia os parâmetros de uma partícula.
%   Calcula os valores das diversas colunas das partículas de um enxame
%   a partir do valor das variáveis.
%
%   Parâmetros de entrada:
%       - exm: enxame;
%       - fc: handle da função de otimização original;
%       - gle: handle da função de restrição de desigualdade;
%       - geq: handle da função de restrição de igualdade;
%       - r: parâmetro de penalidade para a restrição de desigualdade;
%       - s: parâmetro de penalidade para a restrição de igualdade;
%       - colX: colunas das variáveis;
%       - colG: colunas das restrições de desigualdade;
%       - colH: colunas das restrições de igualdade;
%       - colVG: coluna de violações de desigualdade;
%       - colVH: coluna de violações de igualdade;
%       - colP: coluna da função objetivo penalizada.
%       - colF: coluna da função objetivo.
%
%   Parâmetros de saída:
%       - exm: enxame atualizado.

exm(:,colF) = fc(exm(:,colX));  % função objetivo
exm(:,colG) = gle(exm(:,colX)); % restrições de desigualdade
exm(:,colH) = geq(exm(:,colX)); % restrições de igualdade

% Função objetivo penalizada
%exm(:,colP) = exm(:,colF) + r * sum((max(0, exm(:,colG)) .^ 2),2) + ...
%                            s * sum((exm(:,colH) .^ 2),2);

exm(:,colVG) = sum((exm(:,colG) > 0),2);  % violações de desigualdade
%exm(:,colVG) = sum(max(0,exm(:,colG)),2); % violações de desigualdade

exm(:,colVH) = sum(pot10(exm(:,colH)),2); % violações de igualdade
%exm(:,colVH) = sum((exm(:,colH) ~= 0),2); % violações de igualdade
%exm(:,colVH) = sum(abs(exm(:,colH)),2);   % violações de igualdade

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
