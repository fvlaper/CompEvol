% Fl�vio Velloso Laper
%
% Computa��o evolucion�ria
% Prof. Jo�o Ant�nio de Vasconcelos
%
% Avalia��o 1

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
nexm = 200;
nit = floor((ncal / nexm) - 1);
if nit < 50
    nexm = 100;
    nit = floor((ncal / nexm) - 1);
    if nit < 50
        nexm = 50;
        nit = floor((ncal / nexm) - 1);
    end
end

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
%   - Fun��o objetivo penalizada: 1 coluna;
%   - Fun��o objetivo: 1 coluna;
%   - Vari�veis do melhor resultado: nvar colunas;
%   - Restri��es de desigualdade (melhor resultado): nvar colunas;
%   - Restri��es de igualdade (melhor resultado): nvar colunas;
%   - Viola��es de desigualdade (melhor resultado): 1 coluna;
%   - Viola��es de igualdade (melhor resultado): 1 coluna;
%   - Fun��o objetivo penalizada (melhor resultado): 1 coluna;
%   - Fun��o objetivo (melhor resultado): 1 coluna;
%   - Vari�veis do resultado lbest: nvar colunas;
%   - Restri��es de desigualdade (resultado lbest): nvar colunas;
%   - Restri��es de igualdade (resultado lbest): nvar colunas;
%   - Viola��es de desigualdade (resultado lbest): 1 coluna;
%   - Viola��es de igualdade (resultado lbest): 1 coluna;
%   - Fun��o objetivo penalizada (resultado lbest): 1 coluna;
%   - Fun��o objetivo (resultado lbest): 1 coluna;
%   - Velocidade (gradiente): nvar colunas;
%   - Label de vizinhan�a: 1 coluna.
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
colFL  = flo;
colV   = vi:vf;
colVZ  = vz;

% Determina��o da velocidade m�xima das part�culas
vmax = abs((xmax - xmin)) * 0.3;

% Tamanho da vizinhan�a social.
% Cada part�cula interage (em m�dia)
% com 10*p% do enxame.
p = 0.2;

% Estabelecimento do enxame inicial.
exm = exminit(nexm,nc,xmin,xmax,vmax,p,colX,colV,colVZ);

% Execu��o do algoritmo gen�tico para otimiza��o.
[x, f, p] = pso(exm,nit,fc,gle,geq,xmin,xmax,vmax,p, ...
                colX,colG,colH,colVG,colVH,colP,colF,colXB,colGB, ...
                colHB,colVGB,colVHB,colPB,colFB,colXL,colFL,colV,colVZ);

end

function [exm] = exminit (nexm, nc, xmin, xmax, vmax, p, colX, colV, colVZ)
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
%       - p: abrang�ncia da vizinha�a social;
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
        
% Inicializa a vizinhan�a.
exm(:,colVZ) = randi(round(nexm * p),nexm,1);

end

function [ x, f, p ] = pso(exm,nit,fc,gle,geq,xmin,xmax,vmax,pb, ...
           colX,colG,colH,colVG,colVH,colP,colF,colXB,colGB,colHB, ...
           colVGB,colVHB,colPB,colFB,colXL,colFL,colV,colVZ)
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
%       - vmax: valor m�ximo para a velocidade de uma part�cula;
%       - pb: abrang�ncia da vizinha�a social;
%       - colX: colunas das vari�veis;
%       - colG: colunas das restri��es de desigualdade;
%       - colH: colunas das restri��es de igualdade;
%       - colVG: coluna de viola��es de desigualdade;
%       - colVH: coluna de viola��es de igualdade;
%       - colP: coluna da fun��o objetivo penalizada;
%       - colF: coluna da fun��o objetivo;
%       - colXB: colunas das vari�veis (melhor resultado);
%       - colGB: colunas das restri��es de desigualdade (melhor resultado);
%       - colHB: colunas das restri��es de igualdade (melhor resultado);
%       - colVGB: coluna de viola��es de desigualdade (melhor resultado);
%       - colVHB: coluna de viola��es de igualdade (melhor resultado);
%       - colPB: coluna da fun��o objetivo penalizada (melhor resultado);
%       - colFB: coluna da fun��o objetivo (melhor resultado);
%       - colXL: colunas das vari�veis (resultado lbest);
%       - colGL: colunas das restri��es de desigualdade (resultado lbest);
%       - colHL: colunas das restri��es de igualdade (resultado lbest);
%       - colVGL: coluna de viola��es de desigualdade (resultado lbest);
%       - colVHL: coluna de viola��es de igualdade (resultado lbest);
%       - colPL: coluna da fun��o objetivo penalizada (resultado lbest);
%       - colFL: coluna da fun��o objetivo (resultado lbest);
%       - colV: colunas das velocidades;
%       - colVZ: coluna da vizinhan�a.
%
%   Par�metros de sa�da:
%       - x: vetor das vari�veis da melhor part�cula;
%       - f: melhor fun��o objetivo;
%       - p: valor da fun��o de penalidade avaliada em x.

% Faixas de colunas (para compara��es entre part�culas)
colD = min(colX):colF;    % Colunas com os dados da part�cula
colB = min(colXB):colFB;  % Colunas com os melhores dados da part�cula.
colL = min(colXL):colFL;  % Colunas com os dados lbest;

% Par�metros de penalidade para a fun��o objetivo modificada para
% problema sem restri��es. A penalidade pelo desrespeito �s restri��es
% aumenta a cada gera��o.
r = 1 : ((10^8 - 1) / nit) : 10^12;
s = 1 : ((10^8 - 1) / nit) : 10^12;

% Par�metros de peso de in�rcia
w = 0.9 : -((0.9 - 0.4) / nit) : 0.4;

% C�lculos iniciais para a primeira itera��o
exm = avaliacao(exm,fc,gle,geq,r(1),s(1),colX,colG,colH,colVG,colVH,colP,colF);

% O melhor resultado (local) � o �nico resultado
exm(:,colXB)  = exm(:,colX);
exm(:,colGB)  = exm(:,colG);
exm(:,colHB)  = exm(:,colH);
exm(:,colVGB) = exm(:,colVG);
exm(:,colVHB) = exm(:,colVH);
exm(:,colPB)  = exm(:,colP);
exm(:,colFB)  = exm(:,colF);

% Quantidade de part�culas e vari�veis.
npart = size(exm,1);
nvar = max(colX);

% Preenche os dados lbest de todas as part�culas
for p = 1:npart
    exm(p,:) = lbest(p,exm,colB,colL,colF,colVZ);
end

% Loop principal
t1max = 1.5;   % Valores m�ximos dos fatores
t2max = 2.51;   % de aprendizagem.
for it = 1:nit
%for it = 1:1
    
    % Determina os fatores de aprendizagem (randomizados).
    t1 = t1max * rand(npart,nvar); % Corresponde a c1 * rand.
    t2 = t2max * rand(npart,nvar); % Corresponde a c2 * rand.

    % Determina o coeficiente de contra��o
    % de modo a garantir a estabilidade do algoritmo.
    k = (2/(t1max+t2max-2)) * rand(npart,nvar);
    
    % Calcula as novas velocidades
    v = k .* (w(it) * exm(:,colV) + ...
             t1 .* (exm(:,colXB) - exm(:,colX)) + ...
             t2 .* (exm(:,colXL) - exm(:,colX)));
 
    % Mant�m velocidades dentro dos limites
    v = min(vmax, max(-vmax, v));
    
    % Atualiza as posi��es e as velocidades
    exm(:,colX) = exm(:,colX) + v;
    exm(:,colV) = v;
    
    % Mant�m as vari�veis dentro dos limites
    exm(:,colX) = min(xmax, max(xmin, exm(:,colX)));
    
    % Avalia os novos dados
    exm = avaliacao(exm,fc,gle,geq,r(it+1),s(it+1), ...
                    colX,colG,colH,colVG,colVH,colP,colF);

    % Calcula os novos melhores valores
    exm = calcBest(exm,colD,colB,colVG,colVH,colF);

    % Gera novas vizinhan�as
    exm(:,colVZ) = randi(round(npart * pb),npart,1);
    
    % Preenche os dados lbest de todas as part�culas
    for p = 1:npart
        exm(p,:) = lbest(p,exm,colB,colL,colF,colVZ);
    end

end

% Retorna o melhor resultado
exm = sortrows(exm,colP);

best = 1;
x = exm(best,colX);
f = exm(best,colF);
p = sum(max(0,exm(best,colG))) + sum(abs(exm(best,colH)));

end

function part =  lbest(ind,exm,colB,colL,colF,colVZ)
%LBEST Calcula o lbest para uma part�cula.
%   Faz o c�lculo dos valores lbest para uma part�cula. Estes valores
%   s�o diferentes para cada part�cula, dependendo da topologia
%   de vizinhan�a adotada.
%
%   Par�metros de entrada:
%       - ind: �ndice da part�cula;
%       - exm: enxame;
%       - colB: colunas com os melhores valores das part�culas;
%       - colL: colunas com os dados lbest;
%       - colVG: coluna de viola��es de desigualdade;
%       - colVH: coluna de viola��es de igualdade;
%       - colF: coluna da fun��o objetivo;
%       - colVZ: coluna do label de vizinhan�a.
%
%   Par�metros de sa�da:
%       - part: part�cula com valores lbest atualizados.

% Determina a topologia a utilizar.
r = rand;

part = exm(ind,:);

% Determina a topologia da vizinhan�a.
if r < (1/3)
    lexm = vring(ind,exm);
elseif r < (2/3)
    lexm = vmesh(ind,exm);
else
    lexm = vsocial(part,exm,colVZ);
end

% Procura os melhores valores da vizinhan�a
melhor = lexm(1,colB);

for p = 2:size(lexm,1)
    melhor = comparaValor(melhor,lexm(p,colB),colF);
end

part(colL) = melhor;

end

function lexm = vsocial(part,exm,colVZ)
%VSOCIAL Retorna a vizinhan�a social.
%   Retorna a vizinhan�a social da part�cula, baseado no label de
%   vizinhan�a.
%
%   Par�metros de entrada:
%       - part: part�cula.
%       - exm: enxame.
%       - colVZ: coluna do label de vizinhan�a.
%
%   Par�metros de sa�da:
%       - lexm: vizinhan�a da part�cula.

vz = part(colVZ);
lexm = exm(exm(:,colVZ) == vz,:);

end

function lexm = vring(ind,exm)
%VRING Retorna a vizinhan�a (topologia de anel).
%   Retorna a vizinhan�a da part�cula, baseado em uma topologia de anel.
%
%   Par�metros de entrada:
%       - ind: �ndice da part�cula;
%       - exm: enxame.
%
%   Par�metros de sa�da:
%       - lexm: vizinhan�a da part�cula.

k = 5;  % anel tem k part�culas antes e k depois.
nexm = size(exm,1);
aux = [exm(nexm-k+1:nexm,:); exm; exm(1:k,:)];
lexm = aux(ind:ind+2*k,:);

end

function lexm = vmesh(ind,exm)
%VRING Retorna a vizinhan�a (topologia von Neumann).
%   Retorna a vizinhan�a da part�cula, baseado em uma topologia
%   von Neumann.
%
%   Par�metros de entrada:
%       - ind: �ndice da part�cula;
%       - exm: enxame.
%
%   Par�metros de sa�da:
%       - lexm: vizinhan�a da part�cula.

nexm = size(exm,1);
k = floor(sqrt(nexm));  % enxame distribu�do como um quadrado.
norte = ind-k; if norte < 1,    norte = norte + nexm; end;
sul   = ind+k; if sul   > nexm, sul   = sul   - nexm; end
leste = ind+1; if leste > nexm, leste = leste - nexm; end
oeste = ind-1; if oeste < 1,    oeste = oeste + nexm; end
lexm = [exm(norte,:); exm(oeste,:); exm(ind,:); exm(leste,:); exm(sul,:)];

end

function [exm] = calcBest(exm,colD,colB,colVG,colVH,colF)
%CALCBEST Calcula o melhor valor local.
%   Calcula o melhor valor j� encontrado por uma part�cula
%   durante as itera��es.
%
%   Par�metros de entrada:
%       - exm: enxame;
%       - colD: colunas com os valores das part�culas;
%       - colB: colunas com os melhores valores das part�culas;
%       - colVG: coluna de viola��es de desigualdade;
%       - colVH: coluna de viola��es de igualdade;
%       - colF: coluna da fun��o objetivo.
%
%   Par�metros de sa�da:
%       - exm: enxame atualizado.

for p = 1:size(exm,1)
    melhor = comparaFaixa(exm(p,colD),exm(p,colB),colVG,colVH,colF);
    exm(p,colB) = melhor;
end

end

function [melhor] = comparaFaixa(part1,part2,colVG,colVH,colO)
%COMPARA Compara duas part�culas e retorna a melhor.
%   Compara os valores de fun��o objetivo e viola��es de restri��es
%   e retorna os da mais adequada.
%
%   Par�metros de entrada:
%       - part1: valores da part�cula 1;
%       - part2: valores da part�cula 2;
%       - colVG: coluna de viola��es de desigualdade;
%       - colVH: coluna de viola��es de igualdade;
%       - colO: coluna da fun��o objetivo.
%
%   Par�metros de sa�da:
%       - melhor: valores da melhor part�cula.

% Quantidade de viola��es das part�culas.
v1g = round(part1(colVG));
v1h = round(part1(colVH));
v2g = round(part2(colVG));
v2h = round(part2(colVH));


% Seleciona a que viola menos restri��es; em caso de empate,
% a de menor fun��o objetivo. Prioriza as restri��es de
% igualdade.
if v1h < v2h
    melhor = part1;
elseif v2h < v1h
    melhor = part2;
elseif v1g < v2g
    melhor = part1;
elseif v1g > v2g
    melhor = part2;
elseif part1(colO) <= part2(colO)
    melhor = part1;
else
    melhor = part2;
end

end

function [melhor] = comparaValor(part1,part2,colO)
%COMPARA Compara duas part�culas e retorna a melhor.
%   Compara os valores de fun��o objetivo
%   e retorna os da mais adequada.
%
%   Par�metros de entrada:
%       - part1: valores da part�cula 1;
%       - part2: valores da part�cula 2;
%       - colO: coluna da fun��o objetivo.
%
%   Par�metros de sa�da:
%       - melhor: valores da melhor part�cula.

% Compara a fun��o objetivo
if part1(colO) <= part2(colO)
    melhor = part1;
else
    melhor = part2;
end

end

function [exm] = avaliacao(exm,fc,gle,geq,r,s,colX,colG,colH,colVG,colVH,colP,colF)
%AVALIACAO Avalia os par�metros de uma part�cula.
%   Calcula os valores das diversas colunas das part�culas de um enxame
%   a partir do valor das vari�veis.
%
%   Par�metros de entrada:
%       - exm: enxame;
%       - fc: handle da fun��o de otimiza��o original;
%       - gle: handle da fun��o de restri��o de desigualdade;
%       - geq: handle da fun��o de restri��o de igualdade;
%       - r: par�metro de penalidade para a restri��o de desigualdade;
%       - s: par�metro de penalidade para a restri��o de igualdade;
%       - colX: colunas das vari�veis;
%       - colG: colunas das restri��es de desigualdade;
%       - colH: colunas das restri��es de igualdade;
%       - colVG: coluna de viola��es de desigualdade;
%       - colVH: coluna de viola��es de igualdade;
%       - colP: coluna da fun��o objetivo penalizada.
%       - colF: coluna da fun��o objetivo.
%
%   Par�metros de sa�da:
%       - exm: enxame atualizado.

exm(:,colF) = fc(exm(:,colX));  % fun��o objetivo
exm(:,colG) = gle(exm(:,colX)); % restri��es de desigualdade
exm(:,colH) = geq(exm(:,colX)); % restri��es de igualdade

% Fun��o objetivo penalizada
exm(:,colP) = exm(:,colF) + r * sum((max(0, exm(:,colG)) .^ 2),2) + ...
                            s * sum((exm(:,colH) .^ 2),2);

exm(:,colVG) = sum(max(0,exm(:,colG)),2); % viola��es de desigualdade

exm(:,colVH) = sum(abs(exm(:,colH)),2);   % viola��es de igualdade

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
