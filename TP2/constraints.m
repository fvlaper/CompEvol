% Função de restrição para o cálculo utilizando o toolbox do matlab.
% Chamar com optimtool('ga')
% Informar:
%   - Solver: ga - Genetic Algorithm;
%   - Fitness function: @rastriginsfcn;
%   - Number of variables: 2 (ou 5, 10, ...);
%   - Bounds: Lower: -5.12; Upper: 5.12;
%   - Nonlinear constraint function: @constraints;
%   - Selection function: Tournament;
%   - Tournament size: 2;
%   - Elite count: 1;
%   - Crossover fraction: 0.6.

function [c, ceq] = constraints(x)
c = sin(2*pi*x)+0.5;
ceq = cos(2*pi*x)+0.5;
end