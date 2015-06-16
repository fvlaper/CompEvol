
function f = dtlz1(x,m)

% x = indivíduo
% m = número de funções objetivo

x = x(:);

n = length(x);
k = n - m + 1;

%display(['n = ', num2str(n)]);
%display(['m = ', num2str(m)]);
%display(['k = ', num2str(k)]);

s = 0;
for i = m:n
    s = s + (x(i)-0.5)^2 - cos(20*pi*(x(i)-0.5));
end
g = 100*(k+s);

f = zeros(1,m);

f(1) = 0.5 * prod(x(1:m-1)) * (1+g);
for i = 2:m-1
    f(i) = 0.5 * prod(x(1:m-i)) * (1-x(m-i+1)) * (1+g);
end
f(m) = 0.5 * (1-x(1)) * (1+g);

%F = f(:);

end
