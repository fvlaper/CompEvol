
function f = dtlz2(x,m)

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
    s = s + (x(i)-0.5)^2;
end
g = s;

cosx = cos(x*pi/2);
sinx = sin(x*pi/2);

f = zeros(1,m);

f(1) =  (1+g) * prod(cosx(1:m-1));
for i = 2:m-1
    f(i) = (1+g) * prod(cosx(1:m-i)) * sinx(m-i+1);
end
f(m) = (1+g) * sinx(1);

%F = f(:);

end