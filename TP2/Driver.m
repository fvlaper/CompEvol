nvar = 10;
ncalc = 20000;
nexe = 50;


bench = 10 * nvar + nvar * (((-1/3)^2) - 10 * cos(2*pi*(-1/3)));
result = (-1/3) * ones(1,nvar);
s = 0;
r = 0;
f = 0;
for i = 1:nexe
    [a,b,c,d] = Flavio(ncalc,nvar);
    r = r + sum(abs(a-result));
    s = s + abs(bench-b);
    f = f + sum(a);
    display([b, sum(abs(round(a))), sum(c), sum(d)]);
end

display(bench);
display(s / nexe);
display(f / nexe);
display(r / nexe);

%display(bench);
%display(['bench = ', num2str(bench), ' calc = ', num2str(b)]);
%[a,b,c,d] = Flavio(ncalc,nvar)

%Teste(1.1+realmin)-1
%Teste([0 1 2; 1 0 2; 1 2 0])