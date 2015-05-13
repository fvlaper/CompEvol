nvar = 2;
ncalc = 15000;
nexe = 30;
bench = 10 * nvar + nvar * (((-1/3)^2) - 10 * cos(2*pi*(-1/3)));
s = 0;
for i = 1:nexe
%for i = 1:1
    [a,b,c,d] = Flavio(ncalc,nvar);
    display(a);
    s = s + abs(bench-b);
end

display(s / bench);
%display(bench);
%display(['bench = ', num2str(bench), ' calc = ', num2str(b)]);
%[a,b,c,d] = Flavio(ncalc,nvar)

%Teste(1.1+realmin)-1
%Teste([0 1 2; 1 0 2; 1 2 0])