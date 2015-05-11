nvar = 2;
ncalc = 1500;
nexe = 30;
bench = 10 * nvar + nvar * (((-1/3)^2) - 10 * cos(2*pi*(-1/3)));
s = 0;
for i = 1:nexe
    [a,b,c,d] = Flavio(ncalc,nvar);
    s = s + abs(bench-b);
end

display(s / bench);
%display(bench);
%display(['bench = ', num2str(bench), ' calc = ', num2str(b)]);