nvar = 2;
ncalc = 2000;
nexe = 1;


bench = 10 * nvar + nvar * (((-1/3)^2) - 10 * cos(2*pi*(-1/3)));
result = (-1/3) * ones(1,nvar);
deltax = 0;
ft = 0;
pt = 0;
for i = 1:nexe
    [x,f,p] = Flavio(ncalc,nvar);
    %display([x, f, p]);
    deltax = deltax + sum(abs(x-result));
    ft = ft + abs(bench-f);
    pt = pt + p;
end

%display(bench);
%display(['deltax = ', num2str(deltax / nexe), ...
%         ' ft = ', num2str(ft / nexe), ...
%         ' pt = ', num2str(pt / nexe)]);
