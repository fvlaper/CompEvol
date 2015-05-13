function [p] = pot10(n)

i = find(n == 0);
if ~isempty(i)  % Evita logaritmo de 0
    n(i) = realmin;
end

f = log10(abs(n));
p = floor(f) + 1;
z = find((floor(f)-f) == 0); % É potência de 10 exata?
display(z);
if ~isempty(z)
   p(z) = f(z);
end

end

%function p = nextpow10(n)
%NEXTPOW10 Next higher power of 10.
%   NEXTPOW10(N) returns the first P such that 10^P >= abs(N).
%
%   Class support for input N or X:
%      float: double, single
%
%   See also LOG10.

%  Joe Henning - 25 Jul 2007

%f = log10(abs(n));

% Check if n is an exact power of 10.
%p = floor(f) + 1;
%i = find((floor(f)-f) == 0);
%if ~isempty(i)
%   p(i) = f(i);
%end

%end

%function p = nextpowof10(x)
%NEXTPOWOF10 Next power of 10.
%
%   P = NEXTPOWOF10(X) returns the smallest integer P such that 10^P >= abs(X).
%
%   Essentially, NEXTPOWOF10(X) is the same as CEIL(LOG(ABS(X)) / LOG(10)), but
%   special care is taken to catch round-off errors.
%
%   See also PREVPOWOF2, NEXTPOW.

%   Author:      Peter John Acklam
%   Time-stamp:  2003-11-17 11:38:03 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   %error(nargchk(1, 1, nargin));

   %if ~isreal(x)
   %   error('Input must be real.');
   %end

%   x = abs(x);
%   p = ceil(log(x) / log(10));          % estimate
%   k = x <= 10.^(p - 1);
%   p(k) = p(k) - 1;                     % correction
   
%end