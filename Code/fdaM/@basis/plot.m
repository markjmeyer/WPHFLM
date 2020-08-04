function plot(basisobj, nx)
%  Plot a basis object.

%  last modified 6 January 2003

if nargin < 2, nx = 101;  end
 
rangex   = getbasisrange(basisobj);
x        = linspace(rangex(1),rangex(2),nx)';
basismat = full(eval_basis(basisobj, x));
plot (x, basismat, '-');
axis([x(1), x(nx), min(min(basismat)), max(max(basismat))])
