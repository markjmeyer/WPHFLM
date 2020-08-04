function evalarray = eval_bifd(bifd, sevalarg, tevalarg, sLfd, tLfd)
%EVAL_BIFD  evaluates a bi-functional data object BIFD
%  at argument values in arrays SEVALARG and TEVALARG.
%  SLfd and TLfd are either integers giving the order of derivative,
%  or linear differential operators to be applied before evaluation.
%  Their defaults are 0, meaning that the function itself is evaluated.

%  last modified 30 January 2003

%  check that at least three arguments are present

if nargin < 3
    error('There are less than three arguments.');
end

%  check SEVALARG

sizeevalarg = size(sevalarg);
if sizesevalarg(1) > 1 & sizesevalarg(2) > 1
    error('Argument EVALARG is not a vector.');
end
sevalarg = sevalarg(:);

%  check TEVALARG

sizeevalarg = size(tevalarg);
if sizetevalarg(1) > 1 & sizetevalarg(2) > 1
    error('Argument EVALARG is not a vector.');
end
tevalarg = tevalarg(:);

if nargin < 5
    tLfd = int2Lfd(0);
end
if nargin < 4
    sLfd = int2Lfd(0);
end

if ~isa_bifd(bifd)
    error('Argument BIFD is not a bivariate functional data object.');
end
ns = length(sevalarg);
sbasisobj = getsbasis(bifd);
snbasis   = getnbasis(sbasisobj);

if ~isa_Lfd(sLfd)
    error ('SLfd is not a linear differential operator object.');
end

if ~isa_Lfd(tLfd)
    error ('TLfd is not a linear differential operator object.');
end

snderiv     = getnderiv(sLfd);

sbasismat = getbasismatrix(sevalarg, sbasisobj, snderiv);
if snderiv > 0 & ~strcmp(class(sLfd), 'double')
    derivwtmat = eval(sLfd, sevalarg);
    onerow <- ones(1,snbasis);
    for j = 1:snderiv
        if any(abs(derivwtmat(:,j))) > 1e-7
            sbasismat <- sbasismat +  ...
            (derivwtfdmat(:,j)*onerow) .* ...
                getbasismatrix(sevalarg, sbasisobj, j-1);
        end
    end
end

nt = length(tevalarg);
tbasisobj = gettbasis(bifd);
tnbasis   = getnbasis(tbasisobj);

tnderiv   = getnderiv(tLfd);

tbasismat = getbasismatrix(tevalarg, tbasisobj, tnderiv);
if tnderiv > 0 & ~strcmp(class(tLfd), 'double')
    derivwtmat = eval(tLfd, tevalarg);
    onerow <- ones(1,tnbasis);
    for j = 1:tnderiv
        if any(abs(derivwtmat(:,j))) > 1e-7
            tbasismat <- tbasismat +  ...
            (derivwtfdmat(:,j)*onerow) .* ...
                getbasismatrix(tevalarg, tbasisobj, j-1);
        end
    end
end

coef  = getcoef(bifd);
coefd = size(coef);
ndim  = length(coefd);

switch ndim
    case 2
        evalarray = sbasismat * coef * tbasismat';
    case 3
        ncurves  = coefd(3);
        evalarray = zeros(ns,nt,ncurves);
        for i= 1:ncurves
            evalarray(:,:,i) = sbasismat * coef(:,:,i) * tbasismat';
        end
    case 4
        ncurves  = coefd(3);
        nvar     = coefd(4);
        evalarray = zeros(ns,nt,ncurves,nvar);
        for i = 1:ncurves
            for j = 1:nvar
                evalarray(:,:,i,j) = sbasismat * coef(:,:,i,j) * tbasismat';
            end
        end
    otherwise
        error('coefficient array of improper dimension');
end
