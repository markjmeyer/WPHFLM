function basisequal = eq(basis1, basis2)
% EQ assesses whether two bases are equivalent.

%  Last modified 1 March 02

type1   = getbasistype(basis1);
range1  = getbasisrange(basis1);
nbasis1 = getnbasis(basis1);
pars1   = getbasispar(basis1);

type2   = getbasistype(basis2);
range2  = getbasisrange(basis2);
nbasis2 = getnbasis(basis2);
pars2   = getbasispar(basis2);

basisequal = 1;

if ~strcmp(type1,type2)
    basisequal = 0;
    return;
end

if range1(1) ~= range2(1) | range1(2) ~= range2(2)
    basisequal = 0;
    return;
end

if nbasis1 ~= nbasis2
    basisequal = 0;
    return;
end

if ~isequal(pars1, pars2)
    basisequal = 0;
    return;
end



