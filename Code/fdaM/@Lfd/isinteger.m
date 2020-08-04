function intwrd = isinteger(Lfdobj)
% ISINTEGER returns 1 of WFD and AFD are both zero functions.

%  Last modified 14 January 2003

%  check WFDCELL for emptyness or all zero

wfdcell = getwfd(Lfdobj);
wintwrd = 1;
if ~isempty(wfdcell)
    nderiv = Lfdobj.nderiv;
    for j=1:nderiv
        if any(getcoef(wfdcell{j}) ~= 0.0)
            wintwrd = 0;
        end
    end
end

%  check WFDCELL for emptyness or all zero

afdcell = getafd(Lfdobj);
aintwrd = 1;
if ~isempty(afdcell)
    k = max(size(afdcell));
    for i=1:k
        if any(getcoef(afdcell{i}) ~= 0.0)
            aintwrd = 0;
        end
    end
end

intwrd = wintwrd & aintwrd;