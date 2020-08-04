function wfdcell = getwfd(Lfdobj)
%  GETWFD   Extracts the weight functions from LFDOBJ.

%  last modified 12 January 2003

if ~isa_Lfd(Lfdobj)
    error('Argument is not a linear differential operator object');
end

wfdcell = Lfdobj.wfdcell;


