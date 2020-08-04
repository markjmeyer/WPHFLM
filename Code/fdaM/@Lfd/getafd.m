function afdcell = getafd(Lfdobj)
%  GETAFD   Extracts the weight for the forcing function from LFDOBJ.

%  last modified 12 January 2003

if ~isa_Lfd(Lfdobj)
    error('Argument is not a linear differential operator object');
end

afdcell = Lfdobj.afdcell;


