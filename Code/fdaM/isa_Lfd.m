function isaLfd = isa_Lfd(Lfd)
%  ISA_LFD  checks argument is either integer or functional data object.

%  last modified 2 January 2003

isaLfd = 1;
if ~strcmp(class(Lfd), 'Lfd')
    isaLfd = 0;
end

