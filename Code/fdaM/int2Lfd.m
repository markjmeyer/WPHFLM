function Lfdobj = int2Lfd(m)
%INT2LFD converts a nonnegative integer to a linear differential
%  operator object that is equivalent to D^m.  The range of the
%  functional data object in any cell is set to [0,1], and is
%  not actually used when a linear differential operator object
%  of this nature is applied.  
%  In the event that m is already a linear differential operator
%  object, it returns the object immediately.  Thus, INT2LFD can
%  be used to screen whether an object is an integer or not.

%  Last modified 20 January 2003

%  check M

if isa_Lfd(m)
    Lfdobj = m;
    return;
end

if ~isnumeric(m)
    error(['Argument not numeric ', ...
           'and not a linear differential operator.']);
end

if length(m) ~= 1
    error('Argument is not a scalar.');
end
if round(m) ~= m
    error('Argument is not an integer.');
end
if m < 0
    error('Argument is negative.');
end

%  all the checks passed, set up a functional data object

Wfd0 = fd(0,create_constant_basis([0,1]));

%  set up WFDCELL, even if M = 0

wfdcell0{1} = Wfd0;

for j=2:m
    wfdcell0{j} = Wfd0;
end

afdcell0 = {};
ufdcell0 = {};

%  define the Lfd object

Lfdobj = Lfd(m, wfdcell0, afdcell0, ufdcell0);