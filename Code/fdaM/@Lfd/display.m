function display(Lfd)
fprintf(['NDERIV = ', num2str(Lfd.nderiv),'\n']);
fprintf('\nWFD:\n');
display(Lfd.wfdcell);
fprintf('\nAFD:\n');
display(Lfd.afdcell);
fprintf('\nUFD:\n');
display(Lfd.ufdcell);

