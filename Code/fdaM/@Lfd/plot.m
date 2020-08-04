function plot(Lfd)
m       = Lfd.nderiv;
wfdcell = Lfd.wfdcell;
afdcell = Lfd.afdcell;
k       = length(afdcell);
nplot   = m + k;
for j=1:m
    subplot(nplot,1,j)
    plot(wfdcell{j})
    ylabel(['w_',num2str(j-1),'(t)'])
end
for j=1:k
    subplot(nplot,1,m+j)
    plot(afdcell{j})
    ylabel(['a_',num2str(j-1),'(t)'])
end