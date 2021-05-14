clear;
clc;

load mat.out
a = spconvert(mat);
d = eig(full(a));
dlm = eigs(a,160);
[dum,ind] = sort(abs(d));
plot(dlm,'k+')
hold on
plot(d(ind(end-7:end)),'ks')
hold off
legend('eigs(a,160)','eig(full(a))')