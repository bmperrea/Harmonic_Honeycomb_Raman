

I = runit12(1,0);
Ih = runit12(1,.1);

Ev = I{8}; Evh = Ih{8};

hh=figure;%('Position',position);
hold on;
plot(log(Ev),log(I{1}),log(Evh),log(Ih{1}));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for HyperHoneycomb spinons: h-field comparison (log plot)'])
xlabel('E');
ylabel('I(E)');
hold off;
filename = ['3D_Raman_10^7_magneticfield'];
savefig(filename)
print(hh, '-dpng', filename);





