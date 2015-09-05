%The first part of this code sets defaults for the plots to look good
%The second part runs code to make plots.

%initial data
width = 5.1;     % Width in inches
height = 3;    % Height in inches
alw = 1;       % AxesLineWidth
fsz = 18;      % Fontsize
fna = 'Helvetica'; %Fontname
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize
interp = 'tex';

%pixs = get(0,'screensize');
%width = pixs(3)/2; height = pixs(4)/2-50;

%Close the figures
close all;

% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz

set(0,'defaultAxesFontName',fna);
set(0,'defaultAxesFontSize',fsz);
set(0,'defaultTextInterpreter',interp);

% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);


%run code to make plots

%for Jxx = [0,.3,1.43]
   % runit_slab(Jxx,0);
%end

N=50; pbc =0;
N2 = 50;

h=.1; J=1;

Ih_64  = runit_hslab(64,pbc,round(N/sqrt(2)),J,h);
Ih_32  = runit_hslab(32,pbc,N,J,h);
Ih_16  = runit_hslab(16,pbc,round(sqrt(2)*N),J,h);
Ih_8   = runit_hslab(8 ,pbc,2*N,J,h);
Ih_4   = runit_hslab(4 ,pbc,round(2*sqrt(2)*N),J,h);
Ih_2   = runit_hslab(2 ,pbc,4*N,J,h);
Ih_1   = runit_hslab(1 ,pbc,round(4*sqrt(2)*N),J,h);

Ih_0 = runit14_8_f(N2,J,h);


%Ev = 1:18;
hh=figure;%('Position',position);
hold on;
plot(Ih_64{10},Ih_64{11},...
    Ih_32{10},Ih_32{11},Ih_16{10},Ih_16{11},...
Ih_8{10},Ih_8{11},Ih_4{10},Ih_4{11},Ih_2{10},Ih_2{11},...
Ih_1{10},Ih_1{11},Ih_0{10},Ih_0{11});
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['DOS for HH spinons: kappa = 0.1 J in a1 slab'])
xlabel('\omega/J^z');
ylabel('I_{ac}');
legend({'L=64','L=32','L=16','L=8','L=4','L=2','L=1','L=\infty'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_DOS_slab1_10_magneticfieldcomparison','_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

Ev = 1:18;
hh=figure;%('Position',position);
hold on;
plot(log(Ih_64{10}(Ev)),log(Ih_64{11}(Ev)),...
    log(Ih_32{10}(Ev)),log(Ih_32{11}(Ev)),log(Ih_16{10}(Ev)),log(Ih_16{11}(Ev)),...
    log(Ih_8{10}(Ev)),log(Ih_8{11}(Ev)),log(Ih_4{10}(Ev)),log(Ih_4{11}(Ev)),log(Ih_2{10}(Ev)),log(Ih_2{11}(Ev)),...
log(Ih_1{10}),log(Ih_1{11}),log(Ih_0{10}),log(Ih_0{11}));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['DOS for HH spinons: kappa = 0.1 J in a1 slab (log plot)'])
xlabel('log(\omega/J^z)');
ylabel('log(I_{ac})');
legend({'L=64','L=32','L=16','L=8','L=4','L=2','L=1','L=\infty'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_DOS_slab1_10_logmagneticfieldcomparison_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);
% 
% %The edge
% Ev = 1:18;
% hh=figure;%('Position',position);
% hold on;
% plot(log(Ih_64{10}(Ev)),log(Ih_64{11}(Ev,2)),...
%     log(Ih_32{10}(Ev)),log(Ih_32{11}(Ev,2)),log(Ih_16{10}(Ev)),log(Ih_16{11}(Ev,2)),...
%     log(Ih_8{10}(Ev)),log(Ih_8{11}(Ev,2)),log(Ih_4{10}(Ev)),log(Ih_4{11}(Ev,2)),log(Ih_2{10}(Ev)),log(Ih_2{11}(Ev,2)),...
% Ih_1{10},Ih_1{11},Ih_0{10},Ih_0{11});
% %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
% title(['Edge DOS for HH spinons: kappa = 0.1 J in a1 slab (log plot)'])
% xlabel('log(\omega/J^z)');
% ylabel('log(I_{ac}^{edge})');
% legend({'L=64','L=32','L=16','L=8','L=4','L=2','L=1','L=\infty'}, 'Location', 'SouthEast');
% hold off;
% filename = ['H0_DOS_slab1_10_logmagneticfieldcomparison_edge_p_',num2str(pbc)];
% savefig(filename)
% print(hh, '-dpng', filename);
% print(hh, '-depsc', filename);


%Ev = 1:18;
hh=figure;%('Position',position);
hold on;
plot(Ih_64{10},Ih_64{4},...
    Ih_32{10},Ih_32{4},Ih_16{10},Ih_16{4},...
Ih_8{10},Ih_8{4},Ih_4{10},Ih_4{4},Ih_2{10},Ih_2{4},...
Ih_1{10},Ih_1{4},Ih_0{10},Ih_0{4});
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for HH spinons: kappa = 0.1 J in a1 slab'])
xlabel('\omega/J^z');
ylabel('I_{ac}');
legend({'L=64','L=32','L=16','L=8','L=4','L=2','L=1','L=\infty'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_Raman_slab1_10_magneticfieldcomparison','_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

Ev = 1:18;
hh=figure;%('Position',position);
hold on;
plot(log(Ih_64{10}(Ev)),log(Ih_64{4}(Ev)),...
    log(Ih_32{10}(Ev)),log(Ih_32{4}(Ev)),log(Ih_16{10}(Ev)),log(Ih_16{4}(Ev)),...
    log(Ih_8{10}(Ev)),log(Ih_8{4}(Ev)),log(Ih_4{10}(Ev)),log(Ih_4{4}(Ev)),log(Ih_2{10}(Ev)),log(Ih_2{4}(Ev)),...
log(Ih_1{10}),log(Ih_1{4}),log(Ih_0{10}),log(Ih_0{4}));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for HH spinons: kappa = 0.1 J in a1 slab (log plot)'])
xlabel('log(\omega/J^z)');
ylabel('log(I_{ac})');
legend({'L=64','L=32','L=16','L=8','L=4','L=2','L=1','L=\infty'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_Raman_slab1_10_logmagneticfieldcomparison_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);


h=0.03;

Ih2_64  = runit_hslab(64,pbc,round(N/sqrt(2)),J,h);
Ih2_32  = runit_hslab(32,pbc,N,J,h);
Ih2_16  = runit_hslab(16,pbc,round(sqrt(2)*N),J,h);
Ih2_8   = runit_hslab(8 ,pbc,2*N,J,h);
Ih2_4   = runit_hslab(4 ,pbc,round(2*sqrt(2)*N),J,h);
Ih2_2   = runit_hslab(2 ,pbc,4*N,J,h);
Ih2_1   = runit_hslab(1 ,pbc,round(4*sqrt(2)*N),J,h);

Ih2_0 = runit14_8_f(N2,J,h);



%Ev = 1:18;
hh=figure;%('Position',position);
hold on;
plot(Ih2_64{10},Ih2_64{11},...
    Ih2_32{10},Ih2_32{11},Ih2_16{10},Ih2_16{11},...
Ih2_8{10},Ih2_8{11},Ih2_4{10},Ih2_4{11},Ih2_2{10},Ih2_2{11},...
Ih2_1{10},Ih2_1{11},Ih2_0{10},Ih2_0{11});
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['DOS for HH spinons: kappa = 0.03 J in a1 slab'])
xlabel('\omega/J^z');
ylabel('I_{ac}');
legend({'L=64','L=32','L=16','L=8','L=4','L=2','L=1','L=\infty'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_DOS_slab1_03_magneticfieldcomparison_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

Ev = 1:18;
hh=figure;%('Position',position);
hold on;
plot(log(Ih2_64{10}(Ev)),log(Ih2_64{11}(Ev)),...
    log(Ih2_32{10}(Ev)),log(Ih2_32{11}(Ev)),log(Ih2_16{10}(Ev)),log(Ih2_16{11}(Ev)),...
    log(Ih2_8{10}(Ev)),log(Ih2_8{11}(Ev)),log(Ih2_4{10}(Ev)),log(Ih2_4{11}(Ev)),log(Ih2_2{10}(Ev)),log(Ih2_2{11}(Ev)),...
log(Ih2_1{10}),log(Ih2_1{11}),log(Ih2_0{10}),log(Ih2_0{11}));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['DOS for HH spinons: kappa = 0.03 J in a1 slab (log plot)'])
xlabel('log(\omega/J^z)');
ylabel('log(I_{ac})');
legend({'L=64','L=32','L=16','L=8','L=4','L=2','L=1','L=\infty'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_DOS_slab1_03_logmagneticfieldcomparison_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);
% 
% %The edge
% Ev = 1:18;
% hh=figure;%('Position',position);
% hold on;
% plot(log(Ih2_64{10}(Ev)),log(Ih2_64{11}(Ev,2)),...
%     log(Ih2_32{10}(Ev)),log(Ih2_32{11}(Ev,2)),log(Ih2_16{10}(Ev)),log(Ih2_16{11}(Ev,2)),...
%     log(Ih2_8{10}(Ev)),log(Ih2_8{11}(Ev,2)),log(Ih2_4{10}(Ev)),log(Ih2_4{11}(Ev,2)),log(Ih2_2{10}(Ev)),log(Ih2_2{11}(Ev,2)));
% %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
% title(['Edge DOS for HH spinons: kappa = 0.03 J in a1 slab (log plot)'])
% xlabel('log(\omega/J^z)');
% ylabel('log(I_{ac}^{edge})');
% legend({'L=64','L=32','L=16','L=8','L=4','L=2','L=1','L=\infty'}, 'Location', 'SouthEast');
% hold off;
% filename = ['H0_DOS_slab1_03_logmagneticfieldcomparison_edge_p_',num2str(pbc)];
% savefig(filename)
% print(hh, '-dpng', filename);
% print(hh, '-depsc', filename);


%Ev = 1:18;
hh=figure;%('Position',position);
hold on;
plot(Ih2_64{10},Ih2_64{4},...
    Ih2_32{10},Ih2_32{4},Ih2_16{10},Ih2_16{4},...
Ih2_8{10},Ih2_8{4},Ih2_4{10},Ih2_4{4},Ih2_2{10},Ih2_2{4},...
Ih2_1{10},Ih2_1{4},Ih2_0{10},Ih2_0{4});
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for HH spinons: kappa = 0.03 J in a1 slab'])
xlabel('\omega/J^z');
ylabel('I_{ac}');
legend({'L=64','L=32','L=16','L=8','L=4','L=2','L=1','L=\infty'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_Raman_slab1_03_magneticfieldcomparison_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

Ev = 1:18;
hh=figure;%('Position',position);
hold on;
plot(log(Ih2_64{10}(Ev)),log(Ih2_64{4}(Ev)),...
    log(Ih2_32{10}(Ev)),log(Ih2_32{4}(Ev)),log(Ih2_16{10}(Ev)),log(Ih2_16{4}(Ev)),...
    log(Ih2_8{10}(Ev)),log(Ih2_8{4}(Ev)),log(Ih2_4{10}(Ev)),log(Ih2_4{4}(Ev)),log(Ih2_2{10}(Ev)),log(Ih2_2{4}(Ev)),...
log(Ih2_1{10}),log(Ih2_1{4}),log(Ih2_0{10}),log(Ih2_0{4}));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for HH spinons: kappa = 0.03 J in a1 slab (log plot)'])
xlabel('log(\omega/J^z)');
ylabel('log(I_{ac})');
legend({'L=64','L=32','L=16','L=8','L=4','L=2','L=1','L=\infty'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_Raman_slab1_03_logmagneticfieldcomparison_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);




h=0;


I0_64  = runit_hslab(64,pbc,round(N/sqrt(2)),J,h);
I0_32  = runit_hslab(32,pbc,N,J,h);
I0_16  = runit_hslab(16,pbc,round(sqrt(2)*N),J,h);
I0_8   = runit_hslab(8 ,pbc,2*N,J,h);
I0_4   = runit_hslab(4 ,pbc,round(2*sqrt(2)*N),J,h);
I0_2   = runit_hslab(2 ,pbc,4*N,J,h);
I0_1   = runit_hslab(1 ,pbc,round(4*sqrt(2)*N),J,h);

I0_0 = runit14_8_f(N2,J,h);



hh=figure;%('Position',position);
hold on;
plot(I0_64{10},I0_64{11},...
    I0_32{10},I0_32{11},I0_16{10},I0_16{11},...
I0_8{10},I0_8{11},I0_4{10},I0_4{11},I0_2{10},I0_2{11},...
I0_1{10},I0_1{11},I0_0{10},I0_0{11});
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['DOS for HH spinons: kappa = 0 in a1 slab'])
xlabel('\omega/J^z');
ylabel('I_{ac}');
legend({'L=64','L=32','L=16','L=8','L=4','L=2','L=1','L=\infty'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_DOS_slab1_0_magneticfieldcomparison_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

Ev = 1:18;
hh=figure;%('Position',position);
hold on;
plot(log(I0_64{10}(Ev)),log(I0_64{11}(Ev)),...
    log(I0_32{10}(Ev)),log(I0_32{11}(Ev)),log(I0_16{10}(Ev)),log(I0_16{11}(Ev)),...
    log(I0_8{10}(Ev)),log(I0_8{11}(Ev)),log(I0_4{10}(Ev)),log(I0_4{11}(Ev)),log(I0_2{10}(Ev)),log(I0_2{11}(Ev)),...
log(I0_1{10}),log(I0_1{11}),log(I0_0{10}),log(I0_0{11}));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['DOS for HH spinons: kappa = 0 in a1 slab (log plot)'])
xlabel('log(\omega/J^z)');
ylabel('log(I_{ac})');
legend({'L=64','L=32','L=16','L=8','L=4','L=2','L=1','L=\infty'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_DOS_slab1_0_logmagneticfieldcomparison_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

% 
% Ev = 1:18;
% hh=figure;%('Position',position);
% hold on;
% plot(log(I0_64{10}(Ev)),log(I0_64{11}(Ev,2)),...
%     log(I0_32{10}(Ev)),log(I0_32{11}(Ev,2)),log(I0_16{10}(Ev)),log(I0_16{11}(Ev,2)),...
%     log(I0_8{10}(Ev)),log(I0_8{11}(Ev,2)),log(I0_4{10}(Ev)),log(I0_4{11}(Ev,2)),log(I0_2{10}(Ev)),log(I0_2{11}(Ev,2)));
% %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
% title(['Edge DOS for HH spinons: kappa = 0 in a1 slab (log plot)'])
% xlabel('log(\omega/J^z)');
% ylabel('log(I_{ac}^{edge})');
% legend({'L=64','L=32','L=16','L=8','L=4','L=2','L=1','L=\infty'}, 'Location', 'SouthEast');
% hold off;
% filename = ['H0_DOS_slab1_0_logmagneticfieldcomparison_edge_p_',num2str(pbc)];
% savefig(filename)
% print(hh, '-dpng', filename);
% print(hh, '-depsc', filename);



hh=figure;%('Position',position);
hold on;
plot(I0_64{10},I0_64{4},...
    I0_32{10},I0_32{4},I0_16{10},I0_16{4},...
I0_8{10},I0_8{4},I0_4{10},I0_4{4},I0_2{10},I0_2{4},...
I0_1{10},I0_1{4},I0_0{10},I0_0{4});
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for HH spinons: kappa = 0 in a1 slab'])
xlabel('\omega/J^z');
ylabel('I_{ac}');
legend({'L=64','L=32','L=16','L=8','L=4','L=2','L=1','L=\infty'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_Raman_slab1_0_magneticfieldcomparison_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

Ev = 1:18;
hh=figure;%('Position',position);
hold on;
plot(log(I0_64{10}(Ev)),log(I0_64{4}(Ev)),...
    log(I0_32{10}(Ev)),log(I0_32{4}(Ev)),log(I0_16{10}(Ev)),log(I0_16{4}(Ev)),...
    log(I0_8{10}(Ev)),log(I0_8{4}(Ev)),log(I0_4{10}(Ev)),log(I0_4{4}(Ev)),log(I0_2{10}(Ev)),log(I0_2{4}(Ev)),...
log(I0_1{10}),log(I0_1{4}),log(I0_0{10}),log(I0_0{4}));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for HH spinons: kappa = 0 in a1 slab (log plot)'])
xlabel('log(\omega/J^z)');
ylabel('log(I_{ac})');
legend({'L=64','L=32','L=16','L=8','L=4','L=2','L=1','L=\infty'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_Raman_slab1_0_logmagneticfieldcomparison_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);


diary H0_slab_diary_a1
