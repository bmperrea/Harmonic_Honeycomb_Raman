%The first part of this code sets defaults for the plots to look good
%The second part runs code to make plots.

%initial data
width = 5.1;     % Width in inches
height = 3;    % Height in inches
alw = 1;       % AxesLineWidth
fsz = 13;      % Fontsize
fszd = 13;
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
set(0,'defaultAxesFontSize',fszd);
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
N2 = 80;

h=.1; J=1;

% Ih_64  = runit_hslab_8(64,pbc,round(N/sqrt(2)),J,h);
% Ih_32  = runit_hslab_8(32,pbc,N,J,h);
% Ih_16  = runit_hslab_8(16,pbc,round(sqrt(2)*N),J,h);
% Ih_8   = runit_hslab_8(8 ,pbc,2*N,J,h);
% Ih_4   = runit_hslab_8(4 ,pbc,round(2*sqrt(2)*N),J,h);
% Ih_2   = runit_hslab_8(2 ,pbc,4*N,J,h);
% Ih_1   = runit_hslab_8(1 ,pbc,round(4*sqrt(2)*N),J,h);
% 
% Ih_0 = runit14_8_f(64,1,round(N/sqrt(2)),J,h);

%Ih_0 = runit14_8_ac(N2,J,0.1);
%I0_0 = runit14_8_ac(N2,J,0);
%Ih2_0 = runit14_8_ac(N2,J,0.03);


%Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(Ih_4{1}/2,Ih_4{2},Ih_8{1}/2,Ih_8{2},Ih_16{1}/2,Ih_16{2},Ih_32{1}/2,Ih_32{2},Ih_64{1}/2,Ih_64{2} );
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
axis([0 7.5 0 inf])
title(['DOS for HH sinons: kk = 0.1 J in a1 slab'])
xlabel('\omega/J^z');
ylabel('DOS');
legend({'L=4','L=8','L=16','L=32','L=64'}, 'Location', 'NorthEast');
hold off;
filename = ['H0_DOS2_slab1_10_magneticfieldcomparison','_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(Ih_4{1}(Ev)/2),log(Ih_4{2}(Ev)),...
    log(Ih_8{1}(Ev)/2),log(Ih_8{2}(Ev)),log(Ih_16{1}(Ev)/2),log(Ih_16{2}(Ev)),...
    log(Ih_32{1}(Ev)/2),log(Ih_32{2}(Ev)),log(Ih_64{1}(Ev)/2),log(Ih_64{2}(Ev)) );
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['DOS for HH sinons: kk = 0.1 J in a1 slab lod blot'])
xlabel('log(\omega/J^z)');
ylabel('log(DOS)');
legend({'L=4','L=8','L=16','L=32','L=64'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_DOS2_slab1_10_logmagneticfieldcomparison_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

Ev1 = 12; Ev2 = 18;
disp(( log(Ih_64{2}(Ev2)) - log(Ih_64{2}(Ev1)) )/( log(Ih_64{1}(Ev2)) - log(Ih_64{1}(Ev1)) ))

% 
% %The edge
% Ev = 2:18;
% hh=figure;%('Position',position);
% hold on;
% plot(log(Ih_64{1}(Ev)),log(Ih_64{2}(Ev,2)),...
%     log(Ih_32{1}(Ev)),log(Ih_32{2}(Ev,2)),log(Ih_16{1}(Ev)),log(Ih_16{2}(Ev,2)),...
%     log(Ih_8{1}(Ev)),log(Ih_8{2}(Ev,2)),log(Ih_4{1}(Ev)),log(Ih_4{2}(Ev,2)),log(Ih_2{1}(Ev)),log(Ih_2{2}(Ev,2)),...
% Ih_1{1},Ih_1{2},Ih_0{1},Ih_0{2});
% %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
% title(['Edge DOS for HH sinons: kk = 0.1 J in a1 slab lod blot'])
% xlabel('log(\omega/J^z)');
% ylabel('log(I_{ac}^{edge})');
% legend({'L=64','L=32','L=16','L=8','L=4','L=2','L=1','L=\infty'}, 'Location', 'SouthEast');
% hold off;
% filename = ['H0_DOS_slab1_10_logmagneticfieldcomparison_edge_p_',num2str(pbc)];
% savefig(filename)
% print(hh, '-dpng', filename);
% print(hh, '-depsc', filename);
set(0,'defaultAxesFontSize',fsz);

%Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(Ih_4{1},Ih_4{11},Ih_8{1},Ih_8{11},Ih_16{1},Ih_16{11},Ih_32{1},Ih_32{11},Ih_64{1},Ih_64{11} );
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
axis([0 15 0 inf])
title(['Raman Sectrum for HH sinons: kk = 0.1 J in a1 slab'])
xlabel('\omega/J^z');
ylabel('I_{[ac]}');
legend({'L=4','L=8','L=16','L=32','L=64'}, 'Location', 'NorthEast');
hold off;
filename = ['H0_Raman2_slab1_10_magneticfieldcomparison','_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(Ih_4{1}(Ev)),log(Ih_4{11}(Ev)),...
    log(Ih_8{1}(Ev)),log(Ih_8{11}(Ev)),log(Ih_16{1}(Ev)),log(Ih_16{11}(Ev)),...
    log(Ih_32{1}(Ev)),log(Ih_32{11}(Ev)),log(Ih_64{1}(Ev)),log(Ih_64{11}(Ev)) );
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Sectrum for HH sinons: kk = 0.1 J in a1 slab lod blot'])
xlabel('log(\omega/J^z)');
ylabel('log(I_{[ac]})');
legend({'L=4','L=8','L=16','L=32','L=64'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_Raman2_slab1_10_logmagneticfieldcomparison_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);


disp(( log(Ih_64{11}(Ev2)) - log(Ih_64{11}(Ev1)) )/( log(Ih_64{1}(Ev2)) - log(Ih_64{1}(Ev1)) ))

%A different Raman spectrum for comparison

hh=figure;%('Position',position);
hold on;
plot(Ih_4{1},Ih_4{6},Ih_8{1},Ih_8{6},Ih_16{1},Ih_16{6},Ih_32{1},Ih_32{6},Ih_64{1},Ih_64{6} );
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
axis([0 15 0 inf])
title(['Raman Sectrum for HH sinons: kk = 0.1 J in a1 slab'])
xlabel('\omega/J^z');
ylabel('I_{ac}');
legend({'L=4','L=8','L=16','L=32','L=64'}, 'Location', 'NorthEast');
hold off;
filename = ['H0_Raman2_slab1_10_magneticfieldcomparison2','_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(Ih_4{1}(Ev)),log(Ih_4{6}(Ev)),...
    log(Ih_8{1}(Ev)),log(Ih_8{6}(Ev)),log(Ih_16{1}(Ev)),log(Ih_16{6}(Ev)),...
    log(Ih_32{1}(Ev)),log(Ih_32{6}(Ev)),log(Ih_64{1}(Ev)),log(Ih_64{6}(Ev)) );
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Sectrum for HH sinons: kk = 0.1 J in a1 slab lod blot'])
xlabel('log(\omega/J^z)');
ylabel('log(I_{ac})');
legend({'L=4','L=8','L=16','L=32','L=64'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_Raman2_slab1_10_logmagneticfieldcomparison2_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);


disp(( log(Ih_64{6}(Ev2)) - log(Ih_64{6}(Ev1)) )/( log(Ih_64{1}(Ev2)) - log(Ih_64{1}(Ev1)) ))


h=0.03;
% 
% Ih2_64  = runit_hslab_8(64,pbc,round(N/sqrt(2)),J,h);
% Ih2_32  = runit_hslab_8(32,pbc,N,J,h);
% Ih2_16  = runit_hslab_8(16,pbc,round(sqrt(2)*N),J,h);
% Ih2_8   = runit_hslab_8(8 ,pbc,2*N,J,h);
% Ih2_4   = runit_hslab_8(4 ,pbc,round(2*sqrt(2)*N),J,h);
% Ih2_2   = runit_hslab_8(2 ,pbc,4*N,J,h);
% Ih2_1   = runit_hslab_8(1 ,pbc,round(4*sqrt(2)*N),J,h);
% 
% Ih2_0 = runit_hslab_8(64,1,round(N/sqrt(2)),J,h);

set(0,'defaultAxesFontSize',fszd);

%Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(Ih2_4{1}/2,Ih2_4{2},Ih2_8{1}/2,Ih2_8{2},Ih2_16{1}/2,Ih2_16{2},Ih2_32{1}/2,Ih2_32{2},Ih2_64{1}/2,Ih2_64{2} );
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
axis([0 7.5 0 inf])
title(['DOS for HH sinons: kk = 0.03 J in a1 slab'])
xlabel('\omega/J^z');
ylabel('DOS');
legend({'L=4','L=8','L=16','L=32','L=64'}, 'Location', 'NorthEast');
hold off;
filename = ['H0_DOS2_slab1_03_magneticfieldcomparison_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(Ih2_4{1}(Ev)/2),log(Ih2_4{2}(Ev)),...
    log(Ih2_8{1}(Ev)/2),log(Ih2_8{2}(Ev)),log(Ih2_16{1}(Ev)/2),log(Ih2_16{2}(Ev)),...
    log(Ih2_32{1}(Ev)/2),log(Ih2_32{2}(Ev)),log(Ih2_64{1}(Ev)/2),log(Ih2_64{2}(Ev)) );
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['DOS for HH sinons: kk = 0.03 J in a1 slab lod blot'])
xlabel('log(\omega/J^z)');
ylabel('log(DOS)');
legend({'L=4','L=8','L=16','L=32','L=64'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_DOS2_slab1_03_logmagneticfieldcomparison_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);
% 
% %The edge
% Ev = 2:18;
% hh=figure;%('Position',position);
% hold on;
% plot(log(Ih2_64{1}(Ev)),log(Ih2_64{2}(Ev,2)),...
%     log(Ih2_32{1}(Ev)),log(Ih2_32{2}(Ev,2)),log(Ih2_16{1}(Ev)),log(Ih2_16{2}(Ev,2)),...
%     log(Ih2_8{1}(Ev)),log(Ih2_8{2}(Ev,2)),log(Ih2_4{1}(Ev)),log(Ih2_4{2}(Ev,2)),log(Ih2_2{1}(Ev)),log(Ih2_2{2}(Ev,2)));
% %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
% title(['Edge DOS for HH sinons: kk = 0.03 J in a1 slab lod blot'])
% xlabel('log(\omega/J^z)');
% ylabel('log(I_{ac}^{edge})');
% legend({'L=64','L=32','L=16','L=8','L=4','L=2','L=1','L=\infty'}, 'Location', 'SouthEast');
% hold off;
% filename = ['H0_DOS_slab1_03_logmagneticfieldcomparison_edge_p_',num2str(pbc)];
% savefig(filename)
% print(hh, '-dpng', filename);
% print(hh, '-depsc', filename);

disp(( log(Ih2_64{2}(Ev2)) - log(Ih2_64{2}(Ev1)) )/( log(Ih2_64{1}(Ev2)) - log(Ih2_64{1}(Ev1)) ))

set(0,'defaultAxesFontSize',fsz);

%Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(Ih2_4{1},Ih2_4{11},Ih2_8{1},Ih2_8{11},Ih2_16{1},Ih2_16{11},...
    Ih2_32{1},Ih2_32{11},Ih2_64{1},Ih2_64{11} );
axis([0 15 0 inf])
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Sectrum for HH sinons: kk = 0.03 J in a1 slab'])
xlabel('\omega/J^z');
ylabel('I_{[ac]}');
legend({'L=4','L=8','L=16','L=32','L=64'}, 'Location', 'NorthEast');
hold off;
filename = ['H0_Raman2_slab1_03_magneticfieldcomparison_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(Ih2_4{1}(Ev)),log(Ih2_4{11}(Ev)),...
    log(Ih2_8{1}(Ev)),log(Ih2_8{11}(Ev)),log(Ih2_16{1}(Ev)),log(Ih2_16{11}(Ev)),...
    log(Ih2_32{1}(Ev)),log(Ih2_32{11}(Ev)),log(Ih2_64{1}(Ev)),log(Ih2_64{11}(Ev)) );
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Sectrum for HH sinons: kk = 0.03 J in a1 slab lod blot'])
xlabel('log(\omega/J^z)');
ylabel('log(I_{[ac]})');
legend({'L=4','L=8','L=16','L=32','L=64'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_Raman2_slab1_03_logmagneticfieldcomparison_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);



hh=figure;%('Position',position);
hold on;
plot(Ih2_4{1},Ih2_4{6},Ih2_8{1},Ih2_8{6},Ih2_16{1},Ih2_16{6},Ih2_32{1},Ih2_32{6},Ih2_64{1},Ih2_64{6} );
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
axis([0 15 0 inf])
title(['Raman Sectrum for HH sinons: kk = 0.03 J in a1 slab'])
xlabel('\omega/J^z');
ylabel('I_{ac}');
legend({'L=4','L=8','L=16','L=32','L=64'}, 'Location', 'NorthEast');
hold off;
filename = ['H0_Raman2_slab1_03_magneticfieldcomparison2_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(Ih2_4{1}(Ev)),log(Ih2_4{6}(Ev)),... 
    log(Ih2_8{1}(Ev)),log(Ih2_8{6}(Ev)),log(Ih2_16{1}(Ev)),log(Ih2_16{6}(Ev)),...
    log(Ih2_32{1}(Ev)),log(Ih2_32{6}(Ev)),log(Ih2_64{1}(Ev)),log(Ih2_64{6}(Ev)) );
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Sectrum for HH sinons: kk = 0.03 J in a1 slab lod blot'])
xlabel('log(\omega/J^z)');
ylabel('log(I_{ac})');
legend({'L=4','L=8','L=16','L=32','L=64'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_Raman2_slab1_03_logmagneticfieldcomparison2_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

disp(( log(Ih2_64{11}(Ev2)) - log(Ih2_64{11}(Ev1)) )/( log(Ih2_64{1}(Ev2)) - log(Ih2_64{1}(Ev1)) ))
disp(( log(Ih2_64{6}(Ev2)) - log(Ih2_64{6}(Ev1)) )/( log(Ih2_64{1}(Ev2)) - log(Ih2_64{1}(Ev1)) ))


h=0;

% 
% I0_64  = runit_hslab_8(64,pbc,round(N/sqrt(2)),J,h);
% I0_32  = runit_hslab_8(32,pbc,N,J,h);
% I0_16  = runit_hslab_8(16,pbc,round(sqrt(2)*N),J,h);
% I0_8   = runit_hslab_8(8 ,pbc,2*N,J,h);
% I0_4   = runit_hslab_8(4 ,pbc,round(2*sqrt(2)*N),J,h);
% I0_2   = runit_hslab_8(2 ,pbc,4*N,J,h);
% I0_1   = runit_hslab_8(1 ,pbc,round(4*sqrt(2)*N),J,h);
% 
% I0_0  = runit_hslab_8(64,1,round(N/sqrt(2)),J,h);

set(0,'defaultAxesFontSize',fszd);

hh=figure;%('Position',position);
hold on;
plot(I0_4{1}/2,I0_4{2},I0_8{1}/2,I0_8{2},I0_16{1}/2,I0_16{2},I0_32{1}/2,I0_32{2},I0_64{1}/2,I0_64{2} );
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
axis([0 7.5 0 inf])
title(['DOS for HH sinons: kk = 0 in a1 slab'])
xlabel('\omega/J^z');
ylabel('DOS');
legend({'L=4','L=8','L=16','L=32','L=64'}, 'Location', 'NorthEast');
hold off;
filename = ['H0_DOS2_slab1_0_magneticfieldcomparison_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(I0_4{1}(Ev)/2),log(I0_4{2}(Ev)),...
    log(I0_8{1}(Ev)/2),log(I0_8{2}(Ev)),log(I0_16{1}(Ev)/2),log(I0_16{2}(Ev)),...
    log(I0_32{1}(Ev)/2),log(I0_32{2}(Ev)),log(I0_64{1}(Ev)/2),log(I0_64{2}(Ev)) );
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['DOS for HH sinons: kk = 0 in a1 slab lod blot'])
xlabel('log(\omega/J^z)');
ylabel('log(DOS)');
legend({'L=4','L=8','L=16','L=32','L=64'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_DOS2_slab1_0_logmagneticfieldcomparison_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

set(0,'defaultAxesFontSize',fsz);

disp(( log(I0_64{2}(Ev2)) - log(I0_64{2}(Ev1)) )/( log(I0_64{1}(Ev2)) - log(I0_64{1}(Ev1)) ))
% 
% Ev = 2:18;
% hh=figure;%('Position',position);
% hold on;
% plot(log(I0_64{1}(Ev)),log(I0_64{2}(Ev,2)),...
%     log(I0_32{1}(Ev)),log(I0_32{2}(Ev,2)),log(I0_16{1}(Ev)),log(I0_16{2}(Ev,2)),...
%     log(I0_8{1}(Ev)),log(I0_8{2}(Ev,2)),log(I0_4{1}(Ev)),log(I0_4{2}(Ev,2)),log(I0_2{1}(Ev)),log(I0_2{2}(Ev,2)));
% %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
% title(['Edge DOS for HH sinons: kk = 0 in a1 slab lod blot'])
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
plot(I0_4{1},I0_4{11},I0_8{1},I0_8{11},I0_16{1},I0_16{11},I0_32{1},I0_32{11},I0_64{1},I0_64{11} );
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
axis([0 15 0 inf])
title(['Raman Sectrum for HH sinons: kk = 0 in a1 slab'])
xlabel('\omega/J^z');
ylabel('I_{[ac]}');
legend({'L=4','L=8','L=16','L=32','L=64'}, 'Location', 'NorthEast');
hold off;
filename = ['H0_Raman2_slab1_0_magneticfieldcomparison_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);


Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(I0_4{1}(Ev)),log(I0_4{11}(Ev)),log(I0_8{1}(Ev)),log(I0_8{11}(Ev)),...
    log(I0_16{1}(Ev)),log(I0_16{11}(Ev)),log(I0_32{1}(Ev)),log(I0_32{11}(Ev)),log(I0_64{1}(Ev)),log(I0_64{11}(Ev)) );
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Sectrum for HH sinons: kk = 0 J in a1 slab lod blot'])
xlabel('log(\omega/J^z)');
ylabel('log(I_{[ac]})');
legend({'L=4','L=8','L=16','L=32','L=64'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_Raman2_slab1_0_logmagneticfieldcomparison_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

% hh=figure;%('Position',position);
% hold on;
% plot(I0_0{1},I0_0{2},...
%     Ih2_0{1},Ih2_0{2},Ih_0{1},Ih_0{2});
% %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
% title(['Raman Sectrum for HH sinons: TDL bbc'])
% xlabel('\omega/J^z');
% ylabel('I_{[ac]}');
% legend({'\kappa = 0','\kappa = 0.03','\kappa = 0.1'}, 'Location', 'SouthEast');
% hold off;
% filename = ['H0_Raman_magneticfieldcomparison_p_',num2str(pbc)];
% savefig(filename)
% print(hh, '-dpng', filename);
% print(hh, '-depsc', filename);
% 
% Ev = 2:18;
% hh=figure;%('Position',position);
% hold on;
% plot(log(I0_64{1}(Ev)),log(I0_64{11}(Ev)),...
%     log(I0_32{1}(Ev)),log(I0_32{11}(Ev)),log(I0_16{1}(Ev)),log(I0_16{11}(Ev)),...
%     log(I0_8{1}(Ev)),log(I0_8{11}(Ev)),log(I0_4{1}(Ev)),log(I0_4{11}(Ev)),log(I0_2{1}(Ev)),log(I0_2{11}(Ev)),...
% log(I0_1{1}(Ev)),log(I0_1{11}(Ev)));
% %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
% title(['Raman Sectrum for HH sinons: kk = 0 in a1 slab lod blot'])
% xlabel('log(\omega/J^z)');
% ylabel('log(I_{[ac]})');
% legend({'L=4','L=8','L=16','L=32','L=64'}, 'Location', 'SouthEast');
% hold off;
% filename = ['H0_Raman_slab1_0_logmagneticfieldcomparison_p_',num2str(pbc)];
% savefig(filename)
% print(hh, '-dpng', filename);
% print(hh, '-depsc', filename);

disp(( log(I0_64{11}(Ev2)) - log(I0_64{11}(Ev1)) )/( log(I0_64{1}(Ev2)) - log(I0_64{1}(Ev1)) ))


hh=figure;%('Position',position);
hold on;
plot(I0_4{1},I0_4{6},I0_8{1},I0_8{6},I0_16{1},I0_16{6},I0_32{1},I0_32{6},I0_64{1},I0_64{6} );
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
axis([0 15 0 inf])
title(['Raman Sectrum for HH sinons: kk = 0 in a1 slab'])
xlabel('\omega/J^z');
ylabel('I_{ac}');
legend({'L=4','L=8','L=16','L=32','L=64'}, 'Location', 'NorthEast');
hold off;
filename = ['H0_Raman2_slab1_0_magneticfieldcomparison2_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(I0_4{1}(Ev)),log(I0_4{6}(Ev)),...
    log(I0_8{1}(Ev)),log(I0_8{6}(Ev)),log(I0_16{1}(Ev)),log(I0_16{6}(Ev)),...
    log(I0_32{1}(Ev)),log(I0_32{6}(Ev)),log(I0_64{1}(Ev)),log(I0_64{6}(Ev)) );
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Sectrum for HH sinons: kk = 0 in a1 slab lod blot'])
xlabel('log(\omega/J^z)');
ylabel('log(I_{ac})');
legend({'L=4','L=8','L=16','L=32','L=64'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_Raman2_slab1_0_logmagneticfieldcomparison2_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

disp(( log(I0_64{6}(Ev2)) - log(I0_64{6}(Ev1)) )/( log(I0_64{1}(Ev2)) - log(I0_64{1}(Ev1)) ))



%Finally, some plots of the infinite limit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pbc = 1;

disp('%%%%%%%%%%%')

set(0,'defaultAxesFontSize',fszd);

hh=figure;%('Position',position);
hold on;
plot(I0_0{1}/2,I0_0{2},...
    Ih2_0{1}/2,Ih2_0{2},Ih_0{1}/2,Ih_0{2});
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
axis([0 7.5 0 inf])
title(['DOS for HH sinons: TDL bbc'])
xlabel('\omega/J^z');
ylabel('DOS');
legend({'\kappa = 0','\kappa = 0.03','\kappa = 0.1'}, 'Location', 'NorthEast');
hold off;
filename = ['H0_DOS_magneticfieldcomparison_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(I0_0{1}(Ev)/2),log(I0_0{2}(Ev)),...
    log(Ih2_0{1}(Ev)/2),log(Ih2_0{2}(Ev)),log(Ih_0{1}(Ev)/2),log(Ih_0{2}(Ev)));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['DOS for HH sinons: TDL bbc lod blot'])
xlabel('log(\omega/J^z)');
ylabel('log(DOS)');
legend({'\kappa = 0','\kappa = 0.03','\kappa = 0.1'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_DOS_logmagneticfieldcomparison_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

set(0,'defaultAxesFontSize',fsz);

disp(( log(I0_0{2}(Ev2)) - log(I0_0{2}(Ev1)) )/( log(I0_0{1}(Ev2)) - log(I0_0{1}(Ev1)) ))
disp(( log(Ih2_0{2}(Ev2)) - log(Ih2_0{2}(Ev1)) )/( log(Ih2_0{1}(Ev2)) - log(Ih2_0{1}(Ev1)) ))
disp(( log(Ih_0{2}(Ev2)) - log(Ih_0{2}(Ev1)) )/( log(Ih_0{1}(Ev2)) - log(Ih_0{1}(Ev1)) ))

hh=figure;%('Position',position);
hold on;
plot(I0_0{1},I0_0{11},...
    Ih2_0{1},Ih2_0{11},Ih_0{1},Ih_0{11});
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
axis([0 15 0 inf])
title(['Raman Sectrum for HH sinons: TDL bbc'])
xlabel('\omega/J^z');
ylabel('I_{[ac]}');
legend({'\kappa = 0','\kappa = 0.03','\kappa = 0.1'}, 'Location', 'NorthEast');
hold off;
filename = ['H0_Raman_magneticfieldcomparison_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(I0_0{1}(Ev)),log(I0_0{11}(Ev)),...
    log(Ih2_0{1}(Ev)),log(Ih2_0{11}(Ev)),log(Ih_0{1}(Ev)),log(Ih_0{11}(Ev)));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Sectrum for HH sinons: TDL bbc lod blot'])
xlabel('log(\omega/J^z)');
ylabel('log(I_{[ac]})');
legend({'\kappa = 0','\kappa = 0.03','\kappa = 0.1'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_Raman_logmagneticfieldcomparison_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

disp(( log(I0_0{11}(Ev2)) - log(I0_0{11}(Ev1)) )/( log(I0_0{1}(Ev2)) - log(I0_0{1}(Ev1)) ))
disp(( log(Ih2_0{11}(Ev2)) - log(Ih2_0{11}(Ev1)) )/( log(Ih2_0{1}(Ev2)) - log(Ih2_0{1}(Ev1)) ))
disp(( log(Ih_0{11}(Ev2)) - log(Ih_0{11}(Ev1)) )/( log(Ih_0{1}(Ev2)) - log(Ih_0{1}(Ev1)) ))


hh=figure;%('Position',position);
hold on;
plot(I0_0{1},I0_0{6},...
    Ih2_0{1},Ih2_0{6},Ih_0{1},Ih_0{6});
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
axis([0 15 0 inf])
title(['Raman Sectrum for HH sinons: TDL bbc'])
xlabel('\omega/J^z');
ylabel('I_{ac}');
legend({'\kappa = 0','\kappa = 0.03','\kappa = 0.1'}, 'Location', 'NorthEast');
hold off;
filename = ['H0_Raman_magneticfieldcomparison2_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(I0_0{1}(Ev)),log(I0_0{6}(Ev)),...
    log(Ih2_0{1}(Ev)),log(Ih2_0{6}(Ev)),log(Ih_0{1}(Ev)),log(Ih_0{6}(Ev)));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Sectrum for HH sinons: TDL bbc lod blot'])
xlabel('log(\omega/J^z)');
ylabel('log(I_{ac})');
legend({'\kappa = 0','\kappa = 0.03','\kappa = 0.1'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_Raman_logmagneticfieldcomparison2_p_',num2str(pbc)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

disp(( log(I0_0{6}(Ev2)) - log(I0_0{6}(Ev1)) )/( log(I0_0{1}(Ev2)) - log(I0_0{1}(Ev1)) ))
disp(( log(Ih2_0{6}(Ev2)) - log(Ih2_0{6}(Ev1)) )/( log(Ih2_0{1}(Ev2)) - log(Ih2_0{1}(Ev1)) ))
disp(( log(Ih_0{6}(Ev2)) - log(Ih_0{6}(Ev1)) )/( log(Ih_0{1}(Ev2)) - log(Ih_0{1}(Ev1)) ))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Now some comparisons between pbc and obc

%I0_64p runit_hslab_8(64,1,round(N/sqrt(2)),J,0);
%Ih_64p runit_hslab_8(64,1,round(N/sqrt(2)),J,0.1);
%Ih2_64p runit_hslab_8(64,1,round(N/sqrt(2)),J,0.03);


%set(0,'defaultAxesFontSize',fszd);

h=0;

hh=figure;%('Position',position);
hold on;
plot(I0_64p{1}/2,I0_64p{2},I0_64{1}/2,I0_64{2});
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
axis([0 7.5 0 inf])
title(['DOS for HH sinons'])
xlabel('\omega/J^z');
ylabel('DOS');
legend({'periodic','open'}, 'Location', 'NorthEast');
hold off;
filename = ['H0_DOS_po_h_',num2str(100*h)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(I0_64p{1}(Ev)/2),log(I0_64p{2}(Ev)/2),log(I0_64{1}(Ev)/2),log(I0_64{2}(Ev)/2));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['DOS for HH sinons lod'])
xlabel('log(\omega/J^z)');
ylabel('log(DOS)');
legend({'periodic','open'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_DOS_log_po_h_',num2str(100*h)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);


hh=figure;%('Position',position);
hold on;
plot(I0_64p{1},I0_64p{11},I0_64{1},I0_64{11});
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
axis([0 15 0 inf])
title(['Raman for HH sinons'])
xlabel('\omega/J^z');
ylabel('I_{[ac]}');
legend({'periodic','open'}, 'Location', 'NorthEast');
hold off;
filename = ['H0_Raman_po_h_',num2str(100*h)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(I0_64p{1}(Ev)),log(I0_64p{11}(Ev)),log(I0_64{1}(Ev)),log(I0_64{11}(Ev)));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman for HH sinons lod'])
xlabel('log(\omega/J^z)');
ylabel('log(I_{[ac]})');
legend({'periodic','open'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_Raman_log_po_h_',num2str(100*h)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);


h=0.1;

hh=figure;%('Position',position);
hold on;
plot(Ih_64p{1}/2,Ih_64p{2},Ih_64{1}/2,Ih_64{2});
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
axis([0 7.5 0 inf])
title(['DOS for HH sinons'])
xlabel('\omega/J^z');
ylabel('DOS');
legend({'periodic','open'}, 'Location', 'NorthEast');
hold off;
filename = ['H0_DOS_po_h_',num2str(100*h)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(Ih_64p{1}(Ev)/2),log(Ih_64p{2}(Ev)),log(Ih_64{1}(Ev)/2),log(Ih_64{2}(Ev)));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['DOS for HH sinons lod'])
xlabel('log(\omega/J^z)');
ylabel('log(DOS)');
legend({'periodic','open'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_DOS_log_po_h_',num2str(100*h)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);


hh=figure;%('Position',position);
hold on;
plot(Ih_64p{1},Ih_64p{11},Ih_64{1},Ih_64{11});
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
axis([0 15 0 inf])
title(['Raman for HH sinons'])
xlabel('\omega/J^z');
ylabel('I_{[ac]}');
legend({'periodic','open'}, 'Location', 'NorthEast');
hold off;
filename = ['H0_Raman_po_h_',num2str(100*h)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(Ih_64p{1}(Ev)),log(Ih_64p{11}(Ev)),log(Ih_64{1}(Ev)),log(Ih_64{11}(Ev)));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman for HH sinons lod'])
xlabel('log(\omega/J^z)');
ylabel('log(I_{[ac]})');
legend({'periodic','open'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_Raman_log_po_h_',num2str(100*h)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);


h=0.03;

hh=figure;%('Position',position);
hold on;
plot(Ih2_64p{1}/2,Ih2_64p{2},Ih2_64{1}/2,Ih2_64{2});
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
axis([0 7.5 0 inf])
title(['DOS for HH sinons'])
xlabel('\omega/J^z');
ylabel('DOS');
legend({'periodic','open'}, 'Location', 'NorthEast');
hold off;
filename = ['H0_DOS_po_h_',num2str(100*h)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(Ih2_64p{1}(Ev)/2),log(Ih2_64p{2}(Ev)),log(Ih2_64{1}(Ev)/2),log(Ih2_64{2}(Ev)));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['DOS for HH sinons lod'])
xlabel('log(\omega/J^z)');
ylabel('log(DOS)');
legend({'periodic','open'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_DOS_log_po_h_',num2str(100*h)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);


hh=figure;%('Position',position);
hold on;
plot(Ih2_64p{1},Ih2_64p{11},Ih2_64{1},Ih2_64{11});
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
axis([0 15 0 inf])
title(['Raman for HH sinons'])
xlabel('\omega/J^z');
ylabel('I_{[ac]}');
legend({'periodic','open'}, 'Location', 'NorthEast');
hold off;
filename = ['H0_Raman_po_h_',num2str(100*h)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(Ih2_64p{1}(Ev)),log(Ih2_64p{11}(Ev)),log(Ih2_64{1}(Ev)),log(Ih2_64{11}(Ev)));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman for HH sinons lod'])
xlabel('log(\omega/J^z)');
ylabel('log(I_{[ac]})');
legend({'periodic','open'}, 'Location', 'SouthEast');
hold off;
filename = ['H0_Raman_log_po_h_',num2str(100*h)];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);



diary H0_slab8_diary_a1_6
