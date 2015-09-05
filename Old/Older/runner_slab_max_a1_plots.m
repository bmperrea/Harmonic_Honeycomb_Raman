%The first part of this code sets defaults for the plots to look good
%The second part runs code to make plots.

%initial data
width = 5.1;     % Width in inches
height = 3;    % Height in inches
alw = 1;       % AxesLineWidth
fsz = 14;      % Fontsize
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


%Recover some \kappa = 0 plots from fig files

figh = load('3D_Raman_10^7_Jx_100_h_10_2.fig','-mat');
figh2 = load('3D_Raman_10^7_Jx_100_h_3_2.fig','-mat');
fig0 = load('3D_Raman_10^7_Jx_100_h_0_2.fig','-mat');

d1=figh.hgS_070000.children.children;
[d1,d2] = d1.properties;
Evh  = d1.XData;
%Rh = d1.YData;
Rh = d2.YData;
sh = sum(Rh);

d1=figh2.hgS_070000.children.children;
[d1,d2] = d1.properties;
Evh2  = d1.XData;
%R14 = d1.YData;
Rh2 = d2.YData;
sh2 = sum(Rh2);

d1=fig0.hgS_070000.children.children;
[d1,d2] = d1.properties;
Ev0  = d1.XData;
%R03 = d1.YData;
R0 = d2.YData;
s0 = sum(R0);



%run code to make plots

%for Jxx = [0,.3,1.43]
   % runit_slab(Jxx,0);
%end
% 
% Ih_64  = runit_slab_a1_64(1,.1);
% Ih_32  = runit_slab_a1_32(1,.1);
% Ih_16  = runit_slab_a1_16(1,.1);
% Ih_8   = runit_slab_a1_8(1,.1);
% Ih_4   = runit_slab_a1_4(1,.1);
% Ih_2   = runit_slab_a1_2(1,.1);
% %Ih_1   = runit_slab_16(1,.1);

sh_64 = sum(Ih_64{4});
sh_32 = sum(Ih_32{4});
sh_16 = sum(Ih_16{4});
sh_8 = sum(Ih_8{4});
sh_4 = sum(Ih_4{4});
sh_2 = sum(Ih_2{4});

%Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(Evh,Rh/sh,Ih_64{8},Ih_64{4}(:,1)/sh_64(1),...
    Ih_32{8},Ih_32{4}(:,1)/sh_32(1),Ih_16{8},Ih_16{4}(:,1)/sh_16(1),...
Ih_8{8},Ih_8{4}(:,1)/sh_8(1),Ih_4{8},Ih_4{4}(:,1)/sh_4(1),Ih_2{8},Ih_2{4}(:,1)/sh_2(1));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for HH spinons: kappa = 0.1 J in a1 slab'])
xlabel('\omega/J');
ylabel('I_{ac}');
legend({'L=\infty','L=64','L=32','L=16','L=8','L=4','L=2'}, 'Location', 'NorthEast');
hold off;
filename = ['3D_Raman_slab1_10_magneticfieldcomparison2'];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

% height = 4;    % Height in inches
% defsize = [left, bottom, width, height];
% set(0, 'defaultFigurePaperPosition', defsize);

daspect([1.5 1 1])

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(Evh(Ev)),log(Rh(Ev)/sh),log(Ih_64{8}(Ev)),log(Ih_64{4}(Ev,1)/sh_64(1)),...
    log(Ih_32{8}(Ev)),log(Ih_32{4}(Ev,1)/sh_32(1)),log(Ih_16{8}(Ev)),log(Ih_16{4}(Ev,1)/sh_16(1)),...
    log(Ih_8{8}(Ev)),log(Ih_8{4}(Ev,1)/sh_8(1)),log(Ih_4{8}(Ev)),log(Ih_4{4}(Ev,1)/sh_4(1)),log(Ih_2{8}(Ev)),log(Ih_2{4}(Ev,1)/sh_2(1)));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for HH spinons: kappa = 0.1 J in a1 slab (log plot)'])
xlabel('log(\omega/J)');
ylabel('log(I_{ac})');
legend({'L=\infty','L=64','L=32','L=16','L=8','L=4','L=2'}, 'Location', 'NorthWest');
hold off;
filename = ['3D_Raman_slab1_10_logmagneticfieldcomparison2'];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

% %The edge
% Ev = 2:18;
% hh=figure;%('Position',position);
% hold on;
% plot(log(Ih_64{8}(Ev)),log(Ih_64{4}(Ev,2)/sh_64(1)),...
%     log(Ih_32{8}(Ev)),log(Ih_32{4}(Ev,2)/sh_32(1)),log(Ih_16{8}(Ev)),log(Ih_16{4}(Ev,2)/sh_16(1)),...
%     log(Ih_8{8}(Ev)),log(Ih_8{4}(Ev,2)/sh_8(1)),log(Ih_4{8}(Ev)),log(Ih_4{4}(Ev,2)/sh_4(1)),log(Ih_2{8}(Ev)),log(Ih_2{4}(Ev,2)/sh_2(1)));
% %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
% title(['Edge Raman Spectrum for HH spinons: kappa = 0.1 J in a1 slab (log plot)'])
% xlabel('log(\omega/J)');
% ylabel('log(I_{ac}^{edge})');
% legend({'L=64','L=32','L=16','L=8','L=4','L=2'}, 'Location', 'NorthEast');
% hold off;
% filename = ['3D_Raman_slab1_10_logmagneticfieldcomparison_edge2'];
% savefig(filename)
% print(hh, '-dpng', filename);
% print(hh, '-depsc', filename);



% 
% Ih2_64  = runit_slab_a1_64(1,.03);
% Ih2_32  = runit_slab_a1_32(1,.03);
% Ih2_16  = runit_slab_a1_16(1,.03);
% Ih2_8   = runit_slab_a1_8(1,.03);
% Ih2_4   = runit_slab_a1_4(1,.03);
% Ih2_2   = runit_slab_a1_2(1,.03);
% %Ih2_1   = runit_slab_16(1,.03);
sh2_64 = sum(Ih2_64{4});
sh2_32 = sum(Ih2_32{4});
sh2_16 = sum(Ih2_16{4});
sh2_8 = sum(Ih2_8{4});
sh2_4 = sum(Ih2_4{4});
sh2_2 = sum(Ih2_2{4});

height = 3;    % Height in inches
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);

hh=figure;%('Position',position);
hold on;
plot(Evh2,Rh2/sh2,Ih2_64{8},Ih2_64{4}(:,1)/sh2_64(1),...
    Ih2_32{8},Ih2_32{4}(:,1)/sh2_32(1),Ih2_16{8},Ih2_16{4}(:,1)/sh2_16(1),...
Ih2_8{8},Ih2_8{4}(:,1)/sh2_8(1),Ih2_4{8},Ih2_4{4}(:,1)/sh2_4(1),Ih2_2{8},Ih2_2{4}(:,1)/sh2_2(1));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for HH spinons: kappa = 0.03 J in a1 slab'])
xlabel('\omega/J');
ylabel('I_{ac}');
legend({'L=\infty','L=64','L=32','L=16','L=8','L=4','L=2'}, 'Location', 'NorthEast');
hold off;
filename = ['3D_Raman_slab1_03_magneticfieldcomparison2'];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

height = 5;    % Height in inches
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(Evh2(Ev)),log(Rh2(Ev)/sh2),log(Ih2_64{8}(Ev)),log(Ih2_64{4}(Ev,1)/sh2_64(1)),...
    log(Ih2_32{8}(Ev)),log(Ih2_32{4}(Ev,1)/sh2_32(1)),log(Ih2_16{8}(Ev)),log(Ih2_16{4}(Ev,1)/sh2_16(1)),...
    log(Ih2_8{8}(Ev)),log(Ih2_8{4}(Ev,1)/sh2_8(1)),log(Ih2_4{8}(Ev)),log(Ih2_4{4}(Ev,1)/sh2_4(1)),log(Ih2_2{8}(Ev)),log(Ih2_2{4}(Ev,1)/sh2_2(1)));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for HH spinons: kappa = 0.03 J in a1 slab (log plot)'])
xlabel('log(\omega/J)');
ylabel('log(I_{ac})');
legend({'L=\infty','L=64','L=32','L=16','L=8','L=4','L=2'}, 'Location', 'NorthWest');
hold off;
filename = ['3D_Raman_slab1_03_logmagneticfieldcomparison2'];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

% %The edge
% Ev = 2:18;
% hh=figure;%('Position',position);
% hold on;
% plot(log(Ih2_64{8}(Ev)),log(Ih2_64{4}(Ev,2)),...
%     log(Ih2_32{8}(Ev)),log(Ih2_32{4}(Ev,2)),log(Ih2_16{8}(Ev)),log(Ih2_16{4}(Ev,2)),...
%     log(Ih2_8{8}(Ev)),log(Ih2_8{4}(Ev,2)),log(Ih2_4{8}(Ev)),log(Ih2_4{4}(Ev,2)),log(Ih2_2{8}(Ev)),log(Ih2_2{4}(Ev,2)));
% %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
% title(['Edge Raman Spectrum for HH spinons: kappa = 0.03 J in a1 slab (log plot)'])
% xlabel('log(\omega/J)');
% ylabel('log(I_{ac}^{edge})');
% legend({'L=64','L=32','L=16','L=8','L=4','L=2'}, 'Location', 'NorthEast');
% hold off;
% filename = ['3D_Raman_slab1_03_logmagneticfieldcomparison_edge'];
% savefig(filename)
% print(hh, '-dpng', filename);
% print(hh, '-depsc', filename);


% 
% I0_64  = runit_slab_a1_64(1,0);
% I0_32  = runit_slab_a1_32(1,0);
% I0_16  = runit_slab_a1_16(1,0);
% I0_8   = runit_slab_a1_8(1,0);
% I0_4   = runit_slab_a1_4(1,0);
% I0_2   = runit_slab_a1_2(1,0);
% %Ih_1   = runit_slab_16(1,.1);


s0_64 = sum(I0_64{4});
s0_32 = sum(I0_32{4});
s0_16 = sum(I0_16{4});
s0_8 = sum(I0_8{4});
s0_4 = sum(I0_4{4});
s0_2 = sum(I0_2{4});

height = 3;    % Height in inches
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);

hh=figure;%('Position',position);
hold on;
plot(Ev0,R0/s0,I0_64{8},I0_64{4}(:,1)/s0_64(1),...
    I0_32{8},I0_32{4}(:,1)/s0_32(1),I0_16{8},I0_16{4}(:,1)/s0_16(1),...
I0_8{8},I0_8{4}(:,1)/s0_8(1),I0_4{8},I0_4{4}(:,1)/s0_4(1),I0_2{8},I0_2{4}(:,1)/s0_2(1));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for HH spinons: kappa = 0 J in a1 slab'])
xlabel('\omega/J');
ylabel('I_{ac}');
legend({'L=\infty','L=64','L=32','L=16','L=8','L=4','L=2'}, 'Location', 'NorthEast');
hold off;
filename = ['3D_Raman_slab1_0_magneticfieldcomparison2'];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

height = 3;    % Height in inches
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(Ev0(Ev)),log(R0(Ev)/s0),log(I0_64{8}(Ev)),log(I0_64{4}(Ev,1)/s0_64(1)),...
    log(I0_32{8}(Ev)),log(I0_32{4}(Ev,1)/s0_32(1)),log(I0_16{8}(Ev)),log(I0_16{4}(Ev,1)/s0_16(1)),...
    log(I0_8{8}(Ev)),log(I0_8{4}(Ev,1)/s0_8(1)),log(I0_4{8}(Ev)),log(I0_4{4}(Ev,1)/s0_4(1)),log(I0_2{8}(Ev)),log(I0_2{4}(Ev,1)/s0_2(1)));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for HH spinons: kappa = 0 J in a1 slab (log plot)'])
xlabel('log(\omega/J)');
ylabel('log(I_{ac})');
legend({'L=\infty','L=64','L=32','L=16','L=8','L=4','L=2'}, 'Location', 'NorthWest');
hold off;
filename = ['3D_Raman_slab1_0_logmagneticfieldcomparison2'];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

% 
% Ev = 2:18;
% hh=figure;%('Position',position);
% hold on;
% plot(log(I0_64{8}(Ev)),log(I0_64{4}(Ev,2)),...
%     log(I0_32{8}(Ev)),log(I0_32{4}(Ev,2)),log(I0_16{8}(Ev)),log(I0_16{4}(Ev,2)),...
%     log(I0_8{8}(Ev)),log(I0_8{4}(Ev,2)),log(I0_4{8}(Ev)),log(I0_4{4}(Ev,2)),log(I0_2{8}(Ev)),log(I0_2{4}(Ev,2)));
% %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
% title(['Edge Raman Spectrum for HH spinons: kappa = 0 in a1 slab (log plot)'])
% xlabel('log(\omega/J)');
% ylabel('log(I_{ac}^{edge})');
% legend({'L=64','L=32','L=16','L=8','L=4','L=2'}, 'Location', 'NorthEast');
% hold off;
% filename = ['3D_Raman_slab1_0_logmagneticfieldcomparison_edge'];
% savefig(filename)
% print(hh, '-dpng', filename);
% print(hh, '-depsc', filename);


%diary 3D_slab_diary_a1_plots
