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


%run code to make plots

%for Jxx = [0,.3,1.43]
   % runit_slab(Jxx,0);
%end

Ih_64  = runit_slab_64(1,.1);
Ih_32  = runit_slab_32(1,.1);
Ih_16  = runit_slab_16(1,.1);
Ih_8   = runit_slab_8(1,.1);
Ih_4   = runit_slab_4(1,.1);
Ih_2   = runit_slab_2(1,.1);
%Ih_1   = runit_slab_16(1,.1);



%Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(Ih_64{8},Ih_64{4}(:,1),...
    Ih_32{8},Ih_32{4}(:,1),Ih_16{8},Ih_16{4}(:,1),...
Ih_8{8},Ih_8{4}(:,1),Ih_4{8},Ih_4{4}(:,1),Ih_2{8},Ih_2{4}(:,1));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for HH spinons: kappa = 0.1 J in a slab'])
xlabel('\omega/J');
ylabel('I_{ac}');
legend({'L=64','L=32','L=16','L=8','L=4','L=2'}, 'Location', 'SouthEast');
hold off;
filename = ['3D_Raman_slab_10_magneticfieldcomparison'];
savefig(filename)
print(hh, '-dpng', filename);

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(Ih_64{8}(Ev)),log(Ih_64{4}(Ev,1)),...
    log(Ih_32{8}(Ev)),log(Ih_32{4}(Ev,1)),log(Ih_16{8}(Ev)),log(Ih_16{4}(Ev,1)),...
    log(Ih_8{8}(Ev)),log(Ih_8{4}(Ev,1)),log(Ih_4{8}(Ev)),log(Ih_4{4}(Ev,1)),log(Ih_2{8}(Ev)),log(Ih_2{4}(Ev,1)));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for HH spinons: kappa = 0.1 J in a slab (log plot)'])
xlabel('log(\omega/J)');
ylabel('log(I_{ac})');
legend({'L=64','L=32','L=16','L=8','L=4','L=2'}, 'Location', 'SouthEast');
hold off;
filename = ['3D_Raman_slab_10_logmagneticfieldcomparison'];
savefig(filename)
print(hh, '-dpng', filename);

%The edge
Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(Ih_64{8}(Ev)),log(Ih_64{4}(Ev,2)),...
    log(Ih_32{8}(Ev)),log(Ih_32{4}(Ev,2)),log(Ih_16{8}(Ev)),log(Ih_16{4}(Ev,2)),...
    log(Ih_8{8}(Ev)),log(Ih_8{4}(Ev,2)),log(Ih_4{8}(Ev)),log(Ih_4{4}(Ev,2)),log(Ih_2{8}(Ev)),log(Ih_2{4}(Ev,2)));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Edge Raman Spectrum for HH spinons: kappa = 0.1 J in a slab (log plot)'])
xlabel('log(\omega/J)');
ylabel('log(I_{ac}^{edge})');
legend({'L=64','L=32','L=16','L=8','L=4','L=2'}, 'Location', 'SouthEast');
hold off;
filename = ['3D_Raman_slab_10_logmagneticfieldcomparison_edge'];
savefig(filename)
print(hh, '-dpng', filename);




Ih2_64  = runit_slab_64(1,.03);
Ih2_32  = runit_slab_32(1,.03);
Ih2_16  = runit_slab_16(1,.03);
Ih2_8   = runit_slab_8(1,.03);
Ih2_4   = runit_slab_4(1,.03);
Ih2_2   = runit_slab_2(1,.03);
%Ih2_1   = runit_slab_16(1,.03);



%Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(Ih2_64{8},Ih2_64{4}(:,1),...
    Ih2_32{8},Ih2_32{4}(:,1),Ih2_16{8},Ih2_16{4}(:,1),...
Ih2_8{8},Ih2_8{4}(:,1),Ih2_4{8},Ih2_4{4}(:,1),Ih2_2{8},Ih2_2{4}(:,1));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for HH spinons: kappa = 0.03 J in a slab'])
xlabel('\omega/J');
ylabel('I_{ac}');
legend({'L=64','L=32','L=16','L=8','L=4','L=2'}, 'Location', 'SouthEast');
hold off;
filename = ['3D_Raman_slab_03_magneticfieldcomparison'];
savefig(filename)
print(hh, '-dpng', filename);

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(Ih2_64{8}(Ev)),log(Ih2_64{4}(Ev,1)),...
    log(Ih2_32{8}(Ev)),log(Ih2_32{4}(Ev,1)),log(Ih2_16{8}(Ev)),log(Ih2_16{4}(Ev,1)),...
    log(Ih2_8{8}(Ev)),log(Ih2_8{4}(Ev,1)),log(Ih2_4{8}(Ev)),log(Ih2_4{4}(Ev,1)),log(Ih2_2{8}(Ev)),log(Ih2_2{4}(Ev,1)));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for HH spinons: kappa = 0.03 J in a slab (log plot)'])
xlabel('log(\omega/J)');
ylabel('log(I_{ac})');
legend({'L=64','L=32','L=16','L=8','L=4','L=2'}, 'Location', 'SouthEast');
hold off;
filename = ['3D_Raman_slab_03_logmagneticfieldcomparison'];
savefig(filename)
print(hh, '-dpng', filename);

%The edge
Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(Ih2_64{8}(Ev)),log(Ih2_64{4}(Ev,2)),...
    log(Ih2_32{8}(Ev)),log(Ih2_32{4}(Ev,2)),log(Ih2_16{8}(Ev)),log(Ih2_16{4}(Ev,2)),...
    log(Ih2_8{8}(Ev)),log(Ih2_8{4}(Ev,2)),log(Ih2_4{8}(Ev)),log(Ih2_4{4}(Ev,2)),log(Ih2_2{8}(Ev)),log(Ih2_2{4}(Ev,2)));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Edge Raman Spectrum for HH spinons: kappa = 0.03 J in a slab (log plot)'])
xlabel('log(\omega/J)');
ylabel('log(I_{ac}^{edge})');
legend({'L=64','L=32','L=16','L=8','L=4','L=2'}, 'Location', 'SouthEast');
hold off;
filename = ['3D_Raman_slab_03_logmagneticfieldcomparison_edge'];
savefig(filename)
print(hh, '-dpng', filename);



I0_64  = runit_slab_64(1,0);
I0_32  = runit_slab_32(1,0);
I0_16  = runit_slab_16(1,0);
I0_8   = runit_slab_8(1,0);
I0_4   = runit_slab_4(1,0);
I0_2   = runit_slab_2(1,0);
%Ih_1   = runit_slab_16(1,.1);


hh=figure;%('Position',position);
hold on;
plot(I0_64{8},I0_64{4}(:,1),...
    I0_32{8},I0_32{4}(:,1),I0_16{8},I0_16{4}(:,1),...
I0_8{8},I0_8{4}(:,1),I0_4{8},I0_4{4}(:,1),I0_2{8},I0_2{4}(:,1));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for HH spinons: kappa = 0 in a slab'])
xlabel('\omega/J');
ylabel('I_{ac}');
legend({'L=64','L=32','L=16','L=8','L=4','L=2'}, 'Location', 'SouthEast');
hold off;
filename = ['3D_Raman_slab_0_magneticfieldcomparison'];
savefig(filename)
print(hh, '-dpng', filename);

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(I0_64{8}(Ev)),log(I0_64{4}(Ev,1)),...
    log(I0_32{8}(Ev)),log(I0_32{4}(Ev,1)),log(I0_16{8}(Ev)),log(I0_16{4}(Ev,1)),...
    log(I0_8{8}(Ev)),log(I0_8{4}(Ev,1)),log(I0_4{8}(Ev)),log(I0_4{4}(Ev,1)),log(I0_2{8}(Ev)),log(I0_2{4}(Ev,1)));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for HH spinons: kappa = 0 in a slab (log plot)'])
xlabel('log(\omega/J)');
ylabel('log(I_{ac})');
legend({'L=64','L=32','L=16','L=8','L=4','L=2'}, 'Location', 'SouthEast');
hold off;
filename = ['3D_Raman_slab_0_logmagneticfieldcomparison'];
savefig(filename)
print(hh, '-dpng', filename);


Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(I0_64{8}(Ev)),log(I0_64{4}(Ev,2)),...
    log(I0_32{8}(Ev)),log(I0_32{4}(Ev,2)),log(I0_16{8}(Ev)),log(I0_16{4}(Ev,2)),...
    log(I0_8{8}(Ev)),log(I0_8{4}(Ev,2)),log(I0_4{8}(Ev)),log(I0_4{4}(Ev,2)),log(I0_2{8}(Ev)),log(I0_2{4}(Ev,2)));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Edge Raman Spectrum for HH spinons: kappa = 0 in a slab (log plot)'])
xlabel('log(\omega/J)');
ylabel('log(I_{ac}^{edge})');
legend({'L=64','L=32','L=16','L=8','L=4','L=2'}, 'Location', 'SouthEast');
hold off;
filename = ['3D_Raman_slab_0_logmagneticfieldcomparison_edge'];
savefig(filename)
print(hh, '-dpng', filename);


diary 3D_slab_diary
