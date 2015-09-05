%Here I make derivatives plots for the 2D spectra

%initial data
width = 5.1;     % Width in inches
height = 3;    % Height in inches
alw = 1;       % AxesLineWidth
fsz = 15;      % Fontsize
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


%Get Data from fig files.
fig2D = load('2D_Raman_weights2.fig','-mat');
fig3D = load('3D_Raman_weights2.fig','-mat');

d1=fig2D.hgS_070000.children.children;
[d1,d2] = d1.properties;
Ev2D  = d1.XData;
R2D = d1.YData;
%R102 = d2.YData;
%bins = length(Ev10);
%R2D = [diff(R2D)*bins/Ev10(bins),0];

d1=fig3D.hgS_070000.children.children;
[d1,d2] = d1.properties;
Ev3D  = d1.XData;
R3D = d1.YData;
R3D2 = d2.YData;

%hgS_070000
%hgM_070000



hh=figure;%('Position',position);
hold on;
plot(Ev2D,R2D);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Rel. Raman Spectral Weights for 2D spinons \kappa=0'])
xlabel('J^x/J^z');
ylabel('Rel. Spectral Weight');
legend({'SW_{xx}/SW_{xy}'}, 'Location', 'NorthEast');
hold off;
filename = ['2D_Raman_weights3'];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);


hh=figure;%('Position',position);
hold on;
plot(Ev3D,R3D,Ev3D,R3D2);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
%title(['Relative Raman Spectral Weights for HyperHoneycomb spinons h=0'])
xlabel('J^x/J^z');
ylabel('Rel. Spectral Weight');
legend({'SW_{aa}/SW_{ac}','SW_{ab}/SW_{ac}'}, 'Location', 'NorthEast');
hold off;
filename = ['3D_Raman_weights3'];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);