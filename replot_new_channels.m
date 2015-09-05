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



Ev = Ih_0{1}; Ip = Ih_0;
h = 0.1; Jx = 1; 

hh=figure;%('Position',position);
hold on;
plot(Ev,Ip{3},Ev,Ip{6},Ev,Ip{9},Ev,-Ip{5}/3,Ev,Ip{8}/9);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for H-1 spinons: J_x=',num2str(Jx),', h=',num2str(h)])
xlabel('\omega/J_z');
ylabel('I(\omega)');
legend({'I_{aa}','I_{ac}','I_{ab}','-I_{aa,cc}/3','I_{cc}/9'}, 'Location', 'NorthEast');
hold off;
filename = ['H0_Raman_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100),'_88'];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);



Ev = Ih2_0{1}; Ip = Ih2_0;
h = 0.03; Jx = 1;

hh=figure;%('Position',position);
hold on;
plot(Ev,Ip{3},Ev,Ip{6},Ev,Ip{9},Ev,-Ip{5}/3,Ev,Ip{8}/9);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for H-1 spinons: J_x=',num2str(Jx),', h=',num2str(h)])
xlabel('\omega/J_z');
ylabel('I(\omega)');
legend({'I_{aa}','I_{ac}','I_{ab}','-I_{aa,cc}/3','I_{cc}/9'}, 'Location', 'NorthEast');
hold off;
filename = ['H0_Raman_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100),'_88'];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);



Ev = I0_0{1}; Ip = I0_0;
h = 0; Jx = 1; 

hh=figure;%('Position',position);
hold on;
plot(Ev,Ip{3},Ev,Ip{6},Ev,Ip{9},Ev,-Ip{5}/3,Ev,Ip{8}/9);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for H-1 spinons: J_x=',num2str(Jx),', h=',num2str(h)])
xlabel('\omega/J_z');
ylabel('I(\omega)');
legend({'I_{aa}','I_{ac}','I_{ab}','-I_{aa,cc}/3','I_{cc}/9'}, 'Location', 'NorthEast');
hold off;
filename = ['H0_Raman_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100),'_88'];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);




Ip = Ih_64; Ev = Ip{1};
figure;hold on;
plot(Ev,Ip{3},Ev,-Ip{5}/3,Ev,Ip{8}/9,Ev,3*Ip{4},Ev,-Ip{7},Ev,Ip{6},Ev,Ip{9},Ev,Ip{10}/2);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
%title(['Raman Spectrum for H-1 spinons: J_x=',num2str(Jx),', h=',num2str(h)])
xlabel('\omega/J_z');
ylabel('I(\omega)');
legend({'I_{aa}','-I_{aa,cc}/3','I_{cc}/9','3*I_{aa,ac}','-I_{ac,cc}','I_{ac}','I_{ab}','I_{bc}/2'}, 'Location', 'NorthEast');
hold off;