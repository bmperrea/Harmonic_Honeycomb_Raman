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
%close all;

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

N=75; 
h=.15; Jx=1;

Ih1_00 = runit14_8_ac(N,Jx,h);
I = Ih1_00; Ev = I{1};

hh=figure;%('Position',position);

plot(Ev,I{3},Ev,-I{5}/3,Ev,I{8}/9,Ev,I{4},Ev,-I{7}/3,Ev,I{6},Ev,I{10}/2);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
%title(['Raman Spectrum for HyperHoneycomb Kitaev spinons: J_x=',num2str(Jx),', h=',num2str(h)])
hold on;
xlabel('\omega/J^z');
ylabel('I(\omega)');
legend({'I_{aa}','-I_{aa,cc}/3','I_{cc}/9','I_{aa,ac}','-I_{ac,cc}/3','I_{ac}','I_{bc}/2'}, 'Location', 'NorthEast');
hold off;
filename = ['3D_Raman_2_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100),'_38ac_ex'];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);
