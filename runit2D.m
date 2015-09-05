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




N=207;
nn=10;
bins = 200;%5*round(sqrt(N));

tic

%initialization
dosx=zeros(8,bins,3); dosy = dosx; dosz = dosx;


Jx = 1;

flag = 1;
type = 'xx';
run_script_2D
Ixx = I;
dosx=DE;

%Use all the densities of states to get best estimates
%dost = cat(1 ,dosx,dosxy,dosy );
%dost = cat(1 ,dosx,dosxy );
dos  = dosx;
%dos = mean(dost);

doserr = std(dos,0,1).';

disp(mean(doserr))
disp(max (doserr))

%Plot DOS

fsz = 24;

%position = [pixs(3)/2   20       pixs(3)/2   pixs(4)/2-50];
h=figure;%('Position',position);
hold on;
plot(Ev,dos);
title(['DOS for Honeycomb Kitaev spinons: bins = ',num2str(bins),', N= ',num2str(N)])
xlabel('E');
ylabel('density of states');
hold off;
filename = ['2D_DOS_10^7_Jx_',num2str(round(100*Jx)),'over100'];
savefig(filename)
print(h, '-dpng', filename);

fsz = 18;

%Plot 3 ramans
%position = [0   0       pixs(3)/2   pixs(4)/2];
h=figure;%('Position',position);
hold on;
plot(2*Ev,Ixx);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for Honeycomb Kitaev spinons: bins = ',num2str(bins),', N= ',num2str(N)])
xlabel('E');
ylabel('I(E)');
hold off;
filename = ['2D_Raman_',type,'_10^7_Jx_',num2str(round(100*Jx)),'over100'];
savefig(filename)
print(h, '-dpng', filename);


for Jx = [1.43,0.3]
    
    dosx=zeros(8,bins,3); dosy = dosx; dosz = dosx;

flag = 0;
type = 'xx';
run_script_2D
Ixx = I;
dosx=DE;

% flag=0;
% type = 'yy';
% run_script_2D
% Iyy = I;
% dosy=DE;


flag = 0;
type = 'xy';
run_script_2D
Ixy = I;
dosxy=DE;


%Use all the densities of states to get best estimates
%dost = cat(1 ,dosx,dosxy,dosy );
dost = cat(1 ,dosx,dosxy );
%dost  = dosx
dos = mean(dost);

doserr = std(dost,0,1).';

disp(mean(doserr))
disp(max (doserr))

%Plot DOS
%position = [pixs(3)/2   20       pixs(3)/2   pixs(4)/2-50];
h=figure;%('Position',position);
hold on;
plot(Ev,dos);
title(['DOS for Honeycomb Kitaev: bins = ',num2str(bins),', N= ',num2str(N),' Jx=',num2str(Jx)])
xlabel('E');
ylabel('density of states');
hold off;
filename=['2D_DOS_10^7_Jx_',num2str(round(100*Jx)),'over100'];
savefig(filename)
print(h, '-dpng', filename);

%Plot 3 ramans
%position = [0   0       pixs(3)/2   pixs(4)/2];
h=figure;%('Position',position);
hold on;
plot(2*Ev,Ixx,2*Ev,Ixy);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman for Honeycomb Kitaev: bins=',num2str(bins),', N=',num2str(N),' Jx=',num2str(Jx)])
xlabel('E');
ylabel('I(E)');
hold off;
filename = ['2D_Raman_',type,'_10^7_Jx_',num2str(round(100*Jx)),'over100'];
savefig(filename)
print(h, '-dpng', filename);

end

toc