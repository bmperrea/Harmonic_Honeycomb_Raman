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
   % runit13(Jxx,0);
%end

%%%I0 = runit13(1,0);
%I3 = runit13(.3,0);
%I4 = runit13(1.43,0);

%Ih01 = runit13(1,.01);
%Ih03 = runit13(1,.03);
%%%Ih2  = runit13(1,.2);

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(I0{8}(Ev)),log(I0{4}(Ev)))
plot(log(Ih2{8}(Ev)),log(Ih2{4}(Ev)) , '--');
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
%title(['Raman Spectrum for HyperHoneycomb spinons: h-field comparison (log plot)_3'])
xlabel('log(\omega/J)');
ylabel('log(I_{ac})');
legend({'\kappa=0','\kappa=0.2'}, 'Location', 'SouthEast');
hold off;
filename = ['3D_Raman_10^7_magneticfield_4'];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc2', filename);
print(hh, '-djpeg', filename);

% y0 = log(I0{4}(2));
% x0 = log(I0{8}(2));
% m0 = ( log(I0{4}(10))-log(I0{4}(4)) )/( log(I0{8}(10))-log(I0{8}(4)) );
% 
% Ev = 2:18;
% hh=figure;%('Position',position);
% hold on;
% plot(log(I0{8}(Ev)),log(I0{4}(Ev)))
% plot(log(Ih2{8}(Ev)),log(Ih2{4}(Ev)) , '--');%,...
% %log(I0{8}(Ev)), y0 + (log(I0{8}(Ev)) - x0), log(Ih2{8}(Ev)), log(Ih2{4}(Ev)) + 1*(log(Ih2{8}(Ev)) - log(Ih2{8}(2))));
% %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
% %title(['Raman Spectrum for HyperHoneycomb spinons: h-field comparison (log plot)_3'])
% xlabel('log(\omega/J)');
% ylabel('log(I_{ac})');
% legend({'\kappa=0','\kappa=0.2'}, 'Location', 'SouthEast');
% hold off;
% filename = ['3D_Raman_10^7_magneticfield_42'];
% savefig(filename)
% print(hh, '-dpng', filename);


y0 = log(I0{7}(2));
x0 = log(I0{8}(2));
m0 = ( log(I0{7}(10))-log(I0{7}(4)) )/( log(I0{8}(10))-log(I0{8}(4)) );

Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(log(I0{8}(Ev)),log(I0{7}(Ev)))
plot(log(Ih2{8}(Ev)),log(Ih2{7}(Ev)) , '--')% , log(I0{7}(Ev)), y0 + (log(I0{7}(Ev)) - x0) );
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
%title(['Raman Spectrum for HyperHoneycomb spinons: h-field comparison (log plot)_3'])
xlabel('log(\omega/J)');
ylabel('log(DOS)');
legend({'\kappa=0','\kappa=0.2'}, 'Location', 'SouthEast');
hold off;
filename = ['3D_DOS_10^7_magneticfield_43'];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc2', filename);
print(hh, '-djpeg', filename);



hh=figure;%('Position',position); 
hold on;
plot(I0{8},I0{7});
plot(Ih2{8},Ih2{7},'--');
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
%title(['DOS for HyperHoneycomb Kitaev spinons: J_x=',num2str(Jx),', h=',num2str(h)])
xlabel('\omega/J');
ylabel('density of states');
legend({'\kappa=0','\kappa=0.2'}, 'Location', 'NorthEast');
%set(gca,'XTick',-3:3); 
%set(gca,'YTick',2*(0:5));
hold off;
filename = ['3D_DOS_10^7_4'];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc2', filename);
print(hh, '-djpeg', filename);



hh=figure;%('Position',position); 
hold on;
plot(I0{8},I0{4});
plot(Ih2{8},Ih2{4},'--');
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
%title(['DOS for HyperHoneycomb Kitaev spinons: J_x=',num2str(Jx),', h=',num2str(h)])
xlabel('\omega/J');
ylabel('I_{ac}');
legend({'\kappa=0','\kappa=0.2'}, 'Location', 'NorthEast');
%set(gca,'XTick',-3:3); 
%set(gca,'YTick',2*(0:5));
hold off;
filename = ['3D_Raman_10^7_4'];
savefig(filename)
print(hh, '-dpng', filename);
print(hh, '-depsc2', filename);
print(hh, '-djpeg', filename);



diary 3D_Raman_diary4
