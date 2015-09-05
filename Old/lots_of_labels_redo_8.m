%Here I make derivatives plots for the 2D spectra

fsz = 15;
fsd = 19;

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
set(0,'defaultAxesFontSize',fsd);
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
%fig1 = load('3D_DOS_2_Jx_100_h_0_5.fig','-mat');
%fig2 = load('3D_DOS_2_Jx_100_h_0_2.fig','-mat');
 %fig3 = load('3D_Raman_2_Jx_30_h_0_28.fig','-mat');
 %fig4 = load('3D_Raman_2_Jx_100_h_0_28.fig','-mat');
 fig5 = load('3D_Raman_2_Jx_143_h_0_28.fig','-mat');
% fig6 = load('H1_Raman_Jx_100_h_0_22.fig','-mat');
% fig7 = load('H1_pi_Raman_Jx_100_h_0_22.fig','-mat');
% fig8 = load('H1_DOS_120_Jx_100_2.fig','-mat');
% fig9 = load('H1_pi_DOS_120_Jx_100_2.fig','-mat');

%for fig = [fig1,fig2,fig3,fig4,fig5,fig6,fig7,fig8,fig9]

Jx = 1; h =0; N=120;

%fig1
% d1=fig2.hgS_070000.children.children;
% [d1,d2,d3,d4] = d1.properties;
% Ev  = d1.XData;
% R1 = d1.YData;
% R2 = d2.YData;
% R3 = d3.YData;
% R4 = d4.YData;
% 
% hh=figure;%('Position',position);
% hold on;
% plot(Ev,R1,Ev,R2,Ev,R3,Ev,R4);
%      %title(['DOS for HyperHoneycomb Kitaev spinons: J_x=',num2str(Jx),', h=',num2str(h)])
%      xlabel('\omega/J^z');
%      ylabel('density of states');
%      legend({'\rho_{--}', '\rho_{++}', '\rho_{+-}','\rho_{total}'}, 'Location', 'NorthEast');
%      hold off;
%      filename = ['3D_DOS_2_Jx_',num2str(round(100*Jx)),'_h_',num2str(100*h),'_57'];
%      saveas(hh,filename)
%      print(hh, '-dpng', filename);
%      print(hh, '-depsc', filename);

%fig2
% d1=fig2.hgS_070000.children.children;
% d1 = d1.properties;
% Ev  = d1.XData;
% R1 = d1.YData;
% 
% hh=figure;%('Position',position);
% hold on;
% plot(Ev,R1);
%  %title(['DOS for HyperHoneycomb Kitaev spinons: J_x=',num2str(Jx),', h=',num2str(h)])
%  xlabel('\omega/J^z');
%  ylabel('DOS');
% % legend({'\rho_{--}', '\rho_{++}', '\rho_{+-}','\rho_{total}'}, 'Location', 'NorthEast');
%  hold off;
%  filename = ['3D_DOS_2_Jx_',num2str(round(100*Jx)),'_h_',num2str(100*h),'_27'];
%  saveas(hh,filename)
%  print(hh, '-dpng', filename);
%  print(hh, '-depsc', filename);

% %fig3
% Jx = 0.3;
% d1=fig3.hgS_070000.children.children;
% [d1,d2,d3] = d1.properties;
% Ev  = d1.XData;
% R1 = d1.YData;
% R2 = d2.YData;
% R3 = d3.YData;
% 
% hh=figure;%('Position',position);
% hold on;
% plot(Ev,R1,Ev,R2,Ev,R3);
% xlabel('\omega/J^z');
% ylabel('I(\omega)');
% legend({'I_{aa}','I_{ac}','I_{ab}'}, 'Location', 'NorthEast');
% hold off;
% filename = ['3D_Raman_2_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100),'_282'];
% saveas(hh,filename)
% print(hh, '-dpng', filename);
% print(hh, '-depsc', filename);
% 
% %fig4
% Jx = 1;
% d1=fig4.hgS_070000.children.children;
% [d1,d2,d3] = d1.properties;
% Ev  = d1.XData;
% R1 = d1.YData;
% R2 = d2.YData;
% R3 = d3.YData;
% 
% hh=figure;%('Position',position);
% hold on;
% plot(Ev,R1,Ev,R2,Ev,R3);
% xlabel('\omega/J^z');
% ylabel('I(\omega)');
% legend({'I_{aa}','I_{ac}','I_{ab}'}, 'Location', 'NorthEast');
% hold off;
% filename = ['3D_Raman_2_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100),'_282'];
% saveas(hh,filename)
% print(hh, '-dpng', filename);
% print(hh, '-depsc', filename);


%fig5
Jx=1.43;
d1=fig5.hgS_070000.children.children;
[d1,d2,d3] = d1.properties;
Ev  = d1.XData;
R1 = d1.YData;
R2 = d2.YData;
R3 = d3.YData;

hh=figure;%('Position',position);
hold on;
plot(Ev,R1,Ev,R2,Ev,R3);
xlabel('\omega/J^z');
ylabel('I(\omega)');
legend({'I_{aa}','I_{ac}','I_{ab}'}, 'Location', 'NorthEast');
hold off;
filename = ['3D_Raman_2_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100),'_282'];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);
% 
% %fig6
% Jx=1;
% d1=fig6.hgS_070000.children.children;
% [d1,d2,d3] = d1.properties;
% Ev  = d1.XData;
% R1 = d1.YData;
% R2 = d2.YData;
% R3 = d3.YData;
% 
% hh=figure;%('Position',position);
% hold on;
% plot(Ev,R1,Ev,R2,Ev,R3);
% xlabel('\omega/J^z');
% ylabel('I(\omega)');
% legend({'I_{aa}','I_{ac}','I_{ab}'}, 'Location', 'NorthEast');
% hold off;
% filename = ['H1_Raman_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100),'_227'];
% saveas(hh,filename)
% print(hh, '-dpng', filename);
% print(hh, '-depsc', filename);
% 
% %fig7
% d1=fig7.hgS_070000.children.children;
% [d1,d2,d3] = d1.properties;
% Ev  = d1.XData;
% R1 = d1.YData;
% R2 = d2.YData;
% R3 = d3.YData;
% 
% set(0,'defaultAxesFontSize',fsz);
% hh=figure;%('Position',position);
% hold on;
% plot(Ev,R1,Ev,R2,Ev,R3);
% xlabel('\omega/J^z');
% ylabel('I(\omega)');
% legend({'I_{aa}','I_{ac}','I_{ab}'}, 'Location', 'NorthEast');
% hold off;
% filename = ['H1_pi_Raman_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100),'_227'];
% saveas(hh,filename)
% print(hh, '-dpng', filename);
% print(hh, '-depsc', filename);
% 
% %fig8
% d1=fig8.hgS_070000.children.children;
% [d1,d2] = d1.properties;
% Ev  = d1.XData;
% R = d1.YData;
% 
% set(0,'defaultAxesFontSize',fsd);
% hh=figure;%('Position',position); hold on;
% plot(Ev,R);
% xlabel('\omega/J^z');
% ylabel('DOS');
% filename = ['H1_DOS_',num2str(N),'_Jx_',num2str(round(100*Jx)),'_27'];
% saveas(hh,filename)
% print(hh, '-dpng', filename);
% print(hh, '-depsc', filename);
% 
% %fig9
% d1=fig9.hgS_070000.children.children;
% [d1,d2] = d1.properties;
% Ev  = d1.XData;
% R = d1.YData;
% 
% %Plot DOS
% set(0,'defaultAxesFontSize',fsd);
% hh=figure;%('Position',position); hold on;
% plot(Ev,R);
% xlabel('\omega/J^z');
% ylabel('DOS');
% filename = ['H1_pi_DOS_',num2str(N),'_Jx_',num2str(round(100*Jx)),'_27'];
% saveas(hh,filename)
% print(hh, '-dpng', filename);
% print(hh, '-depsc', filename);


%hgS_070000
%hgM_070000