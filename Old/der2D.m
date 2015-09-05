%Here I make derivatives plots for the 2D spectra

%sets defaults for the plots to look good
width = 5.1;     % Width in inches
height = 3;    % Height in inches
alw = 1;       % AxesLineWidth
fsz = 28;      % Fontsize
fna = 'Helvetica'; %Fontname
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize
interp = 'tex';

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
fig10 = load('2D_Raman_xx_10^7_Jx_100over100_2.fig','-mat');
fig14 = load('2D_Raman_xy_10^7_Jx_143over100_2.fig','-mat');
fig03 = load('2D_Raman_xy_10^7_Jx_30over100_2.fig','-mat');

d1=fig10.hgS_070000.children.children;
[d1,d2] = d1.properties;
Ev10  = d1.XData;
R10 = d1.YData;
%R102 = d2.YData;
bins = length(Ev10);
R10 = [diff(R10)*bins/Ev10(bins),0];

d1=fig14.hgS_070000.children.children;
[d1,d2] = d1.properties;
Ev14  = d1.XData;
R14 = d1.YData;
R142 = d2.YData;

R14 = [diff(R14)*bins/Ev14(bins),0];
R142 = [diff(R142)*bins/Ev14(bins),0];

d1=fig03.hgS_070000.children.children;
[d1,d2] = d1.properties;
Ev0  = d1.XData;
R03 = d1.YData;
R032 = d2.YData;

R03 = [diff(R03)*bins/Ev0(bins),0];
R032 = [diff(R032)*bins/Ev0(bins),0];

%hgS_070000
%hgM_070000

Jx=1; N=8000; type='xx';

%filter
a=1;
wS=3;
b = [.25,.5,.25];
%Esing = 114;
%b(114-wS:114+wS) = ones(1,2*wS+1); %avoids filtering the peak
%leng=size(R10(R10~=0),2);
%b(leng-ws,leng) = ones(1,wS+1);%avoids filtering the right end

R0 = filter(b,a,R10(R10~=0));
R0(size(R0,2))=0;

R4 = filter(b,a,R14(R14~=0));
R4(size(R4,2))=0;
R42 = filter(b,a,R142(R142~=0));
R42(size(R42,2))=0;

R3 = filter(b,a,R03(R03~=0));
R3(size(R3,2))=0;
R32 = filter(b,a,R032(R032~=0));
R32(size(R32,2))=0;

%restore peak height
[m,I] = max(R10);
R0(I) = m;

[m,I] = max(R14);
R4(I) = m;
[m,I] = max(R142);
R42(I) = m;
[m,I] = min(R142);
R42(I) = m;

%[m,I] = max(R03);
%R3(I) = m;
[m,I] = max(R032);
[m2,I2] = max(R32);
R32(I2) = m;


h=figure;%('Position',position);
hold on;
plot(Ev10(R10~=0),R0);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
%title(['Raman Spectrum Derivative for Honeycomb Kitaev spinons: bins = ',num2str(bins),', N= ',num2str(N)])
xlabel('\omega/J^z');
ylabel('I''(\omega)');
hold off;
filename = ['2D_Raman_10^7_Jx_',num2str(round(100*Jx)),'over100_der_2'];
savefig(filename)
print(h, '-dpng', filename);
print(h, '-depsc', filename);

Jx = 1.43;

h=figure;%('Position',position);
hold on;
plot(Ev14(R14~=0),R4,Ev14(R142~=0),R42);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
%title(['Raman Spectrum Derivative for Honeycomb Kitaev spinons: bins = ',num2str(bins),', N= ',num2str(N)])
xlabel('\omega/J^z');
ylabel('I''(\omega)');
%legend({'I_{xx}','I_{xy}'}, 'Location', 'SouthEast');
hold off;
filename = ['2D_Raman_10^7_Jx_',num2str(round(100*Jx)),'over100_der_2'];
savefig(filename)
print(h, '-dpng', filename);
print(h, '-depsc', filename);

Jx = 0.3;

h=figure;%('Position',position);
hold on;
plot(Ev0(R03~=0),R3,Ev0(R032~=0),R32);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
%title(['Raman Spectrum Derivative for Honeycomb Kitaev spinons: bins = ',num2str(bins),', N= ',num2str(N)])
xlabel('\omega/J^z');
ylabel('I''(\omega)');
%legend({'I_{xx}','I_{xy}'}, 'Location', 'SouthEast');
hold off;
filename = ['2D_Raman_10^7_Jx_',num2str(round(100*Jx)),'over100_der_2'];
savefig(filename)
print(h, '-dpng', filename);
print(h, '-depsc', filename);

