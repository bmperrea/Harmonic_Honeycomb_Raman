%The first part of this code sets defaults for the plots to look good
%The second part runs code to make plots.

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


%run code to make plots

did = 0:20; s10=did;s4=did;s7=did;Jx =did;th2=did/10;
for n = (1:21);
    th = atan(th2(n)); %(n-1)*pi/40; th2(n)=th*2/pi;
    Jxx = 3*tan(th)/(1+2*tan(th));
    Jzz = 3/(1+2*tan(th));
    Jx(n) = Jxx;
    Jz(n) = Jzz;    
   
    I0_2 = runit14_82(Jxx,Jzz,0);
    Ev = I0_2{10};
    Jav = (Jzz+2*Jxx)/3;
    dE = Ev(200)/(200);
    
    s10(n) = sum(I0_2{1})*dE;
    s4(n) = sum(I0_2{4})*dE;
    s7(n) = sum(I0_2{7})*dE;
end

s10(1)=0;s4(1)=0;s7(1)=0;

%Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(th2,s10./s4,th2,s7./s4);
xlabel('J^x/J^z');
ylabel('Rel. Spectral Weight');
legend({'SW_{aa}/SW_{ac}','SW_{ab}/SW_{ac}'}, 'Location', 'NorthEast');
hold off;
filename = ['3D_Raman_weights38'];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);



%run code to make plots

did = 0:20; s10=did;s4=did;s7=did;Jx =did;th2=did/10;
for n = (1:21);
    th = atan(th2(n));%(n-1)*pi/40; th2(n)=th*2/pi;
    Jxx = 3*tan(th)/(1+2*tan(th));
    Jzz = 3/(1+2*tan(th));
    Jx(n) = Jxx;
    Jz(n) = Jzz;    
   
    I0_2 = runit14_82(Jxx,Jzz,0);
    Ev = I0_2{10};
    Jav = (Jzz+2*Jxx)/3;
    dE = Ev(200)/(200);
    
    s10(n) = sum(I0_2{1})*dE;
    s4(n) = sum(I0_2{4})*dE;
    s7(n) = sum(I0_2{7})*dE;
end
% I3 = runit14(.3,0);
% I4 = runit14(1.43,0);
% 
% Ih01 = runit14(1,.01);
% Ih03 = runit14(1,.03);
% Ih1  = runit14(1,.1);

s10(1)=0;s4(1)=0;s7(1)=0;

%Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(th2,s10./s4,th2,s7./s4);
xlabel('J^x/J^z');
ylabel('Rel. Spectral Weight');
legend({'SW^-_{aa}/SW^-_{ac}','SW^-_{ab}/SW^-_{ac}'}, 'Location', 'NorthEast');
hold off;
filename = ['3D_Raman_weights_lowband38'];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);


diary 3D_Raman_2_diary_double
