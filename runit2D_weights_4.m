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


N=250;
nn=3;
bins = 200;%5*round(sqrt(N));

tic

tt = 50;

di = 3;

flag=1;
      
did = 0:20; s12=did;s22=did;Jxx =did; Jzz=did;th2=did;
for n = (1:21);
    th = (n-1)*pi/40; th2(n)=th*2/pi;
    Jx = 3*tan(th)/(1+2*tan(th));
    Jz = 3/(1+2*tan(th));
    Jxx(n) = Jx;
    Jzz(n) = Jz;    
   
    type = 'xx'; 
    run_script_2D_2 %Sets Ev, I , DE
    
    Jav = (Jzz+2*Jxx)/3;
    dE = Ev(200)/(200);
    
    s12(n) = sum(I)*dE;
    Ixx = I;
        
    if Jx ~= 1
        type = 'xy'; 
        run_script_2D_2
        s22(n) = sum(I)*dE;
        
        
        %Plot
        hh=figure;%('Position',position); hold on;
        plot(Ev,Ixx,Ev,I);
        %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
        title(['2D Raman: J_x=',num2str(Jx),', h=',num2str(h),', ',num2str(round(s12(n),di)),',',num2str(round(s22(n),di))])
        xlabel('\omega/J^z');
        ylabel('I');
        filename = ['2D_Raman3_Jx_',num2str(round(100*Jx)),'_h_',num2str(100*h)];
        saveas(hh,filename)
        print(hh, '-dpng', filename);
        print(hh, '-depsc', filename);
        
    else
        s22(n) = s12(n);
        
        
        %Plot
        hh=figure;%('Position',position); hold on;
        plot(Ev,Ixx);
        %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
        title(['2D Raman: J_x=',num2str(Jx),', h=',num2str(h)])
        xlabel('\omega/J^z');
        ylabel('I');
        filename = ['2D_Raman3_Jx_',num2str(round(100*Jx)),'_h_',num2str(100*h),', ',num2str(round(s12(n),di)),',',num2str(round(s22(n),di))];
        saveas(hh,filename)
        print(hh, '-dpng', filename);
        print(hh, '-depsc', filename);
        
    end
end
% I3 = runit14(.3,0);
% I4 = runit14(1.43,0);
% 
% Ih01 = runit14(1,.01);
% Ih03 = runit14(1,.03);
% Ih1  = runit14(1,.1);

%Ev = 2:18;
hh=figure;%('Position',position);
hold on;
plot(th2,s12./s22);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Rel. Raman Spectral Weights for 2D spinons \kappa=0'])
xlabel('(2/\pi)arctan(J_x/J_{av})');
ylabel('Rel. Spectral Weight');
legend({'I_{xx}/I_{xy}'}, 'Location', 'NorthWest');
hold off;
filename = ['2D_Raman_weights3'];
saveas(hh,filename)
print(hh, '-dpng', filename);


toc

diary 2D_Raman_weights_diary


