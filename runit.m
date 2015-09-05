N=19;
nn=10;
bins = 200;

tic

%initialization
dosx=zeros(nn,bins,3); dosy = dosx; dosz = dosx;


Jx = 1;

flag = 1;
type = '++';
run_script
Izz = Ipp + Imm + Ipm;
%need to grab dos and average over these three takes.
dosx(:,:,1)=out(:,:,2); 
dosx(:,:,2)=out(:,:,3);
dosx(:,:,3)=out(:,:,4);

flag=0;
type = '--';
run_script
Iyz = Ipp + Imm + Ipm;
dosy(:,:,1)=out(:,:,2); 
dosy(:,:,2)=out(:,:,3);
dosy(:,:,3)=out(:,:,4);

flag = 0;
type = '+-';
run_script
II = Ipp + Imm + Ipm;
dosz(:,:,1)=out(:,:,2); 
dosz(:,:,2)=out(:,:,3);
dosz(:,:,3)=out(:,:,4);

%Use all the densities of states to get best estimates
dospt  = cat(1 ,dosx(:,:,1),dosy(:,:,1),dosz(:,:,1) );
dosmt  = cat(1, dosx(:,:,2), dosy(:,:,2), dosz(:,:,2) );
dospmt = cat(1, dosx(:,:,3), dosy(:,:,3), dosz(:,:,3) );

%dosperr = std(dospt,0,1);
%dosmerr = std(dosmt,0,1);
%dospmerr= std(dospmt,0,1);
doserr = [std(dospt,0,1).', std(dosmt,0,1).', std(dospmt,0,1).'];

disp(mean(doserr))
disp(max (doserr))

dosp = mean(dospt,1);
dosm = mean(dosmt,1);
dospm= mean(dospmt,1);

%Plot DOS
position = [pixs(3)/2   20       pixs(3)/2   pixs(4)/2-50];
h=figure('Position',position);
hold on;
plot(Ev/2,dosp,Ev/2,dosm,Ev/2,dospm);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['DOS for HyperHoneycomb Kitaev spinons: bins=',num2str(bins),', N=',num2str(N)])
xlabel('E');
ylabel('I(E)');
hold off;
filename = ['3D_DOS_10^7_Jx_',num2str(round(100*Jx)),'over100','_2'];
savefig(filename)
print(h, '-dpng', filename);

%Calculate the other 4 Raman, which are symmetry-related up to a factor
%Iyy = Ixx/4;
%Ixy = Ixx/2;
%Iyz = Ixz/2;
%Izz = Ixx/(Jx)^2 ?

%Plot 3 or 6 ramans
%position = [0   0       pixs(3)/2   pixs(4)/2];
h=figure;%('Position',position);
hold on;
plot(Ev,Izz,Ev,Iyz,Ev,II);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for HyperHoneycomb Kitaev spinons: bins=',num2str(bins),', N=',num2str(N)])
xlabel('E');
ylabel('I(E)');
hold off;
filename = ['3D_Raman_10^7_Jx_',num2str(round(100*Jx)),'over100','_2'];
savefig(filename)
print(h, '-dpng', filename);

%disp( mean(Ixx(Izz>0.1)./Izz(Izz>0.1)) )


for Jx = [1.43,0.3]
    
    dosx=zeros(nn,bins,3); dosy = dosx; dosz = dosx;

flag = 0;
type = '++';
run_script
Izz = Ipp + Imm + Ipm;
%need to grab dos and average over these three takes.
dosx(:,:,1)=out(:,:,2); 
dosx(:,:,2)=out(:,:,3);
dosx(:,:,3)=out(:,:,4);

flag=0;
type = '--';
run_script
Iyz = Ipp + Imm + Ipm;
dosy(:,:,1)=out(:,:,2); 
dosy(:,:,2)=out(:,:,3);
dosy(:,:,3)=out(:,:,4);

flag = 0;
type = '+-';
run_script
II = Ipp + Imm + Ipm;
dosz(:,:,1)=out(:,:,2); 
dosz(:,:,2)=out(:,:,3);
dosz(:,:,3)=out(:,:,4);

%Use all the densities of states to get best estimates
dospt  = cat(1 ,dosx(:,:,1),dosy(:,:,1),dosz(:,:,1) );
dosmt  = cat(1, dosx(:,:,2), dosy(:,:,2), dosz(:,:,2) );
dospmt = cat(1, dosx(:,:,3), dosy(:,:,3), dosz(:,:,3) );

%dosperr = std(dospt,0,1);
%dosmerr = std(dosmt,0,1);
%dospmerr= std(dospmt,0,1);
doserr = [std(dospt,0,1).', std(dosmt,0,1).', std(dospmt,0,1).'];

disp(mean(doserr))
disp(max (doserr))

dosp = mean(dospt,1);
dosm = mean(dosmt,1);
dospm= mean(dospmt,1);

%Plot DOS
position = [pixs(3)/2   20       pixs(3)/2   pixs(4)/2-50];
h=figure('Position',position);
hold on;
plot(Ev/2,dosp,Ev/2,dosm,Ev/2,dospm);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['DOS for HyperHoneycomb Kitaev: bins=',num2str(bins),', N=',num2str(N),' Jx=',num2str(Jx)])
xlabel('E');
ylabel('I(E)');
hold off;
filename = ['3D_DOS_10^7_Jx_',num2str(round(100*Jx)),'over100','_2'];
savefig(filename)
print(h, '-dpng', filename);

%Calculate the other 3 Raman, which are symmetry-related up to a factor
% Iyy = Ixx/4;
% Ixy = Ixx/2;
% Iyz = Ixz/2;

%Plot 3 or 6 ramans
%position = [0   0       pixs(3)/2   pixs(4)/2];
h=figure;%('Position',position);
hold on;
plot(Ev,Izz,Ev,Iyz,Ev,II);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman for HyperHoneycomb Kitaev: bins=',num2str(bins),', N=',num2str(N),' Jx=',num2str(Jx)])
xlabel('E');
ylabel('I(E)');
hold off;
filename = ['3D_Raman_10^7_Jx_',num2str(round(100*Jx)),'over100','_2'];
savefig(filename)
print(h, '-dpng', filename);

end


%plot the last one again for 5/8 of the width in E-space
maxEind = round(5/8*bins);

%Plot DOS
position = [pixs(3)/2   20       pixs(3)/2   pixs(4)/2-50];
h=figure('Position',position);
hold on;
plot(Ev(1:maxEind)/2,dosp(1:maxEind),Ev(1:maxEind)/2,dosm(1:maxEind),Ev(1:maxEind)/2,dospm(1:maxEind));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['DOS for HyperHoneycomb Kitaev: bins=',num2str(bins),', N=',num2str(N),' Jx=',num2str(Jx)])
xlabel('E');
ylabel('I(E)');
hold off;
filename = ['3D_DOS_10^7_Jx_',num2str(round(100*Jx)),'over100','_2'];
savefig(filename)
print(h, '-dpng', filename);


%Plot 3 or 6 ramans
%position = [0   0       pixs(3)/2   pixs(4)/2];
h=figure;%('Position',position);
hold on;
plot(Ev(1:maxEind),Izz(1:maxEind),Ev(1:maxEind),Iyz(1:maxEind),Ev(1:maxEind),II(1:maxEind));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman for HyperHoneycomb Kitaev: bins=',num2str(bins),', N=',num2str(N),' Jx=',num2str(Jx)])
xlabel('E');
ylabel('I(E)');
hold off;
filename = ['3D_Raman_10^7_Jx_',num2str(round(100*Jx)),'over100','_2'];
savefig(filename)
print(h, '-dpng', filename);


toc