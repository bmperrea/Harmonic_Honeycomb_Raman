N=28;
nn=10;
bins = 200;

tic

for Jx = [1,1.43,.3]

flag=0;
type = 'yy';
run_script2
Iyy = Ipp + Imm + Ipm;

type = 'xx';
run_script2
Ixx = Ipp + Imm + Ipm;

type = 'xy';
run_script2
Ixy = Ipp + Imm + Ipm;

type = '++';
run_script2
II = Ipp + Imm + Ipm;

type = '--';
run_script2
III = Ipp + Imm + Ipm;


% h=figure; hold on; plot(Ev,II,Ev,Ixx*9/4,Ev,Iyy*9,Ev,Izz)
% title(['Ramans ++,xx,yy,zz scaled for comparison ', 'N=',num2str(N),' Jx=',num2str(Jx)])
% hold off;
% 
% h=figure; hold on; plot(Ev,III,Ev,Iyz)
% title(['Ramans --,yz scaled for comparison ', 'N=',num2str(N),' Jx=',num2str(Jx)])
% hold off;

h=figure; hold on; 
plot(Ev(II>0),II(II>0)./Ixx(II>0),Ev(II>0),II(II>0)./Iyy(II>0),Ev(II>0),II(II>0)./Ixy(II>0))
title(['Ramans xx,yy,zz each over ++ ', 'N=',num2str(N),' Jx=',num2str(Jx)])
hold off;

h=figure; hold on; plot(Ev(III>0),Iyz(III>0)./III(III>0))
title(['Ramans yz over -- ', 'N=',num2str(N),' Jx=',num2str(Jx)])
hold off;

%print some results to the screen
disp(Jx)
disp( [II(II>0)./mean(Ixx(II>0)),mean(II(II>0)./Iyy(II>0)),mean(II(II>0)./Izz(II>0))] )
disp ( mean(III(III>0)./Iyz(III>0)) )

end

toc