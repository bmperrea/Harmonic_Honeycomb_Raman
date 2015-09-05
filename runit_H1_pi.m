%This a script that calls the function that gives Raman spectra for the
%Kitaev Hyperhoneycomb (3D) model
function I = runit_H1(Jx,N)

tic

h=0;

disp(['Jxx= ',num2str(Jx),' h= ',num2str(h)])

%m1 = {1,1,1,2,2,3};
%m2 = {1,2,3,2,3,3};
% 
%N=30;
nn=3;
bins = 200;

%initialization
%dosx=zeros(nn,bins,3); dosy = dosx; dosz = dosx;
Dt = cell(1,nn); Iv = Dt; Ier = Dt; 
It = cell(6,nn); 
I = cell(1,8); Iermax = I; Iermean = I; 

%Runs for different parameters
Jz = 1;
flag = 1;



        
for n = 1:nn
    %The function gets run nn times so we can do statistics
    out = KitaevRaman_H1_pi(N,bins,flag,nn,Jx,Jz,h);
    Ev = out{1}; %Same every time
    DDt{n}  = out{2}; %D = [DD,DD2]
    Dt{n}  = out{3};
%     for m=1:6
%         It{m,n} = out{m+2}; % I{m} = [Ipp{m},Imm{m},Ipm{m}];
%     end
end
% 
% for m=1:6
%     Ivtt = It(m,:);    
%     Iv{m}  = mean(Ivtt,1,'c');
%     Ier{m} = std(Ivtt,1,'c');
%     %Ier2{m}= max(Ivtt,1,'c');
%     Iermean{m} = mean(Ier{m});
%     Iermax{m} = max(Ier{m});
%     disp(Iermean{m})
%     disp(Iermax{m})
%     
%     I{m} = Iv{m};
% end
DD = mean(DDt,1,'c');
DDer = std(DDt,1,'c');
DDer2= max(DDer);
DDer3= mean(DDer);
disp(DDer3)
disp(DDer2)

D = mean(Dt,1,'c');
Der = std(DDt,1,'c');
Der2= max(Der);
Der3= mean(Der);
disp(Der3)
disp(Der2)

I{7} = DD;
I{8} = D;
I{9} = Ev;


%Plot 2-particle DOS
%position = [pixs(3)/2   20       pixs(3)/2-114   pixs(4)/2-50];
hh=figure;%('Position',position); hold on;
plot(Ev,DD);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['2p DOS for H1 (\pi): J_x=',num2str(Jx),', h=',num2str(h)])
xlabel('\omega/J_z');
ylabel('density of states');
%set(gca,'XTick',-3:3); 
%set(gca,'YTick',2*(0:5));
hold off;
filename = ['3D_2p_DOS_H1_pi_Jx_',num2str(round(100*Jx)),'_h_',num2str(100*h)];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

%Plot DOS
%position = [pixs(3)/2   20       pixs(3)/2-114   pixs(4)/2-50];
hh=figure;%('Position',position); hold on;
plot(Ev/4,D);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['DOS for H1 (\pi): J_x=',num2str(Jx),', h=',num2str(h)])
xlabel('\omega/J_z');
ylabel('density of states');
%set(gca,'XTick',-3:3); 
%set(gca,'YTick',2*(0:5));
hold off;
filename = ['DOS_H1_pi_Jx_',num2str(round(100*Jx)),'_h_',num2str(100*h)];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

%Display the total energy with an error estimate

dE = Ev(200)/bins;
EE = (Ev/4).*D*dE/4;


for m = 1:nn
    EEE(m)=sum(Dt{m}.*(Ev/4)*dE/4);
end
    
    format long e
disp(mean(EEE));
    format short
disp(std(EEE));
%disp(sum(D*dE/4));


%Store the mean and std of the total energy
I{1}=mean(EEE);
I{2}=std(EEE);

toc
end