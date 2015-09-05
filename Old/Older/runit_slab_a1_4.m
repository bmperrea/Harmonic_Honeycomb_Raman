%This a script that calls the function that gives Raman spectra for the
%Kitaev Hyperhoneycomb (3D) model
function I = runit_slab(Jxx,h)

tic

disp(['Jxx= ',num2str(Jxx),' h= ',num2str(h)])

%m1 = {1,1,1,2,2,3};
%m2 = {1,2,3,2,3,3};
% 
N=400;
nn=2;
bins = 200;

%initialization
%dosx=zeros(nn,bins,3); dosy = dosx; dosz = dosx;
Dt = cell(1,nn); Iv = Dt; Ier = Dt; 
Iv1 = Iv; Iv2 = Iv; Ier1 = Iv; Ier2 = Iv;
It = cell(6,nn); It = It; 
Ivt1 = It; Ivt2 = It; It1 = It; It2 = It;
I = cell(1,8); Iermax = I; Iermean = I; I1=I; I2=I;
Iermax1 = I; Iermax2 = I; Iermean1=I; Iermean2=I;

%Runs for different parameters
Jz = 1;
flag = 1;


for Jx = [Jxx]

  %  if Jx == 1
  %  end
        
for n = 1:nn
    %The function gets run nn times so we can do statistics
    out = KitaevRaman_a1_4(N,bins,flag,nn,Jx,Jz,h);
    Ev = out{1}; %Same every time
    Dt{n}  = out{2}(:,1); %D = [DD,DD2]
    for m=1:6
        It{m,n} = out{m+2}; % I{m} = [Ipp{m},Imm{m},Ipm{m}];
    end
end

for m=1:6
    Ivtt = It(m,:);    
    Iv{m}  = mean(Ivtt,1,'c');
    Ier{m} = std(Ivtt,1,'c');
    %Ier2{m}= max(Ivtt,1,'c');
    Iermean{m} = mean(Ier{m});
    Iermax{m} = max(Ier{m});
    disp(Iermean{m})
    disp(Iermax{m})
    
    I{m} = Iv{m};
end

DD = mean(Dt,1,'c');
Der = std(Dt,1,'c');
Der2= max(Der);
Der3= mean(Der);
disp(Der3)
disp(Der2)

I{7} = DD;
I{8} = Ev;

pixs = get(0,'screensize');

%Plot DOS
%position = [pixs(3)/2   20       pixs(3)/2-114   pixs(4)/2-50];
hh=figure;%('Position',position); hold on;
plot(Ev,DD);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['DOS for Slab: J_x=',num2str(Jx),', h=',num2str(h)])
xlabel('\omega/J_z');
ylabel('density of states');
%set(gca,'XTick',-3:3); 
%set(gca,'YTick',2*(0:5));
hold off;
filename = ['3D_DOS_Slab_Jx_',num2str(round(100*Jx)),'_h_',num2str(100*h),'_2'];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

%Plot the ramans

hh=figure;%('Position',position);
hold on;
plot(Ev,I{1}(:,1),Ev,I{2}(:,1),Ev,I{3}(:,1),Ev,I{4}(:,1),Ev,I{5}(:,1),Ev,I{6}(:,1));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for Slab: J_x=',num2str(Jx),', h=',num2str(h)])
xlabel('\omega/J_z');
ylabel('I(\omega)');
legend({'I_{aa}', 'I_{aa,ac}', 'I_{aa,cc}','I_{ac}','I_{ac,cc}','I_{cc}'}, 'Location', 'NorthEast');
hold off;
filename = ['3D_Raman_Slab_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100)];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);


hh=figure;%('Position',position);
hold on;
plot(Ev,I{1}(:,1),Ev,-I{3}(:,1)/3,Ev,I{6}(:,1)/9);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for Slab: J_x=',num2str(Jx),', h=',num2str(h)])
xlabel('\omega/J_z');
ylabel('I(\omega)');
legend({'I_{aa}','-I_{aa,cc}/3','I_{cc}/9'}, 'Location', 'NorthEast');
hold off;
filename = ['3D_Raman_Slab_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100),'_3'];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);


hh=figure;%('Position',position);
hold on;
plot(Ev,I{1}(:,1),Ev,I{4}(:,1));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman Spectrum for Slab: J_x=',num2str(Jx),', h=',num2str(h)])
xlabel('\omega/J_z');
ylabel('I(\omega)');
legend({'I_{aa}','I_{ac}'}, 'Location', 'NorthEast');
hold off;
filename = ['3D_Raman_Slab_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100),'_2'];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

%Print the ratios
disp(mean(I{3}(I{1}(:,1)>0)./I{1}(I{1}(:,1)>0)))
disp(mean(I{6}(I{1}(:,1)>0)./I{1}(I{1}(:,1)>0)))
disp(mean(I{5}(I{2}(:,1)>0)./I{2}(I{2}(:,1)>0)))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Plot the ramans

hh=figure;%('Position',position);
hold on;
plot(Ev,I{1}(:,2),Ev,I{2}(:,2),Ev,I{3}(:,2),Ev,I{4}(:,2),Ev,I{5}(:,2),Ev,I{6}(:,2));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Edge1 Raman Spectrum for Slab: J_x=',num2str(Jx),', h=',num2str(h)])
xlabel('\omega/J_z');
ylabel('I(\omega)');
legend({'I_{aa}', 'I_{aa,ac}', 'I_{aa,cc}','I_{ac}','I_{ac,cc}','I_{cc}'}, 'Location', 'NorthEast');
hold off;
filename = ['3D_Raman_edge1_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100)];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);


hh=figure;%('Position',position);
hold on;
plot(Ev,I{1}(:,2),Ev,-I{3}(:,2)/3,Ev,I{6}(:,2)/9);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Edge1 Raman Spectrum for Slab: J_x=',num2str(Jx),', h=',num2str(h)])
xlabel('\omega/J_z');
ylabel('I(\omega)');
legend({'I_{aa}','-I_{aa,cc}/3','I_{cc}/9'}, 'Location', 'NorthEast');
hold off;
filename = ['3D_Raman_edge1_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100),'_3'];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);


hh=figure;%('Position',position);
hold on;
plot(Ev,I{1}(:,2),Ev,I{4}(:,2));
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Edge1 Raman Spectrum for Slab: J_x=',num2str(Jx),', h=',num2str(h)])
xlabel('\omega/J_z');
ylabel('I(\omega)');
legend({'I_{aa}','I_{ac}'}, 'Location', 'NorthEast');
hold off;
filename = ['3D_Raman_edge1_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100),'_2'];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

%Print the ratios
disp(mean(I{3}(I{1}(:,2)>0)./I{1}(I{1}(:,2)>0)))
disp(mean(I{6}(I{1}(:,2)>0)./I{1}(I{1}(:,2)>0)))
disp(mean(I{5}(I{2}(:,2)>0)./I{2}(I{2}(:,2)>0)))






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%Plot the ramans
% 
% hh=figure;%('Position',position);
% hold on;
% plot(Ev,I{1}(:,3),Ev,I{2}(:,3),Ev,I{3}(:,3),Ev,I{4}(:,3),Ev,I{5}(:,3),Ev,I{6}(:,3));
% %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
% title(['Edge2 Raman Spectrum for Slab: J_x=',num2str(Jx),', h=',num2str(h)])
% xlabel('\omega/J_z');
% ylabel('I(\omega)');
% legend({'I_{aa}', 'I_{aa,ac}', 'I_{aa,cc}','I_{ac}','I_{ac,cc}','I_{cc}'}, 'Location', 'NorthEast');
% hold off;
% filename = ['3D_Raman_edge2_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100)];
% saveas(hh,filename)
% print(hh, '-dpng', filename);
% 
% 
% hh=figure;%('Position',position);
% hold on;
% plot(Ev,I{1}(:,3),Ev,-I{3}(:,3)/3,Ev,I{6}(:,3)/9);
% %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
% title(['Edge2 Raman Spectrum for Slab: J_x=',num2str(Jx),', h=',num2str(h)])
% xlabel('\omega/J_z');
% ylabel('I(\omega)');
% legend({'I_{aa}','-I_{aa,cc}/3','I_{cc}/9'}, 'Location', 'NorthEast');
% hold off;
% filename = ['3D_Raman_edge2_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100),'_3'];
% saveas(hh,filename)
% print(hh, '-dpng', filename);
% 
% 
% hh=figure;%('Position',position);
% hold on;
% plot(Ev,I{1}(:,3),Ev,I{4}(:,3));
% %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
% title(['Edge2 Raman Spectrum for Slab: J_x=',num2str(Jx),', h=',num2str(h)])
% xlabel('\omega/J_z');
% ylabel('I(\omega)');
% legend({'I_{aa}','I_{ac}'}, 'Location', 'NorthEast');
% hold off;
% filename = ['3D_Raman_edge2_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100),'_2'];
% saveas(hh,filename)
% print(hh, '-dpng', filename);
% 
% %Print the ratios
% disp(mean(I{3}(I{1}(:,3)>0)./I{1}(I{1}(:,3)>0)))
% disp(mean(I{6}(I{1}(:,3)>0)./I{1}(I{1}(:,3)>0)))
% disp(mean(I{5}(I{2}(:,3)>0)./I{2}(I{2}(:,3)>0)))
% 


end


toc
end