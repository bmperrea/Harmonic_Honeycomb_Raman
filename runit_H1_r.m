%This a script that calls the function that gives Raman spectra for the
%Kitaev Hyperhoneycomb (3D) model
function I = runit_H1_r2(Jxx,h,~)

tic

disp(['Jxx= ',num2str(Jxx),' h= ',num2str(h)])

N=50;
nn=3;
bins = 200;

%initialization
 
Dt = cell(1,nn); DDt=Dt; Iv = Dt; Ier = Dt; 
Ivt = cell(8,nn); It = Ivt;
I = cell(1,11); Iermax = I; Iermean = I;

%Runs for different parameters
Jz = 1;
flag = 1;


for Jx = Jxx

  %  if Jx == 1
  %  end
        
%initialization
Dt = cell(1,nn); Iv = Dt; Ier = Dt;
Ivt = cell(8,nn); 
I = cell(1,11); Iermax = I; Iermean = I;

%Runs for different parameters
Jz = 1;
flag = 1;


for Jx = Jxx

  %  if Jx == 1
  %  end
        
for n = 1:nn
    %The function gets run nn times so we can do statistics
    out = KitaevRaman_H1_r(N,bins,flag,nn,Jx,Jz,h);
    Ev = out{1}; %Same every time
    DDt{n}  = out{2};
    Dt{n} = out{3};
    for m=1:8
        Ivt{m,n} = out{m+3}; 
    end
end

for m=1:8
    Ivtt = Ivt(m,:);    
    Iv{m}  = mean(Ivtt,1,'c');
    Ier{m} = std(Ivtt,1,'c');
    Iermean{m} = mean(Ier{m});
    Iermax{m} = max(Ier{m});
    disp(Iermean{m})
    disp(Iermax{m})
    
    I{m} = sum(Iv{m},2);
end

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



dE = Ev(bins)/bins;
EE = (Ev/4).*D*dE/4;
for m = 1:nn
    EEE(m)=sum(Dt{m}.*(Ev/4)*dE/4);
end
 format long e
disp(mean(EEE));
    format short
disp(std(EEE));


I{9} = DD;
I{11} = D;
I{10} = Ev;
I{12} = [mean(EEE),std(EEE)];

fsz = 15;
fsd = 19;

%Plot DOS
set(0,'defaultAxesFontSize',fsd);
%position = [pixs(3)/2   20       pixs(3)/2-114   pixs(4)/2-50];
hh=figure;%('Position',position); hold on;
plot(Ev,DD);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
%title(['DOS for H-1 spinons: J_x=',num2str(Jx),', h=',num2str(h)])
xlabel('\omega/J_z');
ylabel('DOS');
%set(gca,'XTick',-3:3); 
%set(gca,'YTick',2*(0:5));
%hold off;
filename = ['H1_DOS_',num2str(N),'_Jx_',num2str(round(100*Jx))];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);


%Plot the ramans
set(0,'defaultAxesFontSize',fsz);
hh=figure;%('Position',position);
hold on;
plot(Ev,I{1},Ev,I{2},Ev,-I{3},Ev,I{4},Ev,I{5},Ev,I{6}/3,Ev,I{7},Ev,I{8}/2);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
%title(['Raman Spectrum for H-1 spinons: J_x=',num2str(Jx),', h=',num2str(h)])
xlabel('\omega/J_z');
ylabel('I(\omega)');
legend({'I_{aa}', 'I_{aa,ac}', '-I_{aa,cc}','I_{ac}','I_{ac,cc}','I_{cc}/3','I_{ab}','I_{bc}/2'}, 'Location', 'NorthEast');
hold off;
filename = ['H1_Raman_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100)];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);


hh=figure;%('Position',position);
hold on;
plot(Ev,I{1},Ev,-I{3}/3,Ev,I{6}/9,Ev,I{2},Ev,-I{5}/3,Ev,I{4},Ev,I{8}/2);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
%title(['Raman Spectrum for H-1 spinons: J_x=',num2str(Jx),', h=',num2str(h)])
xlabel('\omega/J_z');
ylabel('I(\omega)');
legend({'I_{aa}','-I_{aa,cc}/3','I_{cc}/9','I_{aa,ac}','-I_{ac,cc}/3','I_{ac}','I_{bc}/2'}, 'Location', 'NorthEast');
hold off;
filename = ['H1_Raman_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100),'_3'];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);


hh=figure;%('Position',position);
hold on;
plot(Ev,I{1},Ev,I{4},Ev,I{7});
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
%title(['Raman Spectrum for H-1 spinons: J_x=',num2str(Jx),', h=',num2str(h)])
xlabel('\omega/J_z');
ylabel('I(\omega)');
legend({'I_{aa}','I_{ac}','I_{ab}'}, 'Location', 'NorthEast');
hold off;
filename = ['H1_Raman_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100),'_2'];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

%Print the ratios
disp(mean(I{3}(I{1}>0)./I{1}(I{1}>0)))
disp(mean(I{6}(I{1}>0)./I{1}(I{1}>0)))
disp(mean(I{5}(I{2}>0)./I{2}(I{2}>0)))
disp(mean(I{8}(I{4}>0)./I{4}(I{4}>0)))



end

toc
end