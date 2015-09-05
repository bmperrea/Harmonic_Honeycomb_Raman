%This a script that calls the function that gives Raman spectra for the
%Kitaev Hyperhoneycomb (3D) model
function I = runit14_8_f(N,Jxx,h)

tic

disp(['Jxx= ',num2str(Jxx),' h= ',num2str(h)])

% 
%N=200;
nn=3;
bins = 200;

%initialization
Dt = cell(1,nn); Iv = Dt; Ier = Dt;
Ivt = cell(8,nn); It = Ivt;
I = cell(1,11); Iermax = I; Iermean = I;

%Runs for different parameters
Jz = 1;
flag = 1;


for Jx = Jxx

  %  if Jx == 1
  %  end
        
for n = 1:nn
    %The function gets run nn times so we can do statistics
    out = KitaevRaman14_8_ac(N,bins,flag,nn,Jx,Jz,h);
    Ev = out{1}; %Same every time
    Dt{n}  = out{2}; %D = [Dpp,Dmm,Dpm]
    for m=3:12
        Ivt{m,n} = out{m+1}; % I{m} = [Ipp{m},Imm{m},Ipm{m}];
        It{m,n}  = Ivt{m,n}(1) + Ivt{m,n}(2) + Ivt{m,n}(3);
    end
end

for m=3:12
    Ivtt = Ivt(m,:);    
    Iv{m}  = mean(Ivtt,1,'c');
    Ier{m} = std(Ivtt,1,'c');
    %Ier2{m}= max(Ivtt,1,'c');
    Iermean{m} = mean(Ier{m});
    Iermax{m} = max(Ier{m});
    disp(Iermean{m})
    disp(Iermax{m})
    
    I{m} = sum(Iv{m},2);
    %I{m} = I{m}.*sign(I{m});
end

D = mean(Dt,1,'c');
Der = std(Dt,1,'c');
Der2= max(Der);
Der3= mean(Der);
disp(Der3)
disp(Der2)


DD = D(:,1)+D(:,2)+2*D(:,3);
dE = Ev(bins)/bins;
for m = 1:nn
    EEE(m)=sum( ( Dt{m}(:,1)+Dt{m}(:,2) ).*(Ev/4)*dE/4 )/2;
end
 format long e
disp(mean(EEE));
    format short
disp(std(EEE));

I{13} = DD;
I{2} = D(:,1)+D(:,2);
I{1} = Ev;
%I{11} = D(:,1)+D(:,2); %[mean(EEE),std(EEE)];

pixs = get(0,'screensize');
% 
% %Plot DOS
% %position = [pixs(3)/2   20       pixs(3)/2-114   pixs(4)/2-50];
% hh=figure;%('Position',position); hold on;
% plot(Ev,DD);
% %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
% %title(['DOS for HyperHoneycomb Kitaev spinons: J_x=',num2str(Jx),', h=',num2str(h)])
% xlabel('\omega/J^z');
% ylabel('DOS');
% %set(gca,'XTick',-3:3); 
% %set(gca,'YTick',2*(0:5));
% %hold off;
% filename = ['3D_DOS_2_Jx_',num2str(round(100*Jx)),'_h_',num2str(100*h),'_28'];
% saveas(hh,filename)
% print(hh, '-dpng', filename);
% print(hh, '-depsc', filename);
% 
% position = [pixs(3)/2   20       pixs(3)/2-114   pixs(4)/2-50];
% hh=figure('Position',position); hold on;
% plot(Ev,D(:,2),Ev,D(:,1),Ev,2*D(:,3),Ev,DD);
% %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
% %title(['DOS for HyperHoneycomb Kitaev spinons: J_x=',num2str(Jx),', h=',num2str(h)])
% xlabel('\omega/J^z');
% ylabel('DOS');
% legend({'\rho_{--}', '\rho_{++}', '\rho_{+-}','\rho_{total}'}, 'Location', 'NorthEast');
% %set(gca,'XTick',-3:3); 
% %set(gca,'YTick',2*(0:5));
% hold off;
% filename = ['3D_DOS_2_Jx_',num2str(round(100*Jx)),'_h_',num2str(100*h),'_8'];
% saveas(hh,filename)
% print(hh, '-dpng', filename);
% print(hh, '-depsc', filename);
% % 
% % if h==0 && Jxx ==1
% %     ptsE = [1.65685, 4, 5.65685, 8, 8.94427, 9.65685, 12];
% %     
% %     hh=figure;%('Position',position); 
% %     hold on;
% %     plot(Ev,D(:,2),Ev,D(:,1),Ev,2*D(:,3),Ev,DD);
% %     %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
% %     %title(['DOS for HyperHoneycomb Kitaev spinons: J_x=',num2str(Jx),', h=',num2str(h)])
% %     xlabel('\omega/J_z');
% %     ylabel('DOS');
% %     legend({'\rho_{--}', '\rho_{++}', '\rho_{+-}','\rho_{total}'}, 'Location', 'NorthEast');
% %     for j=1:7
% %         plot([ptsE(j),ptsE(j)],[0,max(DD)])
% %     end
% %     hold off;
% %     filename = ['3D_DOS_2_Jx_',num2str(round(100*Jx)),'_h_',num2str(100*h),'_5'];
% %     saveas(hh,filename)
% %     print(hh, '-dpng', filename);
% % end
% 
% %Plot the ramans
% % 
% hh=figure;%('Position',position);
% hold on;
% plot(Ev,I{1},Ev,I{2},Ev,-I{3},Ev,I{4},Ev,I{5},Ev,I{6}/3,Ev,I{7},Ev,I{8}/2);
% %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
% %title(['Raman Spectrum for HyperHoneycomb Kitaev spinons: J_x=',num2str(Jx),', h=',num2str(h)])
% xlabel('\omega/J^z');
% ylabel('I(\omega)');
% legend({'I_{aa}', 'I_{aa,ac}', '-I_{aa,cc}','I_{ac}','I_{ac,cc}','I_{cc}/3','I_{ab}','I_{bc}/2'}, 'Location', 'NorthEast');
% hold off;
% filename = ['3D_Raman_2_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100),'_8ac'];
% saveas(hh,filename)
% print(hh, '-dpng', filename);
% print(hh, '-depsc', filename);
% 
% 
% hh=figure;%('Position',position);
% hold on;
% plot(Ev,I{1},Ev,-I{3}/3,Ev,I{6}/9,Ev,I{2},Ev,-I{5}/3,Ev,I{4},Ev,I{8}/2);
% %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
% %title(['Raman Spectrum for HyperHoneycomb Kitaev spinons: J_x=',num2str(Jx),', h=',num2str(h)])
% xlabel('\omega/J^z');
% ylabel('I(\omega)');
% legend({'I_{aa}','-I_{aa,cc}/3','I_{cc}/9','I_{aa,ac}','-I_{ac,cc}/3','I_{ac}','I_{bc}/2'}, 'Location', 'NorthEast');
% hold off;
% filename = ['3D_Raman_2_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100),'_38ac'];
% saveas(hh,filename)
% print(hh, '-dpng', filename);
% print(hh, '-depsc', filename);
% 
% 
% hh=figure;%('Position',position);
% hold on;
% plot(Ev,I{1}.*sign(I{1}),Ev,I{4}.*sign(I{4}),Ev,I{7}.*sign(I{7}));
% %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
% %title(['Raman Spectrum for HyperHoneycomb Kitaev spinons: J_x=',num2str(Jx),', h=',num2str(h)])
% xlabel('\omega/J^z');
% ylabel('I(\omega)');
% legend({'I_{aa}','I_{ac}','I_{ab}'}, 'Location', 'NorthEast');
% hold off;
% filename = ['3D_Raman_2_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100),'_28'];
% saveas(hh,filename)
% print(hh, '-dpng', filename);
% print(hh, '-depsc', filename);
% 
% %Print the ratios
% disp(mean(I{3}(I{1}>0)./I{1}(I{1}>0)))
% disp(mean(I{6}(I{1}>0)./I{1}(I{1}>0)))
% disp(mean(I{5}(I{2}>0)./I{2}(I{2}>0)))
% disp(mean(I{8}(I{4}>0)./I{4}(I{4}>0)))
% 
% % hh=figure;%('Position',position);
% % hold on;
% % plot(log(Ev(Ev<5*h)),log(I{1}(Ev<5*h)),log(Ev(Ev<5*h)),log(I{5}(Ev<5*h)));
% % %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
% % %title(['Raman Spectrum for HyperHoneycomb Kitaev spinons: J_x=',num2str(Jx),', h=',num2str(h)])
% % xlabel('log(\omega/J_z)');
% % ylabel('log()');
% % legend({'I_{aa}', 'I_{aa,ac}', 'I_{aa,cc}','I_{ac}','I_{ac,cc}','I_{cc}'}, 'Location', 'NorthEast');
% % hold off;
% % filename = ['3D_Raman_2_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100),'_3'];
% % saveas(hh,filename)
% % print(hh, '-dpng', filename);
% 
% 
% 
% for m=1:6
%     
%   %  if max(abs(I{m})) > 2*max(Iermax{m})   
%   
%   
% % hh=figure;%('Position',position);
% % hold on;
% % plot(Ev,Iv{m}(1),Ev,Iv{m}(2),Ev,Iv{m}(3),Ev,I{m});
% % %errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
% % %title(['Raman Spectrum for HyperHoneycomb spinons: J_x=',num2str(Jx),', h=',num2str(h)])
% % xlabel('\omega/J_z');
% % ylabel('I(\omega)');
% % legend({'I_{--}','I_{+-}','I_{++}','I_{total}'}, 'Location', 'NorthEast');
% % hold off;
% % filename = ['3D_Raman_2_Jx_',num2str(round(100*Jx)),'_h_',num2str(h*100),...
% %     'pol_',num2str(m1{m}),num2str(m2{m})];
% % saveas(hh,filename)
% % print(hh, '-dpng', filename);
% % %hold off;
% % %clf;
% 
%  %   else
%         
% %        disp([max(abs(I{m})),max(Iermax{m})])
% 
% %    end
% 
% end
% 
% flag=0;
% 
 end

%End New Code

toc
end