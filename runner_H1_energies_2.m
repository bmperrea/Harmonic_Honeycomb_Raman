
M=4;


for jind = (1:(M-1))
    N = jind*10;
    
    temp = runit_H1(1,N);
    EnN(jind) = temp{1};
    stN(jind) = temp{2};
    EvvN{jind} = temp{8};
    DOSsN{jind} = temp{7};
    
    temp = runit_H1_pi(1,N);
    En_piN(jind) = temp{1};
    st_piN(jind) = temp{2};
    Evv_piN{jind} = temp{8};
    DOSs_piN{jind} = temp{7};
end
Jxs = (1:M)*10;

EnN(4) = IH1_1{1};
stN(4) = IH1_1{2};
 EvvN{4} = IH1_1{8};
    DOSsN{4} = IH1_1{7};
    
En_piN(4) = IH1_pi_1{1};
st_piN(4) = IH1_pi_1{2};
 Evv_piN{4} = IH1_pi_1{8};
    DOSs_piN{4} = IH1_pi_1{7};    

hh=figure;%('Position',position); 
% hold on;
% a = errorbar(rand(4),rand(4));
% hold on
% b = errorbar(rand(4),rand(4));
% legend([a; b], {'a', 'b', 'c', 'd', 'a2', 'b2', 'c2', 'd2'});

%plot(Jxs(Jxs>0),En(Jxs>0));
% plot([Jxs(Jxs>0),Jxs(Jxs>0)], [En(Jxs>0),En_pi(Jxs>0)])
% errorbar([Jxs(Jxs>0),Jxs(Jxs>0)], [En(Jxs>0),En_pi(Jxs>0)],[st(Jxs>0),st_pi(Jxs>0)])
 aa = errorbar(Jxs(Jxs>0),EnN(Jxs>0),stN(Jxs>0));
 hold on;
 %plot(Jxs(Jxs>0),En_pi(Jxs>0));
 bb = errorbar(Jxs(Jxs>0),En_piN(Jxs>0),st_piN(Jxs>0));
 
 set(aa                          , ...
  'Color'           , [1 0 0]    );
 set(bb                          , ...
  'Color'           , [0 0 1]    );
 
 hold off
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title('GS energy for H1')
xlabel('N');
ylabel('E_0/J');
legend([aa bb],{'E_0^0','E_0^\pi'}, 'Location', 'NorthEast');
filename = 'GS_energy_H1_N';
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

