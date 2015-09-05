
M=10;


for jind = (1:M)
    Jx = jind/M;
    
    temp = runit_H1(Jx);
    En(jind+1) = temp{1};
    st(jind+1) = temp{2};
    Evv{jind+1} = temp{8};
    DOSs{jind+1} = temp{7};
    
    temp = runit_H1_pi(Jx);
    En_pi(jind+1) = temp{1};
    st_pi(jind+1) = temp{2};
    Evv_pi{jind+1} = temp{8};
    DOSs_pi{jind+1} = temp{7};
end
Jxs = (0:M)/M;

hh=figure;%('Position',position); 
% hold on;
% a = errorbar(rand(4),rand(4));
% hold on
% b = errorbar(rand(4),rand(4));
% legend([a; b], {'a', 'b', 'c', 'd', 'a2', 'b2', 'c2', 'd2'});

%plot(Jxs(Jxs>0),En(Jxs>0));
% plot([Jxs(Jxs>0),Jxs(Jxs>0)], [En(Jxs>0),En_pi(Jxs>0)])
% errorbar([Jxs(Jxs>0),Jxs(Jxs>0)], [En(Jxs>0),En_pi(Jxs>0)],[st(Jxs>0),st_pi(Jxs>0)])
 aa = errorbar(Jxs(Jxs>0),En(Jxs>0),st(Jxs>0));
 hold on;
 %plot(Jxs(Jxs>0),En_pi(Jxs>0));
 bb = errorbar(Jxs(Jxs>0),En_pi(Jxs>0),st_pi(Jxs>0));
 
 set(aa                          , ...
  'Color'           , [1 0 0]    );
 set(bb                          , ...
  'Color'           , [0 0 1]    );
 
 hold off
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title('GS energy for H1')
xlabel('J_x/J_z');
ylabel('E_0/J');
legend([aa bb],{'E_0^0','E_0^\pi'}, 'Location', 'NorthEast');
filename = 'GS_energy_H1';
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

