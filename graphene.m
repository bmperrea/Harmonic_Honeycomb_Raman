function [DE,Evalues]= graphendos(eps,t,a,lv,n)
%input-parameter:
%eps...Energy of the main diagonal elements
%t...Energy of the elements on the diagonals next to the main diagonal
%a...distance to nearest neighbour
%lv...the brillioun zone will be cut out of a rectangle shaped grid with
%     lv^2 points
%n...n=number of energy values / lv

%output-values:
%Evalues...energy values from the lowest to the highest energy occuring
%DE...how many times every energy value occurs


%additional output:
%plot of the densitiy of states



%kx...vector with lv values between the highest and the lows k-value in
%     x-direction
%ky...vector with lv values between the highest and the lows k-value in
%     y-direction
%kx and ky form a rectangle shaped grid around the brillioun zone
kx=linspace(0,2*pi./(sqrt(3)*a),lv);
ky=linspace(0,4/3*pi./a,lv);

E=zeros(lv,lv);


Elevband1=[];
Elevband2=[];


for i=1:lv
  for j=1:lv  
%calulation of the energy values and dismissal of energies calulated from k
%outside the brillioun zone 
  if j<=i
    %Energy Matrix of Band1 (E1) and Band 2 (E2)
    E1(i,j)=eps-t*sqrt(1+4*cos(sqrt(3)*kx(i)*a/2)*cos(ky(j)*a/2)+4*cos(ky(j)*a/2)^2);
    E2(i,j)=eps+t*sqrt(1+4*cos(sqrt(3)*kx(i)*a/2)*cos(ky(j)*a/2)+4*cos(ky(j)*a/2)^2);
    %Energy Vektors with all possible Energy states
    Elevband1=[Elevband1,E1(i,j)];
    Elevband2=[Elevband2,E2(i,j)];
  end
  end
  
end


Evalues=linspace(min(Elevband1),max(Elevband2),lv*n);
%counting Energy values
[DE]=hist(Elevband1,Evalues)+hist(Elevband2,Evalues);

[transpose(Evalues),transpose(DE)]

hold on;
plot(Evalues,DE,'b');
title(['Density of states for square lattices: eps = ',num2str(eps),', t = ', num2str(t), ', a = ', num2str(a),', lv= ',num2str(lv),', n= ', num2str(n)])
xlabel('E');
ylabel('D(E)');

hold off;