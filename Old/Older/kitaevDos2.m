function [DE]= kitaevDos2(bins,t,a,lv,n)
%This version is probably a little slower, but can handle much larger choices of lv
%(probably as close to 10^8), whereas kitaevDos quits after 10^4


%input-parameter:
%bins...number of bins on the energy axis
%t...Energy of the elements on the diagonals next to the main diagonal
%(usually 2)
%a...distance to nearest neighbour (1 is fine)
%lv...the brillioun zone will be cut out of a rectangle shaped grid with
%     lv^2 points
%n...n=number of energy values / lv (1 is good))

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
kx=linspace(0,2*pi/(sqrt(3)*a),lv);
ky=linspace(0,2*pi/(3*a),lv);

num = lv*(lv+1)/2;
%I think it is most efficient here to claim the memory right away rather
%than building an array as we go.
EE=zeros(1,lv);
DE=zeros(1,bins);

Evalues=linspace(0,6,bins);

for i=1:lv
  for j=1:lv  
%calulation of the energy values and dismissal of energies calulated from k
%outside the brillioun zone 
  if j<=i
    %Energy Matrix of Band1 (E1) and Band 2 (E2)
    EE(j)=t*sqrt(1+4*cos(sqrt(3)*kx(i)*a/2)*cos(ky(j)*a/2)+4*cos(ky(j)*a/2)^2);
    %Energy Vektors with all possible Energy states
    %Elevband=[Elevband,EE];
  end
  end
  
  %Add to the histogram separately on each i value to save memory
  DE=DE+hist(EE(1:i),Evalues);
  
end
%Elevband = reshape(EE,1,[]);

%Normalize to lv^2/2 values over 12*lv*n bins 12 for the
%equivalent slices of the BZ and /2 for 2 sites per cell?
DE=DE/(6*num/(bins));

%[transpose(DE)]

hold on;
plot(Evalues,DE,'b');
title(['Density of states for Kitaev spinons: bins = ',num2str(bins),', t = ', num2str(t), ', a = ', num2str(a),', lv= ',num2str(lv),', n= ', num2str(n)])
xlabel('E');
ylabel('D(E)');

hold off;