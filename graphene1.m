function [DE,Evalues]= graphene1(lv,bins)
%Based on a code by Peter Hadley
%http://lampx.tugraz.at/~hadley/ss1/bands/tbtable/2d_graphene_dos.html?

%input-parameter:
%lv...the brillioun zone will be cut out of a rectangle shaped grid with
%     lv^2 points. Try 500.
%bins.. the number of energy values to histogram the energies over. Try 200

%output-values:
%Evalues...energy values from the lowest to the highest energy occuring
           %This is on a scale of -3 to 3 in units of t, the hopping
           %integral
%DE... a histogram of energy values from the k-points that are sampled

%additional output:
%plot of the densitiy of states


%kx...vector with lv values between the highest and the lows k-value in
%     x-direction
%ky...vector with lv values between the highest and the lows k-value in
%     y-direction
%kx and ky form a rectangle shaped grid around the brillioun zone
kx=linspace(0,2*pi/sqrt(3),lv);
ky=linspace(0,2*pi/3,lv);

Es=zeros(lv,lv);
Elevband=[];

for l=1:lv
  for j=1:lv  
%calulation of the energy values and dismissal of energies calulated from k
%outside the brillioun zone 
  if j<=l
    %Energy Matrix of Band1 (E1) and Band 2 (E2)
    E1=abs( 1 + exp(1i*(sqrt(3)*kx(l) + 3*ky(j))/2) + exp(1i*(-sqrt(3)*kx(l) + 3*ky(j))/2) );
    Elevband = [Elevband,E1];
  end
  end
  
end

emin = min(Elevband);
emax = max(Elevband);
Evalues=linspace(emin,emax,bins);
%counting Energy values
[DE]=hist(Elevband,Evalues);
DE = DE*bins/( (emax-emin)*sum(DE) );

figure;
hold on;
plot(Evalues,DE,'b');
title(['Density of states for the honeycomb lattice: lv= ',num2str(lv),', bins = ', num2str(bins)])
xlabel('E');
ylabel('\rho(E)');
hold off;