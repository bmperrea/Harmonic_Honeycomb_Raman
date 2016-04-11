function [DE,Evalues]= graphene5(lv,bins)
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

%x values along the row
kx=repmat(linspace(0,2*pi/sqrt(3),lv), lv,1);
%y values along the column index
ky=repmat(linspace(0,2*pi/3,lv), lv,1).';

Es=abs( 1 + exp(1i*(sqrt(3)*kx + 3*ky)/2) + exp(1i*(-sqrt(3)*kx + 3*ky)/2) );

%A logical determining which points are in the Brillouin Zone within the
%bounding rectangle
% The BZ is a diamond, of which we only consider the upper right corner,
% while the other ones are the same by symmetry.
%The unit vectors are (+sqrt(3),3) and (-sqrt(3),3)
boole = ( sqrt(3)*kx + 3*ky < 2*pi );

%Take only the ones in the BZ
Esvec = Es(boole);

%The energy values to histogram over
emin = min(min(Esvec));
emax = max(max(Esvec));
Evalues=linspace(emin,emax,bins);

%counting Energy values
DE=hist(Esvec(:),Evalues);
DE = DE*bins/( (emax-emin)*sum(DE) );

figure;
hold on;
plot(Evalues,DE,'b');
title(['Density of states for square lattices: lv= ',num2str(lv),', bins = ', num2str(bins)])
xlabel('E');
ylabel('\rho(E)');
hold off;