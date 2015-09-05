function [DE,Evalues]= graphene2(t,a,lv)
%input-parameter:
%t...Energy of the elements on the diagonals next to the main diagonal
%a...distance to nearest neighbour
%lv...the brillioun zone will be cut out of a rectangle shaped grid with
%     lv^2 points

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

%A row vector of x values
kx=linspace(0,2*pi./(sqrt(3)*a),lv);
%A column vector of y values
ky=linspace(0,4/3*pi./a,lv).';

%A row vector of ones for when we don't use x-values
onevec=ones(1,lv);

E=t*sqrt(1+4*cos(ky*a/2)*cos(sqrt(3)*kx*a/2)+4*(cos(ky*a/2).^2)*onevec );

Elevband=reshape(E,1,lv^2);


Evalues=linspace(min(Elevband),max(Elevband),lv);
%counting Energy values
[DE]=hist(Elevband,Evalues)

[transpose(Evalues),transpose(DE)]

hold on;
plot(Evalues,DE,'b');
title(['Density of states for square lattices: ', ',t = ', num2str(t), ', a = ', num2str(a),', lv= ',num2str(lv)])
xlabel('E');
ylabel('D(E)');

hold off;