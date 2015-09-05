function [Evalues,DEp,DEm,DE]= DOS3D4(bins,N)

%input-parameter:
%bins...number of bins on the energy axis
%emax...the max energy for the energy axis
%lv...the brillioun zone will be cut out of a rectangle shaped grid with
%     ~lv^3 points

%output-values:
%Evalues...energy values from the lowest to the highest energy occuring
%DE...how many times every energy value occurs


%additional output:
%plot of the densitiy of states

emax = 2.3; %the max value is sqrt(5) ~ 2.236

%initialization
Evalues=linspace(0,emax,bins);
DEp=zeros(1,bins); DEm=DEp; DE=DEp;
ep=zeros(2*N-1,2*N-1);
em=ep; e2=ep;

pts = (1-N):(N-1); pts=pts/N;
%treat x and y values as corresponding to column
%ans row values respectively
x=repmat(pts,2*N-1,1)*pi/2;
z=repmat(pts.',1,2*N-1)*pi/(3*sqrt(2));

d=abs(x)/sqrt(2)+abs(z)/3;

for y = ((1-N):(N-1))*29*pi/(36*sqrt(2)*N)
    %If this value is not positive I negate the energies out of the hist
    sn = sign( 29*pi/(36*sqrt(2)) - abs(y) - d );
    
%    u = cos(a + y/sqrt(2)); 
%    v = cos(b + 3*y/sqrt(2)); 
%    w = cos(c - 3 *y/sqrt(2)); 
% I only store values for this y until I histogram
    a = 4+2*( cos(3*x/sqrt(2) + y/sqrt(2) - z) + cos(x/sqrt(2)+3*y/sqrt(2)+z) );
    b = cos(sqrt(2)*(x+y)) + cos(x/sqrt(2)-y/sqrt(2)-z);
    b = 1 + 4*b.*(b+ cos(3*x/sqrt(2)-3*y/sqrt(2)+3*z));
    c = (a+2).^2 - 4*b; c=sqrt(c);
    ep= sqrt((a+2+c)/8).*sn;
    em= sqrt((a+2-c)/8).*sn;
    e2 = ep + em;
% Now histogram the values I found on this y slice      
    DEp=DEp+hist(ep(:),Evalues);
    DEm=DEm+hist(em(:),Evalues);
    DE =DE +hist(e2(:),Evalues);
end

A = 18/29*(2*N-1)^3*emax;
%Z = (pi^3)/24; 
DEp=DEp*bins/A;
DEm=DEm*bins/A;
DE =DE*bins/A;
DEp(1)=0; DEm(1)=0; DE(1)=0;

%[transpose(DE)]

clf;

hold on;
plot(Evalues,DEp,Evalues,DEm,Evalues,DE);
title(['DOS for HyperHoneycomb Kitaev spinons: bins = ',num2str(bins),', N= ',num2str(N)])
xlabel('E');
ylabel('D(E)');
hold off;