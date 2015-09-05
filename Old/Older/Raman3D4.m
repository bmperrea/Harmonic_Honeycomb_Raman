function [Evp,DEp,Ipp,Evm,DEm,Imm,Ev,DE,Ipm]= Raman3D4(bins,N)

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

%The Hamiltonian used comes from EK Lee et al. arXiv:1308.6592v2


%input data
Jx = 1; 
Jy = Jx; 
Jz = 1;

hx = 1;
hy = 1;
hz = 1;

%The following can be used to ensure that H is diagonalized by the unitary
%matrix below by checking that the Raman spectra are zero.
hx=Jx;hy=Jy;hz=Jz;

epmin = 0.507722357;
epmax = 3/2; erp = epmax-epmin;
emmin = 0; 
emmax = sqrt(5)/2; erm = emmax-emmin;
emin  = 1;
emax  = sqrt(5); er = emax-emin;

%initialization
%Evalues=linspace(0,emax,bins);

DEp=zeros(bins,1);
DEm=zeros(bins,1);
DE=zeros(bins,1);

Ipp=zeros(bins,1);
Imm=zeros(bins,1);
Ipm=zeros(bins,1);

L=2*N-1;

%z=zeros(L,L); ep=z; em=ep;
r(1)=-1;r(2)=-1;r(3)=1;r(4)=1;
%q=(1:(L)^2)'; 

Evp = epmin + (1:bins)'*erp/bins;
Evm = emmin + (1:bins)'*erm/bins;
Ev  = emin  + (1:bins)'*er/bins;

pts = (1-N):(N-1); pts=pts/N;
%treat x and y values as corresponding to column
%ans row values respectively
x=repmat(pts,L,1)*pi/2;
z=repmat(pts.',1,L)*pi/(3*sqrt(2));

U=zeros(L,L,4,4); V=U; Udag = U; R=U;
en = zeros(L,L,4);

Z=0;

for y = ((1-N):(N-1))*29*pi/(36*sqrt(2)*N)
    %If this value is not positive I negate the energies out of the hist
    sn = sign( 29*pi/(36*sqrt(2)) - abs(y) - abs(x)/sqrt(2)+abs(z)/3 );
  %  [L L2]=size(sn(:)>=0);
    
    %Now I do arithmetic as if x,y, and z were numbers but with x and a
    %matrices and doing only elementwise operations with them.
    
    %I only store values for a given y and then reduce the information into
    %a histogram before looping back to keep from overloading the memory
    
    %My computer handles about 2 to 5 * 10^8 doubles at once in memory so I
    %expect this code to max out the memory around N~3000
    
    
    %The Eigenvalues
 %   ep = 6 + 2*( cos(3*x/sqrt(2) + y/sqrt(2) - z) + cos(x/sqrt(2)+3*y/sqrt(2)+z) );
 %   em = sqrt(ep.^2 - 2*(1+4*( cos(sqrt(2)*(x+y)) + cos(x/sqrt(2)-y/sqrt(2)-z) ) ));
 %   ep = ep + em;

    
    %The Eigenvectors are not real and require some work (as well as
    %complex numbers). We define a number of intermediate variables here
    %This comes from the paper cited above and the Mathematica
    %notebook entitled 3D Kitaev
    A = exp(1i*(4*x - 4*y + 2*sqrt(2)*z)/sqrt(2)) + exp(1i*(x - 5*y + 3*sqrt(2)*z)/sqrt(2)); 
    B = 1 + exp(1i*(-x - 3*y - sqrt(2)* z)/sqrt(2));
    a = 4+2*( cos(3*x/sqrt(2) + y/sqrt(2) - z) + cos(x/sqrt(2)+3*y/sqrt(2)+z) );
    b = cos(sqrt(2)*(x+y)) + cos(x/sqrt(2)-y/sqrt(2)-z);
    b = 1 + 4*b.*(b+ cos(3*x/sqrt(2)-3*y/sqrt(2)+3*z));
    c = (a+2).^2 - 4*b; c=sqrt(c);
    ep= sqrt((a+2+c)/8);
    em= sqrt((a+2-c)/8);

    %The unnormalized eigenvectors
    %I do the matrix product on 4x4  
    en(:,:,1)=-em; en(:,:,2)=em; en(:,:,3)=-ep; en(:,:,4)=ep;
    for k=1:4
        U(:,:,k,1)= abs(A-B).^2 + r(k)*c.*(1-2*en(k));
        U(:,:,k,2)= 2i*imag(B.*conj(A)) + 2*en(k).*(abs(A).^2 - abs(B).^2);
        U(:,:,k,3)= conj(A.*B).*(A-B) + conj(A-B).*(1+2*en(k)).^2;
        U(:,:,k,4)= -r(k)*conj(A+B).*c/2 - (abs(A).^2 - abs(B).^2).*conj(A-B)/2;
        R(:,:,k,k)= -sign(en(k))*hz;
    end
    %Normalize it
    nm = sqrt( sum( abs(U).^2 ,4) ); nm = repmat(nm,1,1,1,4); U=U./nm;

    %Now create the rest of the Hamiltonian for the different A and B
    A = hx*exp(1i*(4*x - 4*y + 2*sqrt(2)*z)/sqrt(2)) + hy*exp(1i*(x - 5*y + 3*sqrt(2)*z)/sqrt(2)); 
    B = hx + hy*exp(1i*(-x - 3*y - sqrt(2)* z)/sqrt(2));
    R(:,:,1,3)=A-B;R(:,:,1,4)=-A-B;R(:,:,2,3)=A+B;R(:,:,2,4)=-A+B;

    for n=1:4
        for m=1:4
            R(:,:,n,m) = (R(:,:,n,m)+conj(R(:,:,m,n)))/4;
            Udag(:,:,n,m)= conj( U(:,:,m,n) );
        end
    end
    
    for n=[1 3]
        for m=1:4
            for k=1:4
                V(:,:,n,m)=Udag(:,:,n,k).*R(:,:,k,m);
            end
        end
    end
    
    for n=[1 3]
        for m=[2 4]
            for k=1:4
                V(:,:,n,m)=V(:,:,n,k).*U(:,:,k,m);
            end
        end
    end
    
    %The Delta weight in the Raman term (turned into column vectors)
    Wpp = 2*pi*abs(V(:,:,1,2)).^2;
    Wmm = 2*pi*abs(V(:,:,3,4)).^2;
    Wpm = 2*pi*(abs(V(:,:,1,4)).^2 + abs(V(:,:,3,2)).^2);
        
    [histw, vint] = histwc(ep(sn>=0),Wpp(sn>=0),epmin,epmax,bins);
    Ipp = Ipp + histw;
    [histw, vint] = histwc(em(sn>=0),Wmm(sn>=0),emmin,emmax,bins);
    Imm = Imm + histw;
    [histw, vint] = histwc(ep(sn>=0)+em(sn>=0),Wpm(sn>=0),emin,emax,bins);
    Ipm = Ipm + histw;
    
    DEp=DEp+hist(ep(:),Evp)';
    DEm=DEm+hist(em(:),Evm)';
    DE =DE +hist(em(:)+ep(:),Ev)';
    
    Z=Z+size(sn(sn(:)>=0),1);
    
end

%Z = 18/29*L^3;
%Z = (pi^3)/24; 
DEp=DEp*bins/(Z*erp);
DEm=DEm*bins/(Z*erm);
DE =DE *bins/(Z*er );

Ipp=Ipp*bins/(Z*erp);
Ipm=Ipm*bins/(Z*er );
Imm=Imm*bins/(Z*erm);

%DEp(1)=0; DEm(1)=0;

%Plot the DOS
clf;
hold on;
plot(Evp,DEp,Evm,DEm,Ev,DE);
title(['DOS for HyperHoneycomb Kitaev spinons: bins = ',num2str(bins),', N= ',num2str(N)])
xlabel('E');
ylabel('D(E)');
hold off;