function out = KitaevRaman10(N,bins,type,flag,nn,Jx)
%This function was written by Brent Perreault
%The basic Hamiltonian can be compared with EK Lee et al. arXiv:1308.6592v2


%output form: out = [Ev,DEp,DEm,DE,Ipp,Imm,Ipm]

% Ev   Energy bins for Raman Intensity. The appropriate bins for DOS are Ev/2
% DEp  Density of states for the upper band
% DEm  DOS for lower band
% DE   Total Density of states
% Ipp  Raman Intensity for two spinons in the upper band
% Imm  Raman intensity for two spinons in the lower band
% Ipm  Raman intensity for two spinons, one from each band

%additional output:
%Time taken so far, and estimated end time
%Plot of the densitiy of states (commented out for now)

%input-parameters:

% 2*N    the number of points in a dimension (2*N ~ 100 takes minutes)
% flag   a boolean telling whether to estimate time on this run
% nn     the number of times that you plan to do a similar calculation (for error estimation) 
% type   a string specifying which Raman operator
% Jx     the value of J_x/J_z , where J_y = J_x
% bins   number of bins on the energy axis (~2*N is a good choice)


%Initial data

%the max energy for the energy axis
emax = 3*max(Jx,1) + .2;

%The values of the Kitaev couplings
Jy = Jx;  %Jy is always same as Jz
Jz = 1;  %Jz is always 1 the way I coded it

%pick the right Raman operator
switch type
    case 'xx'
        hx = 2/3;          hy = 2/3;           hz = 0;
    case 'yy'
        hx = 1/3;          hy = 1/3;           hz = 0; 
    case 'zz'
        hx = 1;          hy = 1;           hz = 2; 
    case 'xy'
        hx = -1/sqrt(9/2);  hy = -1/sqrt(9/2);   hz = 0;    
    case 'xz'
        hx = -1/sqrt(3/2);   hy = 1/sqrt(3/2);     hz = 0;  
    case 'yz'
        hx = -1/sqrt(3);   hy = 1/sqrt(3);   hz = 0;   
    case '++'
        hx = 1; hy = 1;  hz = 0;
    case '--'
        hx = 1; hy = -1; hz = 0;
    case 'zzz'
        hx = 0; hy = 0;  hz =1;
    case '+-'
        hx = 1; hy = 1;  hz = 0;
    case '+z'
        hx = 1; hy = 1;  hz = 0;
    case '-z'
        hx = 1; hy = -1;  hz = 0;
end
hx = hx*Jx; hy = hy*Jy; hz = hz*Jz;

%The following can be used to ensure that H is diagonalized by the unitary
%matrix below by checking that the Raman spectra are zero.
%hx=Jx;hy=Jy;hz=Jz;

%Now pick the correct second Raman operator if differen from the first
hx1 = hx; hy1 = hy; hz1 = hz;
if strcmp(type,'xx-xz')
    hx2 = -1/sqrt(6)*Jx;   hy2 = 1/sqrt(6)*Jy;     hz2 = 0; 
    %factor2 = 1;
elseif strcmp(type,'+-') || strcmp(type,'-+')
    hx2 = 1; hy2 = -1;  hz2 = 0;
    %factor2 = 1;
else
    hx2 = hx1;      hy2 = hy1;       hz2=hz1;
    %factor2 = 2;
end


%A convenient notation
L=2*N;

%Initialization of Arrays

DEp=zeros(bins,1);
DEm=zeros(bins,1);
DE=zeros(bins,1);

Ipp=zeros(bins,1);
Imm=zeros(bins,1);
Ipm=zeros(bins,1);

Ev = (1:bins)'*emax/bins;

r1 = rand; r2 = rand; r3 = rand;

pts = (-N):(N-1);
%We treat x and y values as corresponding to column and row values respectively
x=repmat((pts+r1)/N,L,1)*pi/2;
z=repmat((pts.' + r3)/N,1,L)*pi/(3*sqrt(2));

U=zeros(L,L,4,4); 
R=U;
norm=zeros(L,L,4);

%Z counts the number of k-points falling withing the actual BZ (volume of BZ)
Z=0;

    %Now I do arithmetic as if x,y, and z were numbers but with x and a
    %matrices and doing only elementwise operations with them.
    
    %I only store values for a given y and then reduce the information into
    %a histogram before looping back to keep from overloading the memory
    
    %My computer handles about 2 to 5 * 10^8 doubles at once in memory so I
    %expect this code to max out the memory around N~3000
for ind = pts
    y = (ind + r2)*29*pi/(36*sqrt(2)*N);
    
    
    %I display an estimate of the end time
    aleph = 1 + round( (200/N)^2 ); %iterations for a few seconds of computation
    %Or an iteration, which ever is longer
    if ind == -N+2 && flag   %the first few may be slower due to initialization
        cl1 = clock; 
    end
    if ind == -N+2+aleph && flag
        cl2 = clock;
        time = (cl2-cl1)*(2*N)/aleph *9*nn;
        %cl = cl1 + time*(2*N-2)/(2*N);
        
        format shortg
        disp('approximate time to take:')
        disp( datestr(time(6)/24/3600, 'DD-HH:MM:SS') )
        format
    end
        
        
    %Only values above the line below are in the BZ (which is part of the
    %rectangle I have created in x,y,z space.
    %If this value is not positive I negate the energies out of the hist
    sn = sign( 29*pi/(36*sqrt(2)) - abs(y) - abs(x)/sqrt(2)+abs(z)/3 );
     


    %the following code was converted from Mathematica with ToMatlab
    
    %Some preliminary values
    A = exp(1).^(sqrt(-1).*2.^(-1/2).*(4.*x+(-4).*y+2.*2.^(1/2).*z))...
        + exp(1).^(sqrt(-1).*2.^(-1/2).*(x+(-5).*y+3.*2.^(1/2).*z));
    B = 1+exp(1).^(sqrt(-1).*2.^(-1/2).*((-1).*x+(-3).*y+(-1).*2.^(1/2).*z));
    A = A*Jx; B=B*Jx;
    
    %The Hamiltonian
    H(:,:,1,1) = (Jz/2); H(:,:,1,2) = 0; H(:,:,1,3) = (1/4).*(A+(-1).*B); H(:,:,1,4) = (1/4).*((-1).*A+(-1).*B);
    H(:,:,2,1) = 0; H(:,:,2,2) = (-Jz/2); H(:,:,2,3) = (1/4).*(A+B); H(:,:,2,4) = (1/4).*((-1).*A+B);
    H(:,:,3,1) = (1/4).*(conj(A)+(-1).*conj(B)); H(:,:,3,2) = (1/4).*(conj(A)+conj(B)); H(:,:,3,3) = (Jz/2); H(:,:,3,4) = 0;
    H(:,:,4,1) = (1/4).*((-1).*conj(A)+(-1).*conj(B)); H(:,:,4,2) = (1/4).*((-1).*conj(A)+conj(B)); H(:,:,4,3) = 0; H(:,:,4,4) = (-Jz/2);
    
    %The Raman Operators R and R2
    A = A*hx/Jx; B = B*hx/Jx;
    
    R(:,:,1,1) = (hz/2); R(:,:,1,2) = 0; R(:,:,1,3) = (1/4).*(A+(-1).*B); R(:,:,1,4) = (1/4).*((-1).*A+(-1).*B);
    R(:,:,2,1) = 0; R(:,:,2,2) = (-hz/2); R(:,:,2,3) = (1/4).*(A+B); R(:,:,2,4) = (1/4).*((-1).*A+B);
    R(:,:,3,1) = (1/4).*(conj(A)+(-1).*conj(B)); R(:,:,3,2) = (1/4).*(conj(A)+conj(B)); R(:,:,3,3) = (hz/2); R(:,:,3,4) = 0;
    R(:,:,4,1) = (1/4).*((-1).*conj(A)+(-1).*conj(B)); R(:,:,4,2) = (1/4).*((-1).*conj(A)+conj(B)); R(:,:,4,3) = 0; R(:,:,4,4) = (-hz/2);
    
    hx = hx2;
    hz = hz2;
    A = A*hx/Jx; B = B*hx/Jx;
    
    R2(:,:,1,1) = (hz/2); R2(:,:,1,2) = 0; R2(:,:,1,3) = (1/4).*(A+(-1).*B); R2(:,:,1,4) = (1/4).*((-1).*A+(-1).*B);
    R2(:,:,2,1) = 0; R2(:,:,2,2) = (-hz/2); R2(:,:,2,3) = (1/4).*(A+B); R2(:,:,2,4) = (1/4).*((-1).*A+B);
    R2(:,:,3,1) = (1/4).*(conj(A)+(-1).*conj(B)); R2(:,:,3,2) = (1/4).*(conj(A)+conj(B)); R2(:,:,3,3) = (hz/2); R2(:,:,3,4) = 0;
    R2(:,:,4,1) = (1/4).*((-1).*conj(A)+(-1).*conj(B)); R2(:,:,4,2) = (1/4).*((-1).*conj(A)+conj(B)); R2(:,:,4,3) = 0; R2(:,:,4,4) = (-hz/2);
    
    ep = (1/2).*2.^(-1/2).*(2+a+c.^(1/2)).^(1/2);
    em = (1/2).*2.^(-1/2).*(2+a+(-1).*c.^(1/2)).^(1/2);
    
    Hc = 
    
    
    
%The unitary transformations
    
U(:,:,1,1) = c.^(1/2).*((-1/2)+em)+(1/2).*(A+(-1).*B).*(conj(A)+(-1).*conj(B));
U(:,:,1,2) = (1/2).*(B.*conj(A)+(-1).*A.*conj(B))+(-1).*em.*(A.*conj(A)+(-1).*B.*conj(B));
U(:,:,1,3) = (1/2).*(1+(-2).*em).^2.*(conj(A)+(-1).*conj(B))+(1/2).*(A+(-1).*B).*conj(A).*conj(B);
U(:,:,1,4) = (1/4).*dm;

U(:,:,2,1) = c.^(1/2).*((-1/2)+(-1).*em)+(1/2).*(A+(-1).*B).*(conj(A)+(-1).*conj(B));
U(:,:,2,2) = (1/2).*(B.*conj(A)+(-1).*A.*conj(B))+em.*(A.*conj(A)+(-1).*B.*conj(B));
U(:,:,2,3) = (1/2).*(1+2.*em).^2.*(conj(A)+(-1).*conj(B))+(1/2).*(A+(-1).*B).*conj(A).*conj(B);
U(:,:,2,4) = (1/4).*dm;

U(:,:,3,1) = c.^(1/2).*((1/2)+(-1).*ep)+(1/2).*(A+(-1).*B).*(conj(A)+(-1).*conj(B));
U(:,:,3,2) = (1/2).*(B.*conj(A)+(-1).*A.*conj(B))+(-1).*ep.*(A.*conj(A)+(-1).*B.*conj(B));
U(:,:,3,3) = (1/2).*(1+(-2).*ep).^2.*(conj(A)+(-1).*conj(B))+(1/2).*(A+(-1).*B).*conj(A).*conj(B);
U(:,:,3,4) = (-1/4).*dp;

U(:,:,4,1) = c.^(1/2).*((1/2)+ep)+(1/2).*(A+(-1).*B).*(conj(A)+(-1).*conj(B));
U(:,:,4,2) = (1/2).*(B.*conj(A)+(-1).*A.*conj(B))+ep.*(A.*conj(A)+(-1).*B.*conj(B));
U(:,:,4,3) = (1/2).*(1+2.*ep).^2.*(conj(A)+(-1).*conj(B))+(1/2).*(A+(-1).*B).*conj(A).*conj(B);
U(:,:,4,4) = (-1/4).*dp;
 

%The first Raman Operator

hx = hx1;
hy=  hy1;
hz = hz1;

    A = hx*exp(1).^(sqrt(-1).*2.^(-1/2).*(4.*x+(-4).*y+2.*2.^(1/2).*z))...
        + hy*exp(1).^(sqrt(-1).*2.^(-1/2).*(x+(-5).*y+3.*2.^(1/2).*z));
    B = hx+hy*exp(1).^(sqrt(-1).*2.^(-1/2).*((-1).*x+(-3).*y+(-1).*2.^(1/2).*z));
    
R(:,:,1,1) = (hz/2); R(:,:,1,2) = 0; R(:,:,1,3) = (1/4).*(A+(-1).*B); R(:,:,1,4) = (1/4).*((-1).*A+(-1).*B);
R(:,:,2,1) = 0; R(:,:,2,2) = (-hz/2); R(:,:,2,3) = (1/4).*(A+B); R(:,:,2,4) = (1/4).*((-1).*A+B);
R(:,:,3,1) = (1/4).*(conj(A)+(-1).*conj(B)); R(:,:,3,2) = (1/4).*(conj(A)+conj(B)); R(:,:,3,3) = (hz/2); R(:,:,3,4) = 0;
R(:,:,4,1) = (1/4).*((-1).*conj(A)+(-1).*conj(B)); R(:,:,4,2) = (1/4).*((-1).*conj(A)+conj(B)); R(:,:,4,3) = 0; R(:,:,4,4) = (-hz/2);
 
%Matrix Multiplication
    for m=1:4
        norm(:,:,m) = sqrt( abs(U(:,:,m,1)).^2 + abs(U(:,:,m,2)).^2 + abs(U(:,:,m,3)).^2 + abs(U(:,:,m,4)).^2 );
        
        for n=1:4
            U(:,:,m,n) = U(:,:,m,n)./norm(:,:,m);
        end
    end
    
    %fixing relative phasese of the f and f^\dagger equations
    phase2 = ( U(:,:,1,1)./abs( U(:,:,1,1) ) ) .* ( abs( U(:,:,2,2) )./U(:,:,2,2) );
    phase4 = ( U(:,:,3,3)./abs( U(:,:,3,3) ) ) .* ( abs( U(:,:,4,4) )./U(:,:,4,4) );
    for n=1:4
        U(:,:,2,n) = U(:,:,2,n) .* phase2;
        U(:,:,4,n) = U(:,:,4,n) .* phase4;
    end
    
    %Finding R in diagonal space (of H), call it V 
    V=zeros(L,L,4,4);
    for n=1:4
        for m=1:4
            for k=1:4
                V(:,:,n,m) = V(:,:,n,m) + conj(U(:,:,n,k)).*R(:,:,k,m);
            end
        end
    end
    
     VV=zeros(L,L,4,4);
     for n=1:4
         for m=1:4
             for k=1:4
                 VV(:,:,n,m)=VV(:,:,n,m) + V(:,:,n,k).*U(:,:,m,k);
             end
         end
     end
     
     
     %%%%%Do it again for a second set of params to look for mixing of
     %%%%%irreps
     
hx = hx2;
hy=  hy2;
hz = hz2;


    A = hx*exp(1).^(sqrt(-1).*2.^(-1/2).*(4.*x+(-4).*y+2.*2.^(1/2).*z))...
        + hy*exp(1).^(sqrt(-1).*2.^(-1/2).*(x+(-5).*y+3.*2.^(1/2).*z));
    B = hx+hy*exp(1).^(sqrt(-1).*2.^(-1/2).*((-1).*x+(-3).*y+(-1).*2.^(1/2).*z));
    
R(:,:,1,1) = (hz/2); R(:,:,1,2) = 0; R(:,:,1,3) = (1/4).*(A+(-1).*B); R(:,:,1,4) = (1/4).*((-1).*A+(-1).*B);
R(:,:,2,1) = 0; R(:,:,2,2) = (-hz/2); R(:,:,2,3) = (1/4).*(A+B); R(:,:,2,4) = (1/4).*((-1).*A+B);
R(:,:,3,1) = (1/4).*(conj(A)+(-1).*conj(B)); R(:,:,3,2) = (1/4).*(conj(A)+conj(B)); R(:,:,3,3) = (hz/2); R(:,:,3,4) = 0;
R(:,:,4,1) = (1/4).*((-1).*conj(A)+(-1).*conj(B)); R(:,:,4,2) = (1/4).*((-1).*conj(A)+conj(B)); R(:,:,4,3) = 0; R(:,:,4,4) = (-hz/2);
 
    V=zeros(L,L,4,4);
    for n=1:4
        for m=1:4
            for k=1:4
                V(:,:,n,m) = V(:,:,n,m) + conj(U(:,:,n,k)).*R(:,:,k,m);
            end
        end
    end
    
     VV2=zeros(L,L,4,4);
     for n=1:4
         for m=1:4
             for k=1:4
                 VV2(:,:,n,m)=VV2(:,:,n,m) + V(:,:,n,k).*U(:,:,m,k);
             end
         end
     end
     
     
    
    %The Delta weight in the Raman term (turned into column vectors)
    Wpp = pi*real( VV(:,:,1,2).*VV2(:,:,2,1) + VV(:,:,2,1).*VV2(:,:,1,2) )/2;%/factor2;
    Wmm = pi*real( VV(:,:,3,4).*VV2(:,:,4,3) + VV(:,:,4,3).*VV2(:,:,3,4) )/2;%/factor2;
    Wpm = pi*(real( VV(:,:,1,4).*VV2(:,:,4,1) + VV(:,:,4,1).*VV2(:,:,1,4) ) ...
        + real( VV(:,:,3,2).* VV2(:,:,2,3) + VV(:,:,2,3).* VV2(:,:,3,2) ) )/2;%factor2;
          
    [histw, histv] = histwv(2*ep(sn>=0),Wpp(sn>=0),0,emax,bins);
    Ipp = Ipp + histw;
    DEp = DEp + histv;
    
    [histw, histv] = histwv(2*em(sn>=0),Wmm(sn>=0),0,emax,bins);
    Imm = Imm + histw;
    DEm = DEm + histv;
    
    [histw, histv] = histwv(ep(sn>=0)+em(sn>=0),Wpm(sn>=0),0,emax,bins);
    Ipm = Ipm + histw;
    DE  = DE  + histv;
      
    %Add one to the number of BZ points counted if the point is indeed in
    %the BZ
    Z=Z+size(sn(sn(:)>=0),1);
    %Z should approximately be 18/29*L^3;
end

%Normalize the results by the BZ volume

DEp=DEp*bins/(Z*emax );
DEm=DEm*bins/(Z*emax );
DE =DE *bins/(Z*emax );

Ipp=Ipp*bins/(Z*2*emax );
Ipm=Ipm*bins/(Z*2*emax );
Imm=Imm*bins/(Z*2*emax );

out = [Ev,DEp,DEm,DE,Ipp,Imm,Ipm];

%Note that Intensity should be plotted against Ev while DOS against Ev/2

%DEp(1)=0; DEm(1)=0;

%Plot
%figure;
%hold on;
%plot(Ev,Ipp,Ev,Imm,Ev,Ipm,Ev,Ipp+Imm+Ipm);
%title(['Raman Spectrum for HyperHoneycomb Kitaev spinons: bins = ',num2str(bins),', N= ',num2str(N)])
%xlabel('E');
%ylabel('I(E)');
%hold off;