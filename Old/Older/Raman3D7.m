function out = Raman3D7(N,bins,type,flag,nn,Jx)
%output form: out = [Ev,DEp,DEm,DE,Ipp,Imm,Ipm]


%tic

%input-parameter:
%2N is the number of points in a dimension
%flag is a boolean telling whether to estimate time on this run
%nn is the number of times that you plan to do a similar calculation (for error estimation) 
%type is a string specifying which Raman operator


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

%bins = 2*N;

%input data
%Jx = 1; 
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
end
hx = hx*Jx; hy = hy*Jy; hz = hz*Jz;

L=2*N;

%hhz = repmat(hz,L,L);

%The following can be used to ensure that H is diagonalized by the unitary
%matrix below by checking that the Raman spectra are zero.
%hx=Jx;hy=Jy;hz=Jz;

%epmin = 1/2;%0.507722357;
%epmax = 3/2; erp = epmax-epmin;
%emmin = 0; 
%emmax = sqrt(5)/2; erm = emmax-emmin;
%epmmin  = 1;
%epmmax  = max(1,Jx)*sqrt(5); erpm = epmmax-epmmin;

emax = 3*max(Jx,1) + .2;

%initialization
%Evalues=linspace(0,emax,bins);

DEp=zeros(bins,1);
DEm=zeros(bins,1);
DE=zeros(bins,1);

Ipp=zeros(bins,1);
Imm=zeros(bins,1);
Ipm=zeros(bins,1);


%z=zeros(L,L); ep=z; em=ep;
%r(1)=-1;r(2)=-1;r(3)=1;r(4)=1;
%q=(1:(L)^2)'; 

%Evp = epmin + (1:bins)'*erp/bins;
%Evm = emmin + (1:bins)'*erm/bins;
%Evpm  = epmmin  + (1:bins)'*erpm/bins;
Ev = (1:bins)'*emax/bins;

r1 = rand; r2 = rand; r3 = rand;

pts = (-N):(N-1); %pts=pts/N;
%treat x and y values as corresponding to column
%ans row values respectively
x=repmat((pts+r1)/N,L,1)*pi/2;
z=repmat((pts.' + r3)/N,1,L)*pi/(3*sqrt(2));

U=zeros(L,L,4,4); %V=U; VV=U; Udag = U; 
R=U;
norm=zeros(L,L,4);
%en = zeros(L,L,4);

Z=0;

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
    %    disp(time(6))
   %     disp('approximate end datetime:')
   %     disp(cl)
        format
    end
        
        
           
 %   if minute ~= cl(5)
        
 %       disp( ind )  
 %       minute = cl(5);
 %   end    
    
%    tic
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

    %the following is code converted with ToMatlab
    A = exp(1).^(sqrt(-1).*2.^(-1/2).*(4.*x+(-4).*y+2.*2.^(1/2).*z))...
        + exp(1).^(sqrt(-1).*2.^(-1/2).*(x+(-5).*y+3.*2.^(1/2).*z));
    B = 1+exp(1).^(sqrt(-1).*2.^(-1/2).*((-1).*x+(-3).*y+(-1).*2.^(1/2).*z));
    A = A*Jx; B=B*Jx;
    
    a=A.*conj(A)+B.*conj(B);
    b=1+A.*B.*conj(A).*conj(B)+2.*real(B.*conj(A));
    c=(2+a).^2+(-4).*b;
    ep = (1/2).*2.^(-1/2).*(2+a+c.^(1/2)).^(1/2);
    em = (1/2).*2.^(-1/2).*(2+a+(-1).*c.^(1/2)).^(1/2);
    dp = c.^(1/2).*(conj(A)+conj(B))+(conj(A)+(-1).* ...
        conj(B)).*(A.*conj(A)+(-1).*B.*conj(B));
    dm = c.^(1/2).*(conj(A)+conj(B))+(-1).*(conj(A)+(-1).* ...
        conj(B)).*(A.*conj(A)+(-1).*B.*conj(B));
    
 %   toc 
  %  tic
    
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
 

    A = hx*exp(1).^(sqrt(-1).*2.^(-1/2).*(4.*x+(-4).*y+2.*2.^(1/2).*z))...
        + hy*exp(1).^(sqrt(-1).*2.^(-1/2).*(x+(-5).*y+3.*2.^(1/2).*z));
    B = hx+hy*exp(1).^(sqrt(-1).*2.^(-1/2).*((-1).*x+(-3).*y+(-1).*2.^(1/2).*z));
    
R(:,:,1,1) = (hz/2); R(:,:,1,2) = 0; R(:,:,1,3) = (1/4).*(A+(-1).*B); R(:,:,1,4) = (1/4).*((-1).*A+(-1).*B);
R(:,:,2,1) = 0; R(:,:,2,2) = (-hz/2); R(:,:,2,3) = (1/4).*(A+B); R(:,:,2,4) = (1/4).*((-1).*A+B);
R(:,:,3,1) = (1/4).*(conj(A)+(-1).*conj(B)); R(:,:,3,2) = (1/4).*(conj(A)+conj(B)); R(:,:,3,3) = (hz/2); R(:,:,3,4) = 0;
R(:,:,4,1) = (1/4).*((-1).*conj(A)+(-1).*conj(B)); R(:,:,4,2) = (1/4).*((-1).*conj(A)+conj(B)); R(:,:,4,3) = 0; R(:,:,4,4) = (-hz/2);
 
%toc
%tic

    for m=1:4
        norm(:,:,m) = sqrt( abs(U(:,:,m,1)).^2 + abs(U(:,:,m,2)).^2 + abs(U(:,:,m,3)).^2 + abs(U(:,:,m,4)).^2 );
        
        for n=1:4
            %R(:,:,n,m) = (R(:,:,n,m)+conj(R(:,:,m,n)))/4;
            U(:,:,m,n) = U(:,:,m,n)./norm(:,:,m);
            %Udag(:,:,m,n)= conj( U(:,:,n,m) );
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
    
    %The Delta weight in the Raman term (turned into column vectors)
    Wpp = 2*pi*abs(VV(:,:,1,2)).^2;
    Wmm = 2*pi*abs(VV(:,:,3,4)).^2;
    Wpm = 2*pi*(abs(VV(:,:,1,4)).^2 + abs(VV(:,:,3,2)).^2);
          
    [histw, histv] = histwv(2*ep(sn>=0),Wpp(sn>=0),0,emax,bins);
    Ipp = Ipp + histw;
    DEp = DEp + histv;
    
    [histw, histv] = histwv(2*em(sn>=0),Wmm(sn>=0),0,emax,bins);
    Imm = Imm + histw;
    DEm = DEm + histv;
    
    [histw, histv] = histwv(ep(sn>=0)+em(sn>=0),Wpm(sn>=0),0,emax,bins);
    Ipm = Ipm + histw;
    DE  = DE  + histv;
      
  %  toc
    
    Z=Z+size(sn(sn(:)>=0),1);
    
end

%Z = 18/29*L^3;
%Z = (pi^3)/24; 
DEp=DEp*bins/(Z*emax );
DEm=DEm*bins/(Z*emax );
DE =DE *bins/(Z*emax );

Ipp=Ipp*bins/(Z*2*emax );
Ipm=Ipm*bins/(Z*2*emax );
Imm=Imm*bins/(Z*2*emax );

out = [Ev,DEp,DEm,DE,Ipp,Imm,Ipm];

%Note that I should be plotted against Ev while DE against Ev/2

%DEp(1)=0; DEm(1)=0;

%Plot
%figure;
%hold on;
%plot(Ev,Ipp,Ev,Imm,Ev,Ipm,Ev,Ipp+Imm+Ipm);
%title(['Raman Spectrum for HyperHoneycomb Kitaev spinons: bins = ',num2str(bins),', N= ',num2str(N)])
%xlabel('E');
%ylabel('I(E)');
%hold off;

%toc