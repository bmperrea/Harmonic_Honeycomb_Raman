function out = Raman2D(N,bins,type,flag,nn)
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
Jx = 1; 
Jy = Jx;  %Jy is always same as Jz
Jz = 1;  %Jz is always 1 the way I coded it

%pick the right Raman operator
switch type
    case 'xx'
        hx = 3/4*Jx;          hy = 3/4*Jy ; hz = 0;
    case 'yy'
        hx = 1/2*Jx;          hy = 1/2*Jy;  hz = Jz;
    case 'xy'
        hx = -sqrt(3)/4*Jx;  hy = sqrt(3)/4*Jy;  hz = 0;  
end


L=2*N;

%hhz = repmat(hz,L,L);

%The following can be used to ensure that H is diagonalized by the unitary
%matrix below by checking that the Raman spectra are zero.
%hx=Jx;hy=Jy;hz=Jz;

emin  = 0;
emax  = 6*max(1,Jx); er = emax-emin;

%initialization
%Evalues=linspace(0,emax,bins);

DE=zeros(bins,1);

I=zeros(bins,1);


%z=zeros(L,L); ep=z; em=ep;
%r(1)=-1;r(2)=-1;r(3)=1;r(4)=1;
%q=(1:(L)^2)'; 

Ev = (1:bins)'*emax/(bins*max(1,Jx));

r1 = rand; r2 = rand; %r3 = rand;

pts = (-N):(N-1); %pts=pts/N;
%treat x and y values as corresponding to column
%ans row values respectively
x=(repmat(pts,L,1)+r1)/N*2*pi/sqrt(3);
%z=repmat((pts.' + r3)/N,1,L)*pi/(3*sqrt(2));

%U=zeros(L,4,4); %V=U; VV=U; Udag = U; 
%R=U;
%norm=zeros(L,L,4);
%en = zeros(L,L,4);

Z=0;

%for ind = pts 
    y = (repmat(pts.',1,L) + r2)*2*pi/(sqrt(3)*N);
    
    
    %I display an estimate of the end time
%     aleph = 1 ;%+ round( (200/N)^2 ); %iterations for a few seconds of computation
    %Or an iteration, which ever is longer
%     if ind == -N+2 && flag   %the first few may be slower due to initialization
%         cl1 = clock; 
%     end
%     if ind == -N+2+aleph && flag
%         cl2 = clock;
%         time = (cl2-cl1)*(2*N)/aleph *3*nn;
%         %cl = cl1 + time*(2*N-2)/(2*N);
%         
%         format shortg
%         disp('approximate time to take:')
%         disp( datestr(time(6)/24/3600, 'DD-HH:MM:SS') )
%     %    disp(time(6))
%    %     disp('approximate end datetime:')
%    %     disp(cl)
%         format
%     end
        
        
           
 %   if minute ~= cl(5)
        
 %       disp( ind )  
 %       minute = cl(5);
 %   end    
    
%    tic
    %If this value is not positive I negate the energies out of the hist
    sn = sign( 2*pi/sqrt(3) - abs(y) - abs(x)/sqrt(3) );
  %  [L L2]=size(sn(:)>=0);
    
    %Now I do arithmetic as if x,y, and z were numbers but with x and a
    %matrices and doing only elementwise operations with them.
    
    %I only store values for a given y and then reduce the information into
    %a histogram before looping back to keep from overloading the memory
    
    %My computer handles about 2 to 5 * 10^8 doubles at once in memory so I
    %expect this code to max out the memory around N~10^4
    
    
    %The Eigenvalues
        
 %   toc 
  %  tic
    s = Jz + Jx*exp(1).^(1i*(y+x*sqrt(3))/2) + Jy*exp(1).^(1i*(y-x*sqrt(3))/2) ;
    h = hz + hx*exp(1).^(1i*(y+x*sqrt(3))/2) + hy*exp(1).^(1i*(y-x*sqrt(3))/2) ;
    en = 2*abs(s);
 
    V = (imag(h.*conj(s))./abs(s)).^2;

    %The Delta weight in the Raman term (turned into column vectors)
    W = 4*pi*abs(V).^2;

    [histw, histv] = histwv(en(sn>=0).',W(sn>=0).',0,emax,bins);
    I = I + histw;
    DE = DE + histv;
      
  %  toc
    
    Z=Z+size(sn(sn(:)>=0),1);
    
%end

%Z = 18/29*L^3;
%Z = (pi^3)/24; 
DE =DE *bins/(Z*er );
I=I*bins/(Z*er );

out = [Ev,DE,I];