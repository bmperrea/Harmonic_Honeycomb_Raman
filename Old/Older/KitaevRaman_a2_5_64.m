function out = KitaevRaman_a2_5_64(N,bins,flag,nn,jx,jz,h)
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

%The size of the matrix 
T=64;
M = 4*T;
S = 2*T;

N=N/sqrt(T);

%the max energy for the energy axis
emax = 15*max(jx,.5) + h;

%The values of the Kitaev couplings
jy = jx;  %Jy is always same as Jx
%Jz = 1; 
    %The 'kappa' terms characterizing the B-field term
kx = 1; ky=kx; kz=kx;
kx = kx*h;  ky = ky*h;  kz = kz*h;


%A convenient notation
L=2*N;

%Initialization of Arrays

Ev = (1:bins)'*emax/bins;

r1 = rand; r2 = rand;% r3 = rand;

pts = (-N):(N-1);

a1 = [sqrt(3),0,0]';
%a2 = [1,0,3]';
a3 = [-1/sqrt(3),sqrt(29/3),0]';

%The six are the six symmetric binary products of R^aa, R^ac, R^cc
    %h1 = {aa, aa, aa, ac, ac, cc};
    %h2 = {aa, ac, cc, ac, cc, cc};
m1 = {1,1,1,2,2,3};
m2 = {1,2,3,2,3,3};
%They appear as {+,0,-,+,0,+}
%Note   I5 = -3I2 ;  I6/9 = I1 = I3/(-3)
%3/1 6/1 5/2 = -3, 9 , -3

hx = {1/2,1/2,1/2};
hy = {1/2,-1/2,1/2};
hz = {0,0,2};

%Initializtion 

zer = zeros(2);
ey  = eye(2);
gamma = [zer, ey; ey, zer];

R = cell(1,3); R1=R; R2=R;
Re=R; %RR=R; RR1=R; RR2=R; 

I = cell(1,6); 
zerolist = zeros(size(Ev));
DD = zerolist; DD2 = DD;

W = cell(1,6); 
zeroB = zeros(L,S,S);
enty = zeroB;

RRr = cell(1,3); RRr1=RRr; RRr2=RRr;

for m=1:6
    I{m} = zerolist;  
    W{m} = zeroB;
end
I2 = I; I1 = I;
W2 = W; W1 = W;

    %With matrices larger than 4, diagonalization must be done with an
    %iterative algorithm so that vectorization is not an option. Therefore,
    %we loop in this code.
for xind = pts
    
    
    %I display an estimate of the end time
    aleph = 1 + round( (10/N)^2 ); %iterations for ~1 second of computation, or one iteration, which ever takes longer
    %Or an iteration, which ever is longer
    if (xind == -N+2) && flag   %the first few may be slower due to initialization, start from the third one
        cl1 = clock; 
    end
    if (xind == -N+2+aleph) && flag
        cl2 = clock;
        time = (cl2-cl1)*(2*N)/aleph*nn;
        %cl = cl1 + time*(2*N-2)/(2*N);
        
        format shortg
        disp('approximate time to take:')
        disp( datestr(time(6)/24/3600, 'DD-HH:MM:SS') )
        format
    end
    
    
    for yind = pts
        x = (xind + r1)*pi/(sqrt(3)*N);
        y = sqrt(3/29)*pi*(yind + r2)/N + x/sqrt(29); 
    
    
    
             
    p1 = exp(1i*(x*a1(1)+y*a1(2))); 
    %p2 = exp(-1i*(x*a2(1)+y*a2(2)+z*a2(3))); 
    p3 = exp(1i*(x*a3(1)+y*a3(2)));

    A = conj(p3) * (jx + jy*p1);
    B = jx;
    d = ky*(1-conj(p3)) + kx*( p1*conj(p3) );
    a = -kz*( p1-conj(p1) );
    b = 0;%kz*( p2 - conj(p2) );
    %Note that a and b should be purely imaginary
    B2 = jy;
    d2 = +kx;
    b2 = -kz;
    
%     H = [1i*a-1i*b+2*jz,A+conj(B),1i*(a+b),-A+2i.*d+conj(B);...
%             B+conj(A),-1i*a+1i*b+2*jz,-B+conj(A)-2i*conj(d),1i*(a+b);...
%             1i*(a+b),A+2i*d-conj(B),1i*(a-b+2i*jz),-A-conj(B);...
%             B-conj(A)-2i*conj(d),1i*(a+b),-B-conj(A),1i*(-a+b+2i*jz)];
    
    H1 =    [1i*a-1i*b+2*jz,A+conj(B),1i*(a+b),-A+2i*d+conj(B);...
            B+conj(A),-1i*a+1i*b+2*jz,-B+conj(A)-2i*conj(d),1i*(a+b);...
            1i*(a+b),A+2i*d-conj(B),1i*(a-b+2i*jz),-A-conj(B);...
            B-conj(A)-2i*conj(d),1i*(a+b),-B-conj(A),1i*(-a+b+2i*jz)];
        
   H2 =    [-1i*b2,conj(B2),1i*(b2),2i*d2+conj(B2);...
            B2,1i*b2,-B2-2i*conj(d2),1i*(b2);...
            1i*(b2),2i*d2-conj(B2),-1i*(b2),-conj(B2);...
            B2-2i*conj(d2),1i*(b2),-B2,1i*(b2)];
    
   
    zero = zeros(4);
    
    
     Hd = [H1, H2, zero, zero, zero, zero, zero, zero;...
         H2', H1, H2, zero, zero, zero, zero, zero; ...
         zero, H2', H1, H2, zero, zero, zero, zero; ...
         zero, zero, H2', H1, H2, zero, zero, zero; ...
         zero, zero, zero, H2', H1, H2, zero, zero; ...
         zero, zero, zero, zero, H2', H1, H2, zero; ...
         zero, zero, zero, zero, zero, H2', H1, H2; ...
         zero, zero, zero, zero, zero, zero, H2', H1]/2;
     
     Hd2 = [zero, zero, zero, zero, zero, zero, zero, zero;...
            zero, zero, zero, zero, zero, zero, zero, zero; ...
            zero, zero, zero, zero, zero, zero, zero, zero; ...
            zero, zero, zero, zero, zero, zero, zero, zero; ...
            zero, zero, zero, zero, zero, zero, zero, zero; ...
            zero, zero, zero, zero, zero, zero, zero, zero; ...
            zero, zero, zero, zero, zero, zero, zero, zero; ...
              H2, zero, zero, zero, zero, zero, zero, zero]/2;
     
   z2 = [zero, zero, zero, zero, zero, zero, zero, zero;...
         zero, zero, zero, zero, zero, zero, zero, zero; ...
         zero, zero, zero, zero, zero, zero, zero, zero; ...
         zero, zero, zero, zero, zero, zero, zero, zero; ...
         zero, zero, zero, zero, zero, zero, zero, zero; ...
         zero, zero, zero, zero, zero, zero, zero, zero; ...
         zero, zero, zero, zero, zero, zero, zero, zero; ...
         zero, zero, zero, zero, zero, zero, zero, zero];
     
      H = [Hd, Hd2, z2, z2, z2, z2, z2, z2;...
         Hd2', Hd, Hd2, z2, z2, z2, z2, z2; ...
         z2, Hd2', Hd, Hd2, z2, z2, z2, z2; ...
         z2, z2, Hd2', Hd, Hd2, z2, z2, z2; ...
         z2, z2, z2, Hd2', Hd, Hd2, z2, z2; ...
         z2, z2, z2, z2, Hd2', Hd, Hd2, z2; ...
         z2, z2, z2, z2, z2, Hd2', Hd, Hd2; ...
         z2, z2, z2, z2, z2, z2, Hd2', Hd];
        
    [V,D] = eig(H);
    
  
   
    %Sort them to ascending eigenvalue order for the first half
    [~,II1] = sort(real(diag(D)));
    [~,II2] = sort(real(diag(D)),'descend');
    fh = 1:(size(H,1)/2);
    II = [II2(fh),II1(fh)];    
    V = V(:, II);
    %en = temp(fh);
        
    %Need to normalize V to get U
    U = V./repmat( sqrt(sum(V.*conj(V),1)), [size(V,1),1]);
    %U(~isfinite(U)) = 1; This is for error handling
           
    %Raman operators (They are 3 deep in cell space)
    jzc = jz*hz;
    Ac = conj(p3) * (jx*hx + hy*jy*p1);
    Bc = jx*hy;
    dc = 0*hz;
    ac = 0*hz;
    bc = 0*hz;%kz*( p2 - conj(p2) );
    
    B2c = jy*hy;
    d2c = 0*hz;
    b2c = 0*hz;
    
    zeroc = zero*hz; %makes a cell array of zero matrices
    
    for n =1:3
    
        R1{n} = [1i*ac{n}-1i*bc{n}+2*jzc{n},Ac{n}+conj(Bc{n}),1i*(ac{n}+bc{n}),-Ac{n}+2i.*dc{n}+conj(Bc{n});...
            Bc{n}+conj(Ac{n}),-1i*ac{n}+1i*bc{n}+2*jzc{n},-Bc{n}+conj(Ac{n})-2i*conj(dc{n}),1i*(ac{n}+bc{n});...
            1i*(ac{n}+bc{n}),Ac{n}+2i*dc{n}-conj(Bc{n}),1i*(ac{n}-bc{n}+2i*jzc{n}),-Ac{n}-conj(Bc{n});...
            Bc{n}-conj(Ac{n})-2i*conj(dc{n}),1i*(ac{n}+bc{n}),-Bc{n}-conj(Ac{n}),1i*(-ac{n}+bc{n}+2i*jzc{n})]; 

 R2{n} = [-1i*b2c{n},conj(B2c{n}),1i*(b2c{n}),2i*d2c{n}+conj(B2c{n});...
            B2c{n},1i*b2c{n},-B2c{n}-2i*conj(d2c{n}),1i*(b2c{n});...
            1i*(b2c{n}),2i*d2c{n}-conj(B2c{n}),1i*(b2c{n}),-conj(B2c{n});...
            B2c{n}-2i*conj(d2c{n}),1i*(b2c{n}),-B2c{n},1i*(b2c{n})];

        zero = zeros(4);

        %R{n} =  R1{n}/2;%[R1{n}, R2{n}, zeroc{n}, zeroc{n}; ...
            % R2{n}', R1{n}, R2{n}, zeroc{n}; ...
            % zeroc{n}, R2{n}', R1{n}, R2{n}; ...
            % zeroc{n}, zeroc{n}, R2{n}' , R1{n}]/2;
            
   Hd = [R1{n}, R2{n}, zero, zero, zero, zero, zero, zero;...
       R2{n}', R1{n}, R2{n}, zero, zero, zero, zero, zero; ...
       zero, R2{n}', R1{n}, R2{n}, zero, zero, zero, zero; ...
       zero, zero, R2{n}', R1{n}, R2{n}, zero, zero, zero; ...
       zero, zero, zero, R2{n}', R1{n}, R2{n}, zero, zero; ...
       zero, zero, zero, zero, R2{n}', R1{n}, R2{n}, zero; ...
       zero, zero, zero, zero, zero, R2{n}', R1{n}, R2{n}; ...
       zero, zero, zero, zero, zero, zero, R2{n}', R1{n}]/2;
   
   Hd2 = [zero, zero, zero, zero, zero, zero, zero, zero;...
          zero, zero, zero, zero, zero, zero, zero, zero; ...
          zero, zero, zero, zero, zero, zero, zero, zero; ...
          zero, zero, zero, zero, zero, zero, zero, zero; ...
          zero, zero, zero, zero, zero, zero, zero, zero; ...
          zero, zero, zero, zero, zero, zero, zero, zero; ...
          zero, zero, zero, zero, zero, zero, zero, zero; ...
         R2{n}, zero, zero, zero, zero, zero, zero, zero]/2;
        
     R{n} = [Hd, Hd2, z2, z2, z2, z2, z2, z2;...
         Hd2', Hd, Hd2, z2, z2, z2, z2, z2; ...
         z2, Hd2', Hd, Hd2, z2, z2, z2, z2; ...
         z2, z2, Hd2', Hd, Hd2, z2, z2, z2; ...
         z2, z2, z2, Hd2', Hd, Hd2, z2, z2; ...
         z2, z2, z2, z2, Hd2', Hd, Hd2, z2; ...
         z2, z2, z2, z2, z2, Hd2', Hd, Hd2; ...
         z2, z2, z2, z2, z2, z2, Hd2', Hd];

        %For the slab we can also define an edge Raman operator
        He = [R1{n}, zero, zero, zero, zero, zero, zero, zero;...
          zero, zero, zero, zero, zero, zero, zero, zero; ...
          zero, zero, zero, zero, zero, zero, zero, zero; ...
          zero, zero, zero, zero, zero, zero, zero, zero; ...
          zero, zero, zero, zero, zero, zero, zero, zero; ...
          zero, zero, zero, zero, zero, zero, zero, zero; ...
          zero, zero, zero, zero, zero, zero, zero, zero; ...
          zero, zero, zero, zero, zero, zero, zero, zero]/2;
      
              He2 = [zero, zero, zero, zero, zero, zero, zero, zero;...
          zero, zero, zero, zero, zero, zero, zero, zero; ...
          zero, zero, zero, zero, zero, zero, zero, zero; ...
          zero, zero, zero, zero, zero, zero, zero, zero; ...
          zero, zero, zero, zero, zero, zero, zero, zero; ...
          zero, zero, zero, zero, zero, zero, zero, zero; ...
          zero, zero, zero, zero, zero, zero, zero, zero; ...
          zero, zero, zero, zero, zero, zero, zero, R1{n}]/2;
        
 Re{n} = [He, z2, z2, z2, z2, z2, z2, z2;...
           z2, z2, z2, z2, z2, z2, z2, z2; ...
           z2, z2, z2, z2, z2, z2, z2, z2; ...
           z2, z2, z2, z2, z2, z2, z2, z2; ...
           z2, z2, z2, z2, z2, z2, z2, z2; ...
           z2, z2, z2, z2, z2, z2, z2, z2; ...
           z2, z2, z2, z2, z2, z2, z2, z2; ...
           z2, z2, z2,z2, z2, z2, z2, He2];
        
    end
        
 
    
    %Finding R in diagonal space (of H), call it RR
    RR = U'*R*U;
    RR1 = U'*Re*U;
%    RR2 = U'*Re1*U;
    
    D2 = U'*H*U;
    
%     %Check that we diagonalized the matrix
%     
%     Errs = abs(D2 - D);
%     err = max(max(max(max(Errs))));
%     if err>10^(-6)
%         disp(err)
%     end
        
    %The relevant energies
    en = diag(D2);
    en2 = en(1:S);
    
    if emax/2 < max(en2)
        disp([max(en),emax,jx,jz])
    end
    %Check that the energies are positive
    [me,f]=min(en2);
    if me<0
        disp([me,f,jx,jz])
    end
    
    
    %Take out the relevant chunks of the Raman operator (SC terms)
    for n=1:3
        RRr{n} = RR{n}(S+1:M,1:S);
        RRr1{n} = RR1{n}(S+1:M,1:S);
 %       RRr2{n} = RR2{n}(S+1:M,1:S);
    end
    
    
    ent = repmat(en2,1,S) + repmat(en2',S,1);
    

   
    %Check that emax is going to work
    
    %The six binary combinations from above.
    for m =1:6;
        W{m}(yind+N+1,:,:) = pi*real(conj(RRr{m1{m}}).*RRr{m2{m}} ...
            + conj(RRr{m2{m}}).*RRr{m1{m}});
        W1{m}(yind+N+1,:,:) = pi*real(conj(RRr1{m1{m}}).*RRr1{m2{m}} ...
            + conj(RRr1{m2{m}}).*RRr1{m1{m}});
 %       W2{m}(xind+N+1,yind+N+1,:,:) = pi*real(RRr2{m1{m}}.*RRr2{m2{m}} ...
 %           + RRr2{m2{m}}.*RRr2{m1{m}});
    end
            
    enty(yind+N+1,:,:) = ent;              
   
    end
    
    
    %histogram it
    for m=1:6
    %R
    [histw, histv] = histwv(2*enty,W{m},0,emax,bins);
    I{m} = I{m} + histw;
    DD = DD + histv;
    
    %Re1
    [histw, histv] = histwv(2*enty,W1{m},0,emax,bins);
    I1{m} = I1{m} + histw;
    DD2 = DD2 + histv;
    
    %Re2
  %  [histw, histv] = histwv(2*enty,W2{m},0,emax,bins);
  %  I2{m} = histw;
    %Dpp = histv;
    end
    
end

%Normalize the results by the BZ volume
%Vm = sqrt(29)/2.913;
%Z = S^2*L^2*Vm/(2*pi)^2;

%DD=DD*bins./(Z*emax );

dd = sum(DD)*emax/bins;
DD = DD/dd;
DD2=DD2/dd;

I=I*bins./dd;%(Z*2*emax );
I1=I1*bins./dd;%(Z*2*emax );
%I2=I2*bins./(Z*2*emax );


out = cell(1,8*3-2);

out{1} = Ev;
out{2} = [DD,DD2];
for m=1:6
    out{m+2} = [I{m},I1{m}];
end

%Note that Intensity should be plotted against Ev while DOS against Ev/2
%(otherwise it is the two-particle DOS

%DEp(1)=0; DEm(1)=0;

% Plot
% figure;
% hold on;
% plot(Ev,Ipp,Ev,Imm,Ev,Ipm,Ev,Ipp+Imm+Ipm);
% title(['Raman Spectrum for HyperHoneycomb Kitaev spinons: bins = ',num2str(bins),', N= ',num2str(N)])
% xlabel('E');
% ylabel('I(E)');
% hold off;