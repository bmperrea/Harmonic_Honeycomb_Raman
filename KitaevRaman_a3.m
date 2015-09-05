function out = KitaevRaman_a2(N,bins,flag,nn,jx,jz,h)
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
T=4;
M = 4*T;
S = 2*T;

%the max energy for the energy axis
emax = 7*max(jx,.5) + .2 + h;

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
Re1=R; Re2=R; %RR=R; RR1=R; RR2=R; 

I = cell(1,6); 
zerolist = zeros(size(Ev));
D = zerolist;

W = cell(1,6); 
zeroB = zeros(L,L,S,S);
enty = zeroB;

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
    aleph = 1 + round( (70/N)^2 ); %iterations for ~1 second of computation, or one iteration, which ever takes longer
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
        y = sqrt(3/29)*pi*(yind + r2) + x/sqrt(29); 
    
    
    
             
    p1 = exp(1i*(x*a1(1)+y*a1(2))); 
    %p2 = exp(-1i*(x*a2(1)+y*a2(2)+z*a2(3))); 
    p3 = exp(1i*(x*a3(1)+y*a3(2)));

    A = conj(p3) * (jx + jy*p1);
    B = jx;
    d = ky*(1-conj(p3)) + kx*( p1*conj(p3) );
    a = -kz*( p1-conj(p1) );
    b = 0;%kz*( p2 - conj(p2) );
    
    B2 = jy;
    d2 = -kx;
    b2 = kz;
    
%     H = [2*jz,A+conj(B),2i*a,-A+2i*d+conj(B); ...
%             B+conj(A),2.*jz,-B+conj(A)-2i*conj(d),2i*b;...
%             2i*a,A+2i*d-conj(B),-2*jz,-A-conj(B);...
%             B-conj(A)-2i*conj(d),2i*b,-B-conj(A),-2*jz];
    
    H1 = [2*jz,A+conj(B),2i*a,-A+2i*d+conj(B); ...
            B+conj(A),2.*jz,-B+conj(A)-2i*conj(d),2i*b;...
            2i*a,A+2i*d-conj(B),-2*jz,-A-conj(B);...
            B-conj(A)-2i*conj(d),2i*b,-B-conj(A),-2*jz]; 
        
    H2 = [0,conj(B2),0,2i*d2+conj(B2); ...
            B2,0,-B2-2i*conj(d2),2i*b2;...
            0,2i*d2-conj(B2),-conj(B2);...
            B2-2i*conj(d2),2i*b2,-B2,0];
    
   
    zero = zeros(4);
    
    H = [H1, H2, zero, zero; H2, H1, H2, zero; ...
            zero, H2, H1, H2; zero, zero, H2 , H1]/4;
        
    [V,D] = eig(H);
    
  
   
    %Sort them to ascending eigenvalue order for the first half
    [~,II1] = sort(diag(D));
    [~,II2] = sort(diag(D),'descend');
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
    
        R1{n} = [2*jzc{n},Ac{n}+conj(Bc{n}),2i*ac{n},-Ac{n}+2i*dc{n}+conj(Bc{n}); ...
            Bc{n}+conj(Ac{n}),2.*jzc{n},-Bc{n}+conj(Ac{n})-2i*conj(dc{n}),2i*bc{n};...
            2i*ac{n},Ac{n}+2i*dc{n}-conj(Bc{n}),-2*jzc{n},-Ac{n}-conj(Bc{n});...
            Bc{n}-conj(Ac{n})-2i*conj(dc{n}),2i*bc{n},-Bc{n}-conj(Ac{n}),-2*jzc{n}]; 

        R2{n} = [0,conj(B2c{n}),0,2i*d2c{n}+conj(B2c{n}); ...
            B2c{n},0,-B2c{n}-2i*conj(d2c{n}),2i*b2c{n};...
            0,2i*d2c{n}-conj(B2c{n}),-conj(B2c{n});...
            B2c{n}-2i*conj(d2c{n}),2i*b2c{n},-B2c{n},0];

        zero = zeros(4);

        R{n} =  [R1{n}, R2{n}, zeroc{n}, zeroc{n}; ...
             R2{n}, R1{n}, R2{n}, zeroc{n}; ...
             zeroc{n}, R2{n}, R1{n}, R2{n}; ...
             zeroc{n}, zeroc{n}, R2{n} , R1{n}]/4;

        %For the slab we c{n}an also define an edge Raman operator
        Re1{n} = [R1{n}, zeroc{n}, zeroc{n}, zero; ...
             zeroc{n}, zeroc{n}, zeroc{n}, zeroc{n}; ...
             zeroc{n}, zeroc{n}, zeroc{n}, zeroc{n}; ...
             zeroc{n}, zeroc{n}, zeroc{n}, zeroc{n}]/4;
        Re2{n} = [zeroc{n}, zeroc{n}, zeroc{n}, zeroc{n}; ...
             zeroc{n}, zeroc{n}, zeroc{n}, zeroc{n}; ...
             zeroc{n}, zeroc{n}, zeroc{n}, zeroc{n}; ...
             zeroc{n}, zeroc{n}, zeroc{n}, R1{n}]/4;
        
    end
        
 
    
    %Finding R in diagonal space (of H), call it RR
    RR = U'*R*U;
    RR1 = U'*Re1*U;
    RR2 = U'*Re1*U;
    
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
    ent = repmat(en(1:S),1,S) + repmat(en(1:S)',S,1);
    
    %Take out the relevant chunks of the Raman operator (SC terms)
    for n=1:3
        RRr{n} = RR{n}(S+1:M,1:S);
        RRr1{n} = RR1{n}(S+1:M,1:S);
        RRr2{n} = RR2{n}(S+1:M,1:S);
    end
    
    %Check that emax is going to work
    if emax/2 < max(en)
        disp([max(en),emax,jx,jz])
    end
    %Check that the energies are positive
    if min(min(ent)) < 0
        disp([min(min(ent)),jx,jz])
    end
    
    %The six binary combinations from above.
    for m =1:6;
        W{m}(xind+N+1,yind+N+1,:,:) = pi*real(RRr{m1{m}}.*RRr{m2{m}} ...
            + RRr{m2{m}}.*RRr{m1{m}});
        W1{m}(xind+N+1,yind+N+1,:,:) = pi*real(RRr1{m1{m}}.*RRr1{m2{m}} ...
            + RRr1{m2{m}}.*RRr1{m1{m}});
        W2{m}(xind+N+1,yind+N+1,:,:) = pi*real(RRr2{m1{m}}.*RRr2{m2{m}} ...
            + RRr2{m2{m}}.*RRr2{m1{m}});
    end
            
    enty(xind+N+1,yind+N+1,:,:) = ent;              
   
    end
end

%

%histogram it
for m=1:6
    %R
    [histw, histv] = histwv(2*enty,W{m},0,emax,bins);
    I{m} = histw;
    D = histv;
    
    %Re1
    [histw, histv] = histwv(2*enty,W1{m},0,emax,bins);
    I1{m} = histw;
    %D1 = histv;
    
    %Re2
    [histw, histv] = histwv(2*enty,W2{m},0,emax,bins);
    I2{m} = histw;
    %Dpp = histv;
end

%Normalize the results by the BZ volume
Vm = sqrt(29);
Z = Vm/(2*pi)^2;

D=D*bins./(Z*emax );

I=I*bins./(Z*2*emax );
I1=I1*bins./(Z*2*emax );
I2=I2*bins./(Z*2*emax );


out = cell(1,8*3-2);

out{1} = Ev;
out{2} = D;
for m=1:6
    out{m+2} = [I{m},I1{m},I2{m}];
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