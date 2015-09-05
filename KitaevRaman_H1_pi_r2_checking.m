function out = KitaevRaman_H1_pi_r2(N,bins,flag,nn,jx,jz,~)

r1 = rand; r2 = rand; r3 = rand;

pts = (-N):(N-1);

a1 = [-1,-sqrt(2),0]';
a2 = [-1,sqrt(2),0]';
a3 = [0,0,6]';

h=0;

%In the case considered here of a larger unit cell the BZ is actually
%smaller by a factor of two. However, it is just fine to calculate within
%a doubled BZ, although we may be losing a factor of two in accuracy. The
%only barrier to the more precise calculation is a calculation of a new BZ.

%This function was written by Brent Perreault
%The basic Hamiltonian can be compared with Shaffer et al. PRL 114, 116803

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
T=2;
M = 8*T;
S = 4*T;


%the max energy for the energy axis
emax = (2*jx+jz)*4 + h;

%The values of the Kitaev couplings
jy = jx;  %Jy is always same as Jx
%Jz = 1; 
    %The 'kappa' terms characterizing the B-field term
%kx = 1; ky=kx; kz=kx;
%kx = kx*h;  ky = ky*h;  kz = kz*h;


%A convenient notation
L=2*N;

%Initialization of Arrays

Ev = (1:bins)'*emax/bins;
zeroes = zeros(size(Ev));


%The six are the six symmetric binary products of R^aa, R^ac, R^cc
    %h1 = {aa, aa, aa, ac, ac, cc, ab, bc};
    %h2 = {aa, ac, cc, ac, cc, cc, ab, bc};
m1 = {1,1,1,3,3,5,2,4};
m2 = {1,3,5,3,5,5,2,4};

%Note   I5 = -3I2 ;  I6/9 = I1 = I3/(-3)
%3/1 6/1 5/2 = -3, 9 , -3

mm = 8;

%{aa,ab,ac,bc,cc}
hxa = {1,-sqrt(2), 1,-sqrt(2),1}./4;
hya = {1,-sqrt(2),-1, sqrt(2),1}./4;
hxb = {1, sqrt(2), 1, sqrt(2),1}./4;
hyb = {1, sqrt(2),-1,-sqrt(2),1}./4;
hz =  {0,0,0,0               ,1};
nm = 5;

%Initializtion 
R = cell(1,nm); RRr = R;

I = R; 
zerolist = zeros(size(Ev));
DD = zerolist; DDD=DD;

W = cell(1,mm); 
zeroB = zeros(L,L,S,S);
enty = zeroB;
enty2 = zeros(L,L,S);

for m=1:mm
    I{m} = zerolist;  
    W{m} = zeroB;
end
count=0;
    %With matrices larger than 4, diagonalization must be done with an
    %iterative algorithm so that vectorization is not an option. Therefore,
    %we loop in this code.
for xind = pts
    
    
    %I display an estimate of the end time
    aleph = 1;% + round( (10/N)^2 ); %iterations for ~1 second of computation, or one iteration, which ever takes longer
    %Or an iteration, which ever is longer
    if (xind == -N+1) && flag   %the first few may be slower due to initialization, start from the third one
        cl1 = clock; 
        count = 0;
    end
    if (xind >= -N+1+aleph) && flag && count>0
        cl2 = clock;
        time = (cl2-cl1)*(4*N^3/count)*nn;
        %cl = cl1 + time*(2*N-2)/(2*N);
        
        format shortg
        disp('approximate time to take:')
        disp( datestr(time(6)/24/3600, 'DD-HH:MM:SS') )
        format
        
        flag = 0;
    end
    
    
    for yind = pts
        for zind = pts
            
%         x = (xind + r1)*pi*3/4 *1/N;
%         y = (yind + r2)*3*pi/sqrt(32)*1/N; 
%         z = (zind + r3)*pi/6 *1/N;

        x = (xind + r1)*9*pi/16 *1/N;
        y = (yind + r2)*7*pi/(8*sqrt(2))*1/N; 
        z = (zind + r3)*pi/6 *1/N;
    
        sn = ((4.*abs(x+(-1).*2.^(-1/2).*y))<(3.*pi))&((4.*abs(2.*x+2.^(1/2).* ...
  y))<(3.*pi))&((4.*abs((-6).*x+2.^(1/2).*y))<(19.*pi))&((4.*abs(2.* ...
  x+(-3).*2.^(1/2).*y))<(11.*pi));
             
        if sn
    count = count+1;
            
    p1 = exp(1i*(x*a1(1)+y*a1(2)+z*a1(3))); 
    p2 = exp(1i*(x*a2(1)+y*a2(2)+z*a2(3))); 
    p3 = exp(1i*(x*a3(1)+y*a3(2)+z*a3(3)));

    %switch so that jx + jy + jz = bandwidth = 3 
%     if jx>0
%         alp = jx/jz;
%         jz = 3/(1+2*alp);
%         jx = 3/(2+1/alp);     
%         jy = jx;
%     elseif jx ==0
%         jz = 3; jy =0;
%     else
%         disp(jx);
%     end

%     if jx+jy+jz ~= 3
%        disp([jx,jy,jz]) 
%     end
    
    jxa = jx; jxb = jx; jya = jy; jyb = jy;
    
    Ay = jya + jxa/p1;
    Ax = jxa + jya*p1;
    
    Aym = -jya + jxa/p1;
    Axm = -jxa + jya*p1;
   % Bp = jyb + jxb*p2;
   % Bp2 = jxb + jyb*p2;
    
    %Bp = jyb + jxb*p2;
    %Bm = jyb - jxb*p2;
    %Bp2 = jxb + jyb*p2;
    %Bm2 = jxb - jyb*p2;
    
    
    F = 1i*[jz,0,0,Ay/p3,0,0,0,0;...
        Ax,jz,0,0,0,0,0,0;...
        0,jyb,jz,0,0,jxb/p2^2,0,0; ...
  0,0,jxb,jz,0,0,jyb,0;...
  0,0,0,0,jz,0,0,Aym/p3;...
  0,0,0,0,Axm,jz,0,0;...
  0,jxb,0,0,0,jyb,jz,0;...
  0,0,jyb*p2^2,0,0,0,jxb,jz];
    
    zero = 0*F;

    Hf = [zero , F; F', zero]; 
    
    oney = eye(size(zero));
    WW = [oney, oney; -1i*oney, 1i*oney];
    H = WW'*Hf*WW/2;
        
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
    D2 = U'*H*U;
     en = diag(D2);
     en2 = en(1:S);
     ent = repmat(en2,1,S) + repmat(en2',S,1);
     enty(yind+N+1,zind+N+1,:,:) = ent;
     enty2(yind+N+1,zind+N+1,:) = en2;
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %Raman operators (They are nm deep in cell space)
     jxac = jxa*hxa; jxbc = jxa*hxb; jyac = jya*hya; jybc = jya*hyb; jzc = jz*hz;
    
    %Apc = jyac + jxac*p1;
    %Ap2c = jxac + jyac*p1;
    
    %Bpc = jybc + jxbc*p2;
    %Bp2c = jxbc + jybc*p2;
    
    Ayc = jyac + jxac./p1;
    Axc = jxac + jyac.*p1;
    
    Aymc = -jyac + jxac./p1;
    Axmc = -jxac + jyac.*p1;
    
  %  zeroc = 0*hz;
   
    for n =1:nm
        %F = 1i*[jz,0,0,Ay/p3,0,0,0,0;...
        %Ax,jz,0,0,0,0,0,0;...
        %0,jyb,jz,0,0,jxb/p2^2,0,0; ...
        %0,0,jxb,jz,0,0,jyb,0;...
        %0,0,0,0,jz,0,0,Aym/p3;...
        %0,0,0,0,Axm,jz,0,0;...
        %0,jxb,0,0,0,jyb,jz,0;...
        %0,0,jyb*p2^2,0,0,0,jxb,jz];
RF = 1i*[jzc{n},0,0,Ayc{n}./p3,0,0,0,0;...
        Axc{n},jzc{n},0,0,0,0,0,0;...
        0,jybc{n},jzc{n},0,0,jxbc{n}/p2^2,0,0; ...
  0,0,jxbc{n},jzc{n},0,0,jybc{n},0;...
  0,0,0,0,jzc{n},0,0,Aymc{n}./p3;...
  0,0,0,0,Axmc{n},jzc{n},0,0;...
  0,jxbc{n},0,0,0,jybc{n},jzc{n},0;...
  0,0,jybc{n}.*p2^2,0,0,0,jxbc{n},jzc{n}];
                           
    Rff = [zero , RF; RF', zero];   

    R{n} = WW'*Rff*WW/2;         
    end
        
    %Finding R in diagonal space (of H), call it RR
    RR = U'*R*U;
    
    if emax/2 < max(en2)
        disp([max(en2),emax,jx,jz])
    end
    %Check that the energies are positive
    [me,f]=min(en2);
    if me<0
        disp([me,f,jx,jz])
    end
        
    %Take out the relevant chunks of the Raman operator (SC terms)
    for n=1:nm
        RRr{n} = RR{n}(S+1:M,1:S);
    end
    
    %The six binary combinations from above.
    for m =1:mm;
        W{m}(yind+N+1,zind+N+1,:,:) = pi*real(conj(RRr{m1{m}}).*RRr{m2{m}} ...
            + conj(RRr{m2{m}}).*RRr{m1{m}});
    end
            
 %   enty(yind+N+1,:,:) = ent;           
    
%     if h>0
%         ecut = 0.4*h*(8/sqrt(T));
%     else
%         ecut = (.015)*(8/sqrt(T));
%     end
    
%     if min(en2)<ecut
%         q = q+1;
%         ks(q,1) = x;
%         ks(q,2) = y;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     
     
        else %If sn<0 then the energy is put in as -1.
                 enty(yind+N+1,zind+N+1,:,:) = -ones(S);
                enty2(yind+N+1,zind+N+1,:) = -ones(S,1);
        end
            
        end 
    end


   %histogram it
   % for m=1:mm
    %R
%     if min(min(min(min(enty))))<0
%         disp(xind)
%         disp(enty)
%     end
    
    %disp(xind)
    
    [histw, histv] = histwv(2*enty(enty>0),W{1}(enty>0),0,emax,bins);
    DD = DD + histv;
    I{1} = I{1} + histw;
    
    if mm>1
    for m=2:mm
    [histw, ~] = histwv(2*enty(enty>0),W{m}(enty>0),0,emax,bins);
    I{m} = I{m} + histw;
    end
    end
        
    %Store the one-particle DOS as well
    [~, histv] = histwv(enty2(enty2>0),0*enty2(enty2>0),0,emax/4,bins);
    DDD = DDD + histv;
    
    zeroes = zeroes + sum(sum(sum(enty2==0)));

end
%



%Normalize the results (takes care of infinitesimal volume element and size
%of actual BZ within the square that was sampled)
ddd = sum(DDD)*emax/bins *1/4;
disp(sum(DDD)/(8*(2*N)^3))
%disp(zeroes/(4*(2*N)^3))
DDD = DDD./ddd;
DD = DD./(ddd);
I = I./(ddd);

out = cell(1,3+mm);

out{1} = Ev;
out{2} = DD;
out{3} = DDD;
for m=4:(3+mm)
    out{m} = I{m-3};
end

%Note that Intensity should be plotted against Ev while DOS against Ev/2
%(otherwise it is the two-particle DOS

end