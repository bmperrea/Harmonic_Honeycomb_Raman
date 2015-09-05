function out = KitaevRaman_H1(N,bins,flag,nn,jx,jz,h)

r1 = rand; r2 = rand; r3 = rand;

pts = (-N):(N-1);

a1 = [-1,-sqrt(2),0]';
a2 = [-1,sqrt(2),0]';
a3 = [0,0,6]';

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
emax = 12 + h;

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

R = cell(1,3); 

I = cell(1,6); 
zerolist = zeros(size(Ev));
DD = zerolist; DDD=DD;

W = cell(1,6); 
zeroB = zeros(L,L,S,S);
enty = zeroB;
enty2 = zeros(L,L,S);

RRr = cell(1,3); 

for m=1:6
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
    if jx>0
        alp = jx/jz;
        jz = 3/(1+2*alp);
        jx = 3/(2+1/alp);     
        jy = jx;
    elseif jx ==0
        jz = 3; jy =0;
    else
        disp(jx);
    end
    
    jxa = jx; jxb = jx; jya = jy; jyb = jy;
    
    Ap = jya + jxa*p1;
    Ap2 = jxa + jya*p1;
    
    Bp = jyb + jxb*p2;
    Bp2 = jxb + jyb*p2;
    
      Bp = jyb + jxb*p2;
    Bm = jyb - jxb*p2;
    Bp2 = jxb + jyb*p2;
    Bm2 = jxb - jyb*p2;

    
    
    F = [jz,jxa,0,0,0,jya,0,0;0,jz,Bp,0,0,0,0,0;0,0,jz,conj(Bp2),0,0,0,0; ...
  jya.*conj(p3),0,0,jz,jxa.*p1.^2.*conj(p3),0,0,0;0,jya.*p1.^(-2),0, ...
  0,jz,jxa,0,0;0,0,0,0,0,jz,Bm,0;0,0,0,0,0,0,jz,conj(Bm2);jxa.*conj( ...
  p3),0,0,0,jya.*conj(p3),0,0,jz];
    
    zero = 0*F;

    H = [zero , F; F', zero];    
        
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
%            
%     %Raman operators (They are 3 deep in cell space)
%     jzc = jz*hz;
%     Ac = conj(p3) * (jx*hx);%+jy*p1);
%     Bc = jx*hx+jy*p2*hy;
%     dc = 0*hz;%ky*(1-conj(p3)) + kx*( -conj(p2) );% +  p1*conj(p3) );
%     ac = 0*hz;%-kz*( p1-conj(p1) );
%     bc = 0*hz;%kz*( p2 - conj(p2) );
%     %Note that a and b should be purely imaginary
%    % B2c = jy;
%    % d2c = kx;
%    % b2c = kz;
%         
%     zeroc = zero*hz; %makes a cell array of zero matrices
%     
%     for n =1:3
%     
%         R1{n} = [1i*ac{n}-1i*bc{n}+2*jzc{n},Ac{n}+conj(Bc{n}),1i*(ac{n}+bc{n}),-Ac{n}+2i*dc{n}+conj(Bc{n});...
%             Bc{n}+conj(Ac{n}),-1i*ac{n}+1i*bc{n}+2*jzc{n},-Bc{n}+conj(Ac{n})-2i*conj(dc{n}),1i*(ac{n}+bc{n});...
%             1i*(ac{n}+bc{n}),A+2i*dc{n}-conj(Bc{n}),1i*(ac{n}-bc{n}+2i*jzc{n}),-Ac{n}-conj(Bc{n});...
%             Bc{n}-conj(Ac{n})-2i*conj(dc{n}),1i*(ac{n}+bc{n}),-Bc{n}-conj(Ac{n}),1i*(-ac{n}+bc{n}+2i*jzc{n})];
% 
%     R2{n} = [0,jy*hy{n}*conj(p3),0,(-jy)*hy{n}*conj(p3);...
%             0,0,0,0;...
%             0,(jy)*hy{n}*conj(p3),0,-jy*hy{n}*conj(p3);...
%             0,-0,0,0];
% 
%         zero = zeros(4);
% 
%         %R{n} =  R1{n}/2;%[R1{n}, R2{n}, zeroc{n}, zeroc{n}; ...
%             % R2{n}', R1{n}, R2{n}, zeroc{n}; ...
%             % zeroc{n}, R2{n}', R1{n}, R2{n}; ...
%             % zeroc{n}, zeroc{n}, R2{n}' , R1{n}]/2;
%             
%    Hd = [R1{n}, R2{n}, zero, zero, zero, zero, zero, zero;...
%        R2{n}', R1{n}, R2{n}, zero, zero, zero, zero, zero; ...
%        zero, R2{n}', R1{n}, R2{n}, zero, zero, zero, zero; ...
%        zero, zero, R2{n}', R1{n}, R2{n}, zero, zero, zero; ...
%        zero, zero, zero, R2{n}', R1{n}, R2{n}, zero, zero; ...
%        zero, zero, zero, zero, R2{n}', R1{n}, R2{n}, zero; ...
%        zero, zero, zero, zero, zero, R2{n}', R1{n}, R2{n}; ...
%        zero, zero, zero, zero, zero, zero, R2{n}', R1{n}]/2;
%    
%    Hd2 = [zero, zero, zero, zero, zero, zero, zero, zero;...
%           zero, zero, zero, zero, zero, zero, zero, zero; ...
%           zero, zero, zero, zero, zero, zero, zero, zero; ...
%           zero, zero, zero, zero, zero, zero, zero, zero; ...
%           zero, zero, zero, zero, zero, zero, zero, zero; ...
%           zero, zero, zero, zero, zero, zero, zero, zero; ...
%           zero, zero, zero, zero, zero, zero, zero, zero; ...
%          R2{n}, zero, zero, zero, zero, zero, zero, zero]/2;
%         
%      R{n} = [Hd, Hd2, z2, z2, z2, z2, z2, Hd2';...
%          Hd2', Hd, Hd2, z2, z2, z2, z2, z2; ...
%          z2, Hd2', Hd, Hd2, z2, z2, z2, z2; ...
%          z2, z2, Hd2', Hd, Hd2, z2, z2, z2; ...
%          z2, z2, z2, Hd2', Hd, Hd2, z2, z2; ...
%          z2, z2, z2, z2, Hd2', Hd, Hd2, z2; ...
%          z2, z2, z2, z2, z2, Hd2', Hd, Hd2; ...
%          Hd2, z2, z2, z2, z2, z2, Hd2', Hd];
% 
%         %For the slab we can also define an edge Raman operator
%         He = [R1{n}, zero, zero, zero, zero, zero, zero, zero;...
%           zero, zero, zero, zero, zero, zero, zero, zero; ...
%           zero, zero, zero, zero, zero, zero, zero, zero; ...
%           zero, zero, zero, zero, zero, zero, zero, zero; ...
%           zero, zero, zero, zero, zero, zero, zero, zero; ...
%           zero, zero, zero, zero, zero, zero, zero, zero; ...
%           zero, zero, zero, zero, zero, zero, zero, zero; ...
%           zero, zero, zero, zero, zero, zero, zero, zero]/2;
%       
%               He2 = [zero, zero, zero, zero, zero, zero, zero, zero;...
%           zero, zero, zero, zero, zero, zero, zero, zero; ...
%           zero, zero, zero, zero, zero, zero, zero, zero; ...
%           zero, zero, zero, zero, zero, zero, zero, zero; ...
%           zero, zero, zero, zero, zero, zero, zero, zero; ...
%           zero, zero, zero, zero, zero, zero, zero, zero; ...
%           zero, zero, zero, zero, zero, zero, zero, zero; ...
%           zero, zero, zero, zero, zero, zero, zero, R1{n}]/2;
%         
%  Re{n} = [He, z2, z2, z2, z2, z2, z2, z2;...
%            z2, z2, z2, z2, z2, z2, z2, z2; ...
%            z2, z2, z2, z2, z2, z2, z2, z2; ...
%            z2, z2, z2, z2, z2, z2, z2, z2; ...
%            z2, z2, z2, z2, z2, z2, z2, z2; ...
%            z2, z2, z2, z2, z2, z2, z2, z2; ...
%            z2, z2, z2, z2, z2, z2, z2, z2; ...
%            z2, z2, z2,z2, z2, z2, z2, He2];
%         
%     end
%         
%  
%     
%     %Finding R in diagonal space (of H), call it RR
%     RR = U'*R*U;
%     RR1 = U'*Re*U;
% %    RR2 = U'*Re1*U;
%     
%     D2 = U'*H*U;
%     
% %     %Check that we diagonalized the matrix
% %     
% %     Errs = abs(D2 - D);
% %     err = max(max(max(max(Errs))));
% %     if err>10^(-6)
% %         disp(err)
% %     end
%         
%     %The relevant energies
%     en = diag(D2);
%     en2 = en(1:S);
%     
%     if emax/2 < max(en2)
%         disp([max(en),emax,jx,jz])
%     end
%     %Check that the energies are positive
%     [me,f]=min(en2);
%     if me<0
%         disp([me,f,jx,jz])
%     end
%     
%     
%     %Take out the relevant chunks of the Raman operator (SC terms)
%     for n=1:3
%         RRr{n} = RR{n}(S+1:M,1:S);
%         RRr1{n} = RR1{n}(S+1:M,1:S);
%  %       RRr2{n} = RR2{n}(S+1:M,1:S);
%     end
%     
%     
%     ent = repmat(en2,1,S) + repmat(en2',S,1);
%     
% 
%    
%     %Check that emax is going to work
%     
%     %The six binary combinations from above.
%     for m =1:6;
%         W{m}(yind+N+1,:,:) = pi*real(conj(RRr{m1{m}}).*RRr{m2{m}} ...
%             + conj(RRr{m2{m}}).*RRr{m1{m}});
%         W1{m}(yind+N+1,:,:) = pi*real(conj(RRr1{m1{m}}).*RRr1{m2{m}} ...
%             + conj(RRr1{m2{m}}).*RRr1{m1{m}});
%  %       W2{m}(xind+N+1,yind+N+1,:,:) = pi*real(RRr2{m1{m}}.*RRr2{m2{m}} ...
%  %           + RRr2{m2{m}}.*RRr2{m1{m}});
%     end
%             
%     enty(yind+N+1,:,:) = ent;           
%     
%     if h>0
%         ecut = 0.4*h*(8/sqrt(T));
%     else
%         ecut = (.015)*(8/sqrt(T));
%     end
%     
%     if min(en2)<ecut
%         q = q+1;
%         ks(q,1) = x;
%         ks(q,2) = y;

     D2 = U'*H*U;
     en = diag(D2);
     en2 = en(1:S);
     ent = repmat(en2,1,S) + repmat(en2',S,1);
     enty(yind+N+1,zind+N+1,:,:) = ent;
     enty2(yind+N+1,zind+N+1,:) = en2;
     
        else %If sn<0 then the energy is put in as -1.
                 enty(yind+N+1,zind+N+1,:,:) = -ones(S);
                enty2(yind+N+1,zind+N+1,:) = -ones(S,1);
        end
            
        end 
    end


    %histogram it
   % for m=1:6
    %R
%     if min(min(min(min(enty))))<0
%         disp(xind)
%         disp(enty)
%     end
    
    %disp(xind)
    
    [histw, histv] = histwv(2*enty(enty>0),W{m}(enty>0),0,emax,bins);
   % I{m} = I{m} + histw;
    DD = DD + histv;
 %   end
    
     [histw, histv] = histwv(enty2(enty2>0),0*enty2(enty2>0),0,emax/4,bins);
   % I{m} = I{m} + histw;
    DDD = DDD + histv;

end
%



%Normalize the results by the BZ volume
%Vm = sqrt(29)/2.913;
%Z = S^2*L^2*Vm/(2*pi)^2;

%DD=DD*bins./(Z*emax );

dd = sum(DD)*emax/bins;
DD = DD/dd;

ddd = sum(DDD)*emax/bins *1/8;
disp(sum(DDD)/(8*(2*N)^3))
DDD = DDD/ddd;

%I=I*bins./dd;%(Z*2*emax );

out = cell(1,9);

out{1} = Ev;
out{2} = DD;
out{3} = DDD;
% for m=1:6
%     out{m+3} = I{m};
% end

%Note that Intensity should be plotted against Ev while DOS against Ev/2
%(otherwise it is the two-particle DOS

end