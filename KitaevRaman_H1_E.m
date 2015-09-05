function energy = KitaevRaman_H1_E(N,~,flag,nn,jx,jz,~)

r1 = rand; r2 = rand; r3 = rand;

pts = (-N):(N-1);

a1 = [-1,-sqrt(2),0]';
a2 = [-1,sqrt(2),0]';
a3 = [0,0,6]';



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
T=1;
M = 8*T;
S = 4*T;

h=0;
%the max energy for the energy axis
emax = 12 + h;

%The values of the Kitaev couplings
jy = jx;  %Jy is always same as Jx
%Jz = 1; 
    %The 'kappa' terms characterizing the B-field term
% kx = 1; ky=kx; kz=kx;
% kx = kx*h;  ky = ky*h;  kz = kz*h;


%A convenient notation
L=2*N;

count=0;
energy=0;
    %With matrices larger than 4, diagonalization must be done with an
    %iterative algorithm so that vectorization is not an option. Therefore,
    %we loop in this code.
for xind = pts
    
    
    %I display an estimate of the end time
    aleph = 1;% + round( (10/N)^2 ); %iterations for ~1 second of computation, or one iteration, which ever takes longer
    %Or an iteration, which ever is longer
    if (xind == -N+1) && flag   %the first few may be slower due to initialization, start from the third one
        cl1 = clock;
        count=0;
    end
    if (xind == -N+1+aleph) && flag
        cl2 = clock;
        time = (cl2-cl1)*(4*N^3/count)*nn;
        %cl = cl1 + time*(2*N-2)/(2*N);
        
        format shortg
        disp('approximate time to take:')
        disp( datestr(time(6)/24/3600, 'DD-HH:MM:SS') )
        format
    end
    
    
    for yind = pts
        for zind = pts
            
        x = (xind + r1)*pi*3/4 *1/N;
        y = (yind + r2)*pi/sqrt(2) *1/N; 
        z = (zind + r3)*pi/6 *1/N;
    
        sn = sign( 3*pi/4 - abs(x) - abs(y)/sqrt(2) );
     
        if sn>0
          %  count = count+1;
        
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
    
    jxa = jx; jxb = jx; jya = jy; jyb = jy;
    
    Ap = jya + jxa*p1;
    Ap2 = jxa + jya*p1;
    
    Bp = jyb + jxb*p2;
    Bp2 = jxb + jyb*p2;
    
    F = [jz,conj(Ap2),0,0;0,jz,Bp,0;0,0,jz,conj(Bp2);Ap.*conj(p3),0,0,jz];
    
    zero = 0*F;

    H = [zero , F; F', zero];    
        
    P = eig(H);
    
    energy=energy+sum(abs(P));
    count=count+1;
   
        end
        
        end 
    end
end


energy = energy/(2*M*count);


end