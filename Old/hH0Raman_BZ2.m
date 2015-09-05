function out = hH0Raman(N,bins,flag,nn,T,pbc,jx,jz,h)
%This function computes Raman spectra for the H0 lattice
%   These calculations are based off of the paper by B Perreault, K Knolle,
%   NB Perkins, and FJ Burnell, with a three spin term added (NNN Majorana
%   spinon hopping)
%   Some of the background algebra was done in a notebook by Perreault
%   called H0_finite_cell_and_BZ.nb
% We choose finite in a1.


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Inputs
%N - 2*N is the number of sampling points in each dimension (2 here)
%bins - the number of bins on the energy axis for the plots
%flag - whether to estimate time
%nn - the number of times this same code will be run for this time estimate
%T - the number of layers in the finite direction
%pbc - whether to do periodic boundary conditions
%jx,jz - the spin-exchange or hopping parameters (jy=jx)
%h - the coefficient of the three-spin term (assumed isotropic)
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Unit vectors
a2 = [sqrt(3),0]; %%%%%%%%%%
a3 = [1,-sqrt(29)]/sqrt(3);%%%%%%%%%%

%The size of the matrices 
M = 4*T;
S = 2*T;

%Coupling choice consistent with symmetry
jy = jx;
jxa=jx; jya=jy; jxb=jx; jyb=jy;
kxa=h; kxb=h; kya=h; kyb=h; kz=h;

%the max energy for the energy axis
emax = 4*(jx+jy+jz+6*h);  % This is the exact bandwidth

%A convenient notation for sampling points per dimension
L=2*N;

%Initialization
r1 = rand; r2 = rand; %r3 = rand;
pts = (-N):(N-1);

Ev = (1:bins)'*emax/bins;

%The six are the six symmetric binary products of R^aa, R^ac, R^cc
    %h1 = {aa, aa, aa, ac, ac, cc, ab, bc};
    %h2 = {aa, ac, cc, ac, cc, cc, ab, bc};
m1 = {1,1,1,3,3,5,2,4};
m2 = {1,3,5,3,5,5,2,4};
mm = size(m1,2);

%The different terms given a certain combination
%{aa,ab,ac,bc,cc}
hxa = {1,-sqrt(2), 1,-sqrt(2),1}./4;
hya = {1,-sqrt(2),-1, sqrt(2),1}./4;
hxb = {1, sqrt(2), 1, sqrt(2),1}./4;
hyb = {1, sqrt(2),-1,-sqrt(2),1}./4;
hz =  {0,0,0,0               ,1};
nm = 5;

%Initialize the matrices so their memory is fixed.
RRR = cell(1,nm); RRr = RRR;

I = RRR; 
zerolist = zeros(size(Ev));
Dd = zerolist; Ddd=Dd;

W = cell(1,mm); 
zeroB = zeros(L,S,S);
enty = zeroB;
enty2 = zeros(L,S);
zeromat = zeros(S,S);
zeromatD = zeros(T,T);

for m=1:mm
    I{m} = zerolist;  
    W{m} = zeroB;
end
count=0;
zeroes=0;
count1=Inf;
clockit = true;
clockit2 = true;

%Loop over the Brillouin Zone
for xind = pts 
    %count
   % count1
    %I display an estimate of the end time
    aleph = 1;% + round( (10/N)^2 ); %iterations for ~1 second of computation, or one iteration, which ever takes longer
    %Or an iteration, which ever is longer
    if (count >= 1) && flag && clockit   %the first few may be slower due to initialization, start from the third one
        cl1 = clock;
        count1=count;
        clockit = false;
    end
    if (count >= count1+aleph) && flag && clockit2
        cl2 = clock;
        time = (cl2-cl1)*((4*N^2)/(count-count1))*nn*.75;
        %cl = cl1 + time*(2*N-2)/(2*N);
        clockit2=false;
        
        format shortg
        disp('approximate time to take:')
        disp( datestr(time(6)/24/3600, 'DD-HH:MM:SS') )
        format
    end
    
    for yind = pts


kx = (xind + r1)*pi*sqrt(1/3) *1/N;%%%%%%%%%%%%%%%
ky = (yind + r2)*pi*sqrt(3/29) *1/N; %%%%%%%%%
            
k = [kx,ky].';

%Check that this point is in the BZ, otherwise go to next point
sn = ( abs((a2-a3)*k) < pi && abs(a2*k) < pi && abs(a3*k) < pi );
     
if sn>0
%     if ylast<yind
     count = count+1;
%     ylast = yind;
%     end

%Phases
p1 = 0;
p2 = exp(1i.*(a2*k));
p3 = exp(1i.*(a3*k));

%Build the matrix
DD = [jz,(jxa + jya.*p1)./p3 ;
      (jxb + jyb.*p2),jz];
alph = kyb-kya./p3-kxb*p2+kxa*p1/p3;
FF = [kz*(conj(p1)-p1), alph; -conj(alph), kz*(conj(p2)-p2)];
%F24 = [-kz*(conj(p1)-p1), -alph; conj(alph), -kz*(conj(p2)-p2)];


%"p1=1" "p1*=0"
D1 = [0,jya./p3 ; 0,0];
alph1 = kxa./p3;
F1 = [-kz, alph1; 0, 0];
F2 = [0,0;0,kz];

%H0 = 1i*[FF , DD; -DD', -FF]; 
%H1 = 1i*[F1 , D1; zero, F2];
%%%%%%%%%%%%%%%%%%%%%% The BZ must be wrong. The Hamilonia are identical

%From left to right goes in the positive a1 direction

DH=zeromat;
if T>1
for t=0:T-2
    DH(2*t+1:2*t+2,2*t+1:2*t+2) = DD;
    %DH(2*(t-1)+1:2*(t-1)+2,2*t+1:2*t+2) = D1;
    DH(2*t+1:2*t+2,2*(t+1)+1:2*(t+1)+2) = D1;
end
    DH(2*(T-1)+1:2*(T-1)+2,2*(T-1)+1:2*(T-1)+2) = DD;
    %DH(1:2,1:2) = D1;
    %DH(2*(T-2)+1:2*(T-2)+2,2*(T-1)+1:2*(T-1)+2) = D1;
end

if T==1
    if pbc == true 
        DH = DD + D1;
    else
        DH = DD;
    end
end

%Periodic B.C.'s?
if pbc == true && T>1
    %DH(2*(T-1)+1:2*(T-1)+2,1:2) = H1';
    DH(2*(T-1)+1:2*(T-1)+2,1:2) = D1;    
end

zero = 0*DH;

H = 1i*[zero, DH ; -DH',zero];

%Diagonalize it
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
  
%Store the energies
 D2 = U'*H*U;
 en = diag(D2);
 en2 = en(1:S);
 ent = repmat(en2,1,S) + repmat(en2',S,1);
 enty(yind+N+1,:,:) = ent;
 enty2(yind+N+1,:) = en2;
 
 
 
%Compute Raman operators
%Note that we are making frequent use of methods added
%manually to the cell data type to allow arithmetic on them
%This is convenient because the different polarization combinations are
%stored in cell space (making for fewer matrix/tensor indices).
jxa = jx; jxb = jx; jya = jy; jyb = jy;
jxac = jxa*hxa; jxbc = jxa*hxb; jyac = jya*hya; 
jybc = jya*hyb; jzc = jz*hz; 
for n =1:nm

RDD = [jzc{n},(jxac{n} + jyac{n}.*p1)./p3 ;
      (jxbc{n} + jybc{n}.*p2),jzc{n}];

%"p1=1" "p1*=0"
RD1 = [0,jyac{n}./p3 ; 0,0];

%H0 = 1i*[FF , DD; -DD', -FF]; 
%H1 = 1i*[F1 , D1; zero, F2];
%%%%%%%%%%%%%%%%%%%%%% The BZ must be wrong. The Hamilonia are identical

%From left to right goes in the positive a1 direction

RDH=zeromat;
if T>1
for t=0:T-2
    RDH(2*t+1:2*t+2,2*t+1:2*t+2) = RDD;
    %DH(2*(t-1)+1:2*(t-1)+2,2*t+1:2*t+2) = D1;
    RDH(2*t+1:2*t+2,2*(t+1)+1:2*(t+1)+2) = RD1;
end
    RDH(2*(T-1)+1:2*(T-1)+2,2*(T-1)+1:2*(T-1)+2) = RDD;
    %DH(1:2,1:2) = D1;
    %DH(2*(T-2)+1:2*(T-2)+2,2*(T-1)+1:2*(T-1)+2) = D1;
end

if T==1
    if pbc == true 
        RDH = RDD + RD1;
    else
        RDH = RDD;
    end
end

%Periodic B.C.'s?
if pbc == true && T>1
    %DH(2*(T-1)+1:2*(T-1)+2,1:2) = H1';
    RDH(2*(T-1)+1:2*(T-1)+2,1:2) = D1;    
end

R = 1i*[zero, RDH ; -RDH',zero];

RRR{n} = R;
end
RR = U'*RRR*U;

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
    W{m}(yind+N+1,:,:) = pi*real(conj(RRr{m1{m}}).*RRr{m2{m}} ...
        + conj(RRr{m2{m}}).*RRr{m1{m}});
end

else %If sn<0 then the energy is put in as -1.
     enty(yind+N+1,:,:) = -ones(S);
     enty2(yind+N+1,:) = -ones(S,1);
end
    end
    
    if sum(sum(sum(enty>0,1)))>0
    %We histogram for each kx-loop 
    [histw, histv] = histwv(2*enty(enty>0),W{1}(enty>0),0,emax,bins);
    Dd = Dd + histv;
    I{1} = I{1} + histw;
    
    if mm>1
    for m=2:mm
    [histw, ~] = histwv(2*enty(enty>0),W{m}(enty>0),0,emax,bins);
    I{m} = I{m} + histw;
    end
    end
        
    %Store the one-particle DOS as well
    [~, histv] = histwv(enty2(enty2>0),0*enty2(enty2>0),0,emax/4,bins);
    Ddd = Ddd + histv;
    
    zeroes = zeroes + sum(sum(sum(enty2==0)));
    
    end
end

%Normalize the results (takes care of infinitesimal volume element and size
%of actual BZ within the square that was sampled)
ddd = sum(Ddd)*emax/bins *1/4;
disp(sum(Ddd)/(4*(2*N)^2))
disp(zeroes/(4*(2*N)^2))
disp(count/((2*N)^2))
Ddd = Ddd./ddd;
Dd = Dd./(ddd);
I = I./(ddd);

out = cell(1,3+mm);

out{1} = Ev;
out{2} = Dd;
out{3} = Ddd;
for m=4:(3+mm)
    out{m} = I{m-3};
end

%Note that Intensity should be plotted against Ev while DOS against Ev/2
%(otherwise it is the two-particle DOS

end