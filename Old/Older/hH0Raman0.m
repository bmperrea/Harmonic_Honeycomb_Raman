function out = hH0Raman(N,bins,flag,nn,T,jx,jz,h)
%This function computes Raman spectra for the H0 lattice
%   These calculations are based off of the paper by B Perreault, j Knolle,
%   N Perkins, and Fj Burnell, with a three spin term added (NNN Majorana
%   spinon hopping)
% We choose finite in a1.


%I still need:
%1. unit vectors , BZ
%2. Hopping matrix (add 3-spin terms)
%  and Generalization of Hm to finite system of T layers 


% Unit vectors
a1 = [-1,-sqrt(2)]; %%%%%%%%%%
a2 = [-1,sqrt(2)];%%%%%%%%%%

%The size of the matrices 
M = 8*T;
S = 4*T;

%Coupling choice consistent with symmetry
jy = jx;

%the max energy for the energy axis
emax = 4*(jx+jy+jz+6*h);  % This is the exact bandwidth

%A convenient notation for sampling points per dimension
L=2*N;

%Initialization
r1 = rand; r2 = rand; r3 = rand;
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
R = cell(1,nm); RRr = R;

I = R; 
zerolist = zeros(size(Ev));
DD = zerolist; DDD=DD;

W = cell(1,mm); 
zeroB = zeros(L,S,S);
enty = zeroB;
enty2 = zeros(L,S);

for m=1:mm
    I{m} = zerolist;  
    W{m} = zeroB;
end
count=0;

%Loop over the Brillouin Zone
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


kx = (xind + r1)*pi*3/4 *1/N;%%%%%%%%%%%%%%%
ky = (yind + r2)*pi/sqrt(2) *1/N; %%%%%%%%%
            
k = [kx,ky];

%Check that this point is in the BZ, otherwise go to next point
sn = sign( 3*pi/4 - abs(x) - abs(y)/sqrt(2) );
     
if sn>0
    count = count+1;

%Phases
p1 = e^(1i.*(a1*k));
p2 = e^(1i.*(a2*k));

%Build the matrix
DD = [jz,(jya + jxa.*p1)./p3 ;
      (jxb + jyb.*p2)./p3,jz];
zero = 0*DD;
H = 1i*[zero , DD; -DD', zero]; 

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
jybc = jya*hyb; jzc = jz*hz; zeroc = 0*hz;
for n =1:nm
RD = [jzc{n},(jyac{n} + jxac{n}.*p1)./p3 ;
      (jxbc{n} + jybc{n}.*p2)./p3,jzc{n}];
R{n} = 1i*[zero , RD; -RD', zero]; 
end
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

else %If sn<0 then the energy is put in as -1.
     enty(yind+N+1,zind+N+1,:,:) = -ones(S);
     enty2(yind+N+1,zind+N+1,:) = -ones(S,1);
end
    end
    %We histogram for each kx-loop 
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

%Normalize the results (takes care of infinitesimal volume element and size
%of actual BZ within the square that was sampled)
ddd = sum(DDD)*emax/bins *1/4;
disp(sum(DDD)/(4*(2*N)^3))
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