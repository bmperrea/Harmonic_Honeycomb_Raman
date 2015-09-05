function out = Raman_eig(N,bins,jx,jz,h) 
%This function computes Raman spectra for the Kitaev model on the
%Hyperhoneycomb lattice with a simple NNN perturbation due to a magnetic
%field (all in the zero flux sector). The algorithm diagonalizes the 
%Hamiltonian at ~10^6 k-points and computes the energies and Raman 
%vertices at those points, binning their results to their respective energy
%values, ultimately corresponding to integrals over a delta function as in
%the density of states I(w) = int_k delta(w_k - w) f_k. 

%To improve speed one might choose instead of looping to vectorize the
%steps. However, this requires an algorthm (nearly) without logical
%operations to allow vectorization, which rules out the use of the 'eig'
%function. The way to get around this would be to solve for the eigenpairs
%of the 4x4 H explicitly. This requires solving the quartic
%characteristic equation, then finding the system of 3 coupled equations
%for each eigenvalue and solving that system for the corresponding
%eigenvalue. For psuedo-random values of k we can expect that the
%eigenvalues will be non-degenerate so that gram-schmidt would not be
%required.

%Inputs


%Output


%Initial data
    %unit vectors
a1 = [0,sqrt(2),3]';
a2 = [1,0,3]';
a3 = [1,sqrt(2),0]';
    %The gamma matrix that implements the f^dagger - f symmetry
g=[0 0 1 0 ; 0 0 0 1; 1 0 0 0; 0 1 0 0];
    %The 'kappa' terms characterizing the B-field term
kx = 1; ky=kx; kz=kx;
kx = kx*h;  ky = ky*h;  kz = kz*h;
    %the max energy for the energy axis
emax = 3*max(jx,1) + .2;
    %energy bins
Ev = (1:bins)'*emax/bins;
    %The values of the Kitaev couplings
    % Jx and Jz are given
jy = jx;  %Jy is always same as Jx
    %Set the effective couplings for the Raman operators
    %number of different polarization combinations
%pols = 6;
    %The six are the six symmetric binary products of R^aa, R^ac, R^cc
    %aa = [1, 1, 0]; ac = [1,-1,0]; cc = [1/2,1/2,2];
hx = {1,1,1/2};
hy = {1,-1,1/2};
hz = {0,0,2};

%h1 = {aa, aa, aa, ac, ac, cc};
%h2 = {aa, ac, cc, ac, cc, cc};
%one = {1,1,1,1,1,1}

%psuedo-random shift of the equally spaced mesh 
rx = rand; ry = rand; rz = rand;
L=2*N;
pts = 1:L;

%Initialization
out = cell(8,1);
R=cell(3,1);
Ipp = cell(6,1);
Imm = Ipp; Ipm = Ipp;
Wpp = cell(6,1); Wmm = Wpp; Wpm = Wpp;

enp = zeros(L,L);
enm = enp;

for m =1:6
    Wpp{m} = enp;
    Wmm{m} = enp;
    Wpm{m} = enp;
end

Dpp=zeros(bins,1);
Dmm=zeros(bins,1);
Dpm=zeros(bins,1);


%The six are the six symmetric binary products of R^aa, R^ac, R^cc
    %h1 = {aa, aa, aa, ac, ac, cc};
    %h2 = {aa, ac, cc, ac, cc, cc};
m1 = {1,1,1,2,2,3};
m2 = {1,2,3,1,3,3};

    %counter for states included
Z=0;

%Loop over sites
for ix = pts
    
    
    %Timing estimate
    aleph = 1 + round( (200/N)^2 ); %iterations for a few seconds of computation
    %Or an iteration, which ever is longer
    if ix == 3 %&& flag   %the first few may be slower due to initialization
        cl1 = clock; 
    end
    if ix == 3+aleph %&& flag
        cl2 = clock;
        time = (cl2-cl1)*(2*N)/aleph;
        %cl = cl1 + time*(2*N-2)/(2*N);
        
        format shortg
        disp('approximate time to take:')
        disp( datestr(time(6)/24/3600, 'DD-HH:MM:SS') )
        format
    end
    
    for iz = pts
        for iy = 1:L
            
            x = (ix+rx-N-1)*(29*pi/36)/L; 
            y = (iy-ry)*(pi/sqrt(2))/L; 
            z = (iz+rz-N-1)*(pi/3)/L;
            
            %Is it a point in the BZ (half are and half aren't in the box
            if (29*pi/36 > abs(z)/3 + abs(y)/sqrt(2) + abs(x))
            
            Z=Z+1;
            
            r=[x y z];
            p1 = exp(-1i*r*a1); p2 = exp(-1i*r*a2); p3 = exp(-1i*r*a3);
            A = jx*p2 + jy*p1;
            B = jx + jy*p3;
            delta = ky + kx*p2 - ky*p1 - kx*p3 - kz*p3 + kz*conj(p3);
            
            alpha = [-real(A-B), -1i*imag(A+B); 1i*imag(A+B), real(A-B)];
            beta  = [1i*(imag(A-B)-2*jz), real(A+B)-1i*delta; ...
                        -real(A+B)+1i*conj(delta), -1i*(imag(A-B)+2*jz)];
            
            H = [alpha, beta; beta', -alpha]/2;
            
            %Diagonalize the Hamiltonian
            [U,d] = eig(H);
            %Sort the eigenvalues (and vectors correspondingly) in
            %ascending order
            if ~issorted(diag(d))
                [d,I] = sort(diag(d));
                U = U(:,I);
            end        
           
            %Now do the same phyisics for -k
            r=-[x y z];
            q1 = exp(-1i*r*a1); q2 = exp(-1i*r*a2); q3 = exp(-1i*r*a3);
            A = jx*q2 + jy*q1;
            B = jx + jy*q3;
            delta = ky + kx*q2 - ky*q1 - kx*q3 - kz*q3 + kz*conj(q3);
            
            alpha = [-real(A-B), -1i*imag(A+B); 1i*imag(A+B), real(A-B)];
            beta  = [1i*(imag(A-B)-2*jz), real(A+B)-1i*delta; ...
                        -real(A+B)+1i*conj(delta), -1i*(imag(A-B)+2*jz)];
            
            Hm = [alpha, beta; beta', -alpha]/2;
            
            %Diagonalize the Hamiltonian
            [Um,dm] = eig(Hm);
            %Sort the eigenvalues (and vectors correspondingly) in
            %ascending order
            if ~issorted(diag(dm))
                [dm,Im] = sort(diag(dm));
                Um = Um(:,Im);
            end             
            
            %enforce f^dagger is the dagger of f
            U(:,3) = g*conj(Um(:,1));
            U(:,4) = g*conj(Um(:,2));
            Um(:,3) = g*conj(U(:,1));
            Um(:,4) = g*conj(U(:,2));            
            %swap the order of the second half of the eigenvalues
            temp = d(3,3);
            d(3,3) = d(4,4);
            d(4,4) = temp;
            
            Ud = U'; %Conjugate transpose
            Udm= Um';
            %en = diagonal(d);
            
            %Check that it diagonalizes the matrix
            err = sum( sum( abs( Ud*H*U - d ) ) );
            if err > 10^(-5)
                disp(r)
                disp(err)
            end
            %Check that \gamma U_k \gamma = U_{-k}^*
            
            %Compute the Raman operators in the diagonal basis
                %Note that these operations on cell arrays require overriden
                %functions of each name for this data type. The benefit is
                %the very easy coding that you find below (which is makes
                %errors much less likely).
            Ac = jx*hx*p2 + jy*hy*p1;
            Bc = jx*hx + jy*hy*p3;
            %delta = ky + kx*p2 - ky*p1 - kx*p3 - kz*p3 + kz*conj(p3);
            for n=1:3
                alphac = [-real(Ac{n}-Bc{n}), -1i*imag(Ac{n}+Bc{n}); ...
                    1i*imag(Ac{n}+Bc{n}), real(Ac{n}-Bc{n})];
                betac = [1i*(imag(Ac{n}-Bc{n})-2*jz*hz{n}), ...
                    real(Ac{n}+Bc{n})-1i*delta; ...
                        -real(Ac{n}+Bc{n})+1i*conj(delta), ...
                        -1i*(imag(Ac{n}-Bc{n})+2*jz*hz{n})];
                %The three Raman operators    
                R{n} = [alphac, betac; betac', -alphac]/2;
            end
            
            %Change bases
            RR = Ud*R*U;
            %The six binary combinations from above.
            for m=1:6
                RRR = RR{m1{m}}*RR{m2{m}};
                        %Store the relevant information
                Wpp{m}(iy,iz) = pi*abs(RRR(1,3))^2;
                Wmm{m}(iy,iz) = pi*abs(RRR(2,4))^2;
                Wpm{m}(iy,iz) = pi*abs(RRR(1,4)-RRR(2,3))^2;
            end
            enp(iy,iz) = d(3,3)-d(1,1);
            enm(iy,iz) = d(4,4)-d(2,2); 
            
            %Again for -k
            RR = Udm*R*Um; %This does not work to use R here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            %The six binary combinations from above.
            for m=1:6
                RRR = RR{m1{m}}*RR{m2{m}};
                        %Store the relevant information
                Wpp{m}(iy+L,iz) = pi*abs(RRR(1,3))^2;
                Wmm{m}(iy+L,iz) = pi*abs(RRR(2,4))^2;
                Wpm{m}(iy+L,iz) = pi*abs(RRR(1,4)-RRR(2,3))^2;
            end
            enp(iy+L,iz) = d(3,3);
            enm(iy+L,iz) = d(4,4); 
            
            else
                %What to do if the point is not in the BZ
            end
        end
    end
    
    %histogram at each x so that the storage matrices don't use too much
    %memory
    for m=1:6
        [wt,vt] = histwv(2*enp(:),Wpp{m}(:),0,emax,bins);
        Ipp{m} = Ipp{m} + wt;
        Dpp = Dpp + vt;
        
        [wt,vt] = histwv(2*enm(:),Wmm{m}(:),0,emax,bins);
        Imm{m} = Imm{m} + wt;
        Dmm = Dmm + vt;        

        [wt,vt] = histwv( ( enp(:)+enp{m}(:) ) ,Wpm{m}(:),0,emax,bins);
        Ipm{m} = Ipm{w} + wt;
        Dpm = Dpm + vt;        
    end
    
end

%Normalize by number of sites counted
out{1} = Ev;
out{2} = [Dpp,Dmm,Dpm]/Z;
for m=3:8
    out{m} = [Ipp{m},Imm{m},Ipm{m}]/Z;
end


end
