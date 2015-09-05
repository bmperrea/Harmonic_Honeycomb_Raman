function out = hH0Raman_tdl_hsl_special(N,flag,nn,TT,jx,jz,h)
%This function computes Raman spectra for the H0 lattice
%   These calculations are based off of the paper by B Perreault, K Knolle,
%   NB Perkins, and FJ Burnell, with a three spin term added (NNN Majorana
%   spinon hopping)
%   Some of the background algebra was done in a notebook by Perreault
%   called H0_finite_cell_and_BZ.nb
% We choose finite in a1.

T=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Inputs
%N - 2*N is the number of sampling points in each dimension (3 here)
%bins - the number of bins on the energy axis for the plots
%flag - whether to estimate time
%nn - the number of times this same code will be run for this time estimate
%T - the number of layers in the finite direction
%jx,jz - the spin-exchange or hopping parameters (jy=jx)
%h - the coefficient of the three-spin term (assumed isotropic)
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Unit vectors
a1 = [-1,-sqrt(2),0];
a2 = [-1,sqrt(2),0]; %%%%%%%%%%
a3 = [-1,0,3];%%%%%%%%%%

%The size of the matrices 
M = 4*T;
S = 2*T;

%Coupling choice consistent with symmetry
jy = jx;
jxa=jx; jya=jy; jxb=jx; jyb=jy;
kxa=h; kxb=h; kya=h; kyb=h; kz=h;

%the max energy for the energy axis
%emax = 4*(jx+jy+jz+6*h);  % This is the exact bandwidth

%A convenient notation for sampling points per dimension
L=2*N;

%Initialization
randy=0;%rand;
pts = (-TT):(TT-1);
%r1 = randy; r2 = randy; r3 = randy;
%pts = 1:N;

%Ev = (1:bins)'*emax/bins;

%The six are the six symmetric binary products of R^aa, R^ac, R^cc
    %h1 = {aa, aa, aa, ac, ac, cc, ab, bc};
    %h2 = {aa, ac, cc, ac, cc, cc, ab, bc};
%m1 = {1,1,1,3,3,5,2,4};
%m2 = {1,3,5,3,5,5,2,4};
%mm = size(m1,2);

%The different terms given a certain combination
%{aa,ab,ac,bc,cc}
% hxa = {1,-sqrt(2), 1,-sqrt(2),1}./4;
% hya = {1,-sqrt(2),-1, sqrt(2),1}./4;
% hxb = {1, sqrt(2), 1, sqrt(2),1}./4;
% hyb = {1, sqrt(2),-1,-sqrt(2),1}./4;
% hz =  {0,0,0,0               ,1};
% nm = 5;

%Initialize the matrices so their memory is fixed.
% RRR = cell(1,nm); RRr = RRR;
% 
% I = RRR; 
% zerolist = zeros(size(Ev));
% Dd = zerolist; Ddd=Dd;
% 
% W = cell(1,mm); 
zeroB = zeros(L,L,S,S);
enty = zeroB;
enty2 = zeros(L,L,S);
%zeromat = zeros(S,S);
% 
% for m=1:mm
%     I{m} = zerolist;  
%     W{m} = zeroB;
% end
count=0;
zeroes=0;
count1=Inf;
clockit = true;
clockit2 = true;

%High Symmetry Points:  \[CapitalGamma]-Y-T-Z-\[CapitalGamma]-X-A1-Y | T-X1 | X-A-Z   |  L-\[CapitalGamma]
HSPs = [0, 0, 0;  0, pi/sqrt(2), 0;  0, pi/sqrt(2), pi/3; ...
    0, 0, pi/3; 0, 0, 0; (29*pi)/36, 0, 0; (11*pi)/36, pi/sqrt(2), 0; ...
    0, pi/sqrt(2), 0; 0, pi/sqrt(2), pi/3; (7*pi)/36, pi/sqrt(2), pi/3; ...
    (29*pi)/36, 0, 0; (25*pi)/36, 0, pi/3; 0, 0, pi/3; ...
    pi/2, pi/(2*sqrt(2)), pi/6; 0, 0, 0];

dHSPs = diff(HSPs);

normPs = sqrt(sum(abs(dHSPs.').^2,1)).';
%onepath = sum(normPs);

%length = sum(normPs);

path = []; %onepath = [];
onepath = []; temp2 = 0;
for n = 1:size(dHSPs,1)
    tpath = [ linspace( HSPs(n,1),HSPs(n+1,1), N).', linspace( HSPs(n,2),HSPs(n+1,2), N).', linspace( HSPs(n,3),HSPs(n+1,3), N).'];
    path = vertcat(path,tpath);
    temp3 = temp2 + normPs(n);
    tonepath = linspace(temp2,temp2 + normPs(n),N).';
    onepath = vertcat(onepath,tonepath);
    temp2 = temp3;
end

tout = zeros( S^2, 2*TT);
out = zeros( 2*TT*S^2 ,size(path,1));

%Loop over the Brillouin Zone
for ind = 1:size(path,1) 
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
        time = (cl2-cl1)*((N*13*TT)/(count-count1))*nn;
        %cl = cl1 + time*(2*N-2)/(2*N);
        clockit2=false;
        
        format shortg
        disp('approximate time to take:')
        disp( datestr(time(6)/24/3600, 'DD-HH:MM:SS') )
        format
    end
    
    
    
    
    
kc = path(ind,3);%%%%%%%%%%%%%%%
kb = path(ind,2); %%%%%%%%%
ka = path(ind,1); %%%%%%%%%
            
k = [ka,kb,kc].';

%Check that this point is in the BZ, otherwise go to next point
%sn = sign(29*pi/36 - abs(ka) - abs(kb)/sqrt(2) - abs(kc)/3 );
     
if 0;%sn<0
    %disp(k)
else
%     if ylast<yind
     count = count+1;
%     ylast = yind;
%     end

%Phases

p2 = exp(1i.*(a2*k));
p3 = exp(1i.*(a3*k));

for p1ind = pts
p1 = exp(1i.*(p1ind+randy)*pi/TT);%exp(1i.*(a1*k));

%Build the matrix
DD = [jz,(jxa + jya.*p1)./p3 ;
      (jxb + jyb.*p2),jz];
alph = kyb-kya./p3-kxb*p2+kxa*p1/p3;
FF = [kz*(conj(p1)-p1), alph; -conj(alph), kz*(conj(p2)-p2)];
%F24 = [-kz*(conj(p1)-p1), -alph; conj(alph), -kz*(conj(p2)-p2)];
%zero = 0*DD;
H = 1i*[FF , DD; -DD', -FF]; 

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
 %enty(yind+N+1,zind+N+1,:,:) = ent;
 %enty2(yind+N+1,zind+N+1,:) = en2;
 
 
 tout(:,p1ind+TT+1) = ent(:);
 
end

out(:,ind) = tout(:);

end

end


%Note that Intensity should be plotted against Ev while DOS against Ev/2
%(otherwise it is the two-particle DOS

si = size(ent(:),1);


hh=figure;%('Position',position); hold on;
hold on;
% for mm = 1:si
%     plot(onepath,real(out(mm,:)));
% end

first = 1:7*N;
second = 8*N:9*N;
third = 10*N:12*N;
fourth = 13*N:14*N;

%plot(x,y);
%  pl1 = real(out(:,first.'));
%  pl2 = real(out(:,second.':));
%  pl3 =  real(out(:,third.'));
%  pl4 = real(out(:,fourth.'));

plot(onepath(first),real(out(:,first.')).');
plot(onepath(second)+1,real(out(:,second.')).');
plot(onepath(third),real(out(:,third.')).');
plot(onepath(fourth),real(out(:,fourth.')).');

%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
%title(['DOS for HyperHoneycomb Kitaev spinons: J_x=',num2str(Jx),', h=',num2str(h)])
xlabel('q');
ylabel('Spectrum');
%legend({'\rho_{--}', '\rho_{++}', '\rho_{+-}','\rho_{total}'}, 'Location', 'NorthEast');
%set(gca,'XTick',-3:3); 
%set(gca,'YTick',2*(0:5));
hold off;
filename = ['3D_Spectrum_Jx_',num2str(round(100*jx)),'_h_',num2str(100*h),'_special'];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

end