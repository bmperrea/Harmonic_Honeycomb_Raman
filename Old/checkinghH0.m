kx = (xind + r1)*pi*sqrt(1/3) *1/N;%%%%%%%%%%%%%%%
ky = (yind + r2)*pi*sqrt(3/29) *1/N; %%%%%%%%%

jx=1;jz=1;
jy = jx;
jxa=jx; jya=jy; jxb=jx; jyb=jy;
%kxa=h; kxb=h; kya=h; kyb=h; kz=h;

%Phases
p1 = 0;
p2 = exp(1i.*(a2*k));
p3 = exp(1i.*(a3*k));

%Build the matrix
DD = [jz,(jya + jxa.*p1)./p3 ;
      (jxb + jyb.*p2),jz];
alph = kyb-kya./p3-kxb*p2+kxa*p1/p3;
FF = [kz*(conj(p1)-p1), alph; -conj(alph), kz*(conj(p2)-p2)];
%F24 = [-kz*(conj(p1)-p1), -alph; conj(alph), -kz*(conj(p2)-p2)];
zero = 0*DD;
H0 = 1i*[FF , DD; -DD', -FF]/2; 

D1 = [0,jxa./p3 ; 0,0];
alph1 = kxa./p3;
F1 = [-kz, alph1; 0, 0];
H1 = 1i*[F1 , D1; -D1', -F1]/2;

H=zeromat;
if T>2
for t=1:T-2
    H(4*t+1:4*t+4,4*t+1:4*t+4) = H0;
    H(4*(t-1)+1:4*(t-1)+4,4*t+1:4*t+4) = H1';
    H(4*(t+1)+1:4*(t+1)+4,4*t+1:4*t+4) = H1;
end
end
H(1:4,1:4) = H0;
if T>1
    H(4*(T-1)+1:4*(T-1)+4,4*(T-1)+1:4*(T-1)+4) = H0;
    H(5:8,1:4) = H1;
    H(4*(T-2)+1:4*(T-2)+4,4*(T-1)+1:4*(T-1)+4) = H1';
end

%Periodic B.C.'s?
if pbc == true && T>1
    H(4*(T-1)+1:4*(T-1)+4,1:4) = H1';
    H(1:4,4*(T-1)+1:4*(T-1)+4) = H1;    
end

%Diagonalize it
[V,D] = eig(H);