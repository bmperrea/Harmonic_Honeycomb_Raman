A1 = [-1, -sqrt(2), 0]; A2 = [-1, sqrt(2), 0]; 
A3 = [-1, 0, 3]; 

Jx = 1; Jy = 1; Jz = 1;


NN=100;
p2s = (1:NN);
%zt = rand;
n=20;

p1s = (1:n);



u = 4;
Ezs = zeros(n,u);
Ezs2 = zeros(NN,n*u);
Es = zeros(NN,n*u);


for ind3 = p2s
phi3 = -pi + ind3*2*pi/NN;

for ind2 = p2s;
phi2 = -pi + ind2*2*pi/NN;
%k = [kx, ky, kz];    
    
    for ind1 = p1s;

        phi1 = 2*pi*ind1/n;
        
        Hinf = [0,0,1i.*Jz,1i.*(exp(-1i.*phi3).*Jx+exp(-1i.*(-phi1+phi3)).*Jy);...
            0,0,1i.*(Jx+exp(1i.*phi2).*Jy),1i.*Jz;...
            -1i.*Jz,1i.*((-1).*Jx+(-1).*exp(-1i.*phi2).*Jy),0,0;...
            -1i.*(exp(1i.*phi3).*Jx + exp(1i.*(-phi1+phi3)).*Jy),-1i.*Jz,0,0];

        %Store the eigenvalues
        Ezs(ind1,:) = eig(Hinf);
    end
    Ezs2(ind2,:) = sort(Ezs(:));


Hm = [0,0,0,1i.*exp(-1i.*phi3).*Jy;0,0,0,0; 0,0,0,0;0,0,0,0];

Hp = [0,0,0,0;0,0,0,0;0,0,0,0;-1i.*exp(1i.*phi3).*Jy,0,0,0];

H0 = [0,0,1i.*Jz,1i.*exp((1i*(-1)).*phi3).* ...
  Jx;0,0,1i.*(Jx+exp(1i.*phi2).*Jy),1i.* ...
  Jz;(1i*(-1)).*Jz,1i.*((-1).*Jx+(-1).*exp(((-1i)).*phi2).*Jy),0,0;(-1i).*exp(1i ...
  .*phi3).*Jx,(1i*(-1)).*Jz,0,0];

% MakeH0Layer
M = zeros(n*u);
M(1:u,1:u) = H0;
M(1:u,u*(n-1)+1:u*n) = Hm;
M(1:u,u+1:2*u) = M(1:u,u+1:2*u) + Hp;
for j=2:n-1
    M(u*(j-1)+1 : u*j, u*(j-2)+1:u*(j-1))=Hm;
    M(u*(j-1)+1 : u*j, u*(j-1)+1:u*(j))=H0;
    M(u*(j-1)+1 : u*j, u*(j)+1:u*(j+1))=Hp;
end
M(u*(n-1)+1:u*n,1:u) = Hp;
M(u*(n-1)+1:u*n,u*(n-2)+1:u*(n-1))=M(u*(n-1)+1:u*n,u*(n-1)+1:u*(n))+Hm;
M(u*(n-1)+1:u*n,u*(n-1)+1:u*(n)) = H0;

%Store the slab eigenvalues
Es(ind2,:) = sort(eig(M));
end

end

figure; plot(p2s,Ezs2);
figure; plot(p2s,Es);

test = Ezs2 - Es;
diffs = zeros(1,NN);
for j=qxs
    diffs(j) = sqrt(dot(test(j,:),test(j,:)));
end
disp(max(diffs));
