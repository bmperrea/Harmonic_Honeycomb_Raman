%N=10;
s=1:N;
C=zeros(N,N);
h = zeros(1,N);

tic

for x = s;
 for y = s;
  for z = s;
     H = [sin(x) cos(y) tan(z) 1;
         pi x z y;
         sqrt(y+x)/10 1/z z+y x+z;
         1/x y^2/100 x*z/10 x*y/10 ];
     [V,D] = eig(H);
     B = ctranspose(V)*H*V - D;
     C(y,z) = sqrt(B(1,2)); 
  end
 end
 h=h+hist(C(:),1:N);
end

toc

tic

y=s; z=s;
for x=s;
   C = sqrt(sin(x) + cos(y) + tan(z) + sqrt(y+x)/10 + x.*z/10 + x.*y/10 )/N^3;
   h=h+hist(C(:),1:N);
end

toc