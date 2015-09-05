%Initialization

tic

    Jx=.8; Jz = 1;
    
    x = repmat(1:N,N,1);
    z = repmat((1:N)',1,N);
    y=.2;

    A = exp(1).^(sqrt(-1).*2.^(-1/2).*(4.*x+(-4).*y+2.*2.^(1/2).*z))...
        + exp(1).^(sqrt(-1).*2.^(-1/2).*(x+(-5).*y+3.*2.^(1/2).*z));
    B = 1+exp(1).^(sqrt(-1).*2.^(-1/2).*((-1).*x+(-3).*y+(-1).*2.^(1/2).*z));
    A = A*Jx; B=B*Jx;
    
    %The Hamiltonian
    H = zeros(4,4,N,N);
    
    H(1,1,:,:) = (Jz/2); H(1,2,:,:) = 0; H(1,3,:,:) = (1/4).*(A+(-1).*B); H(1,4,:,:) = (1/4).*((-1).*A+(-1).*B);     H(2,1,:,:) = 0; H(2,2,:,:) = (-Jz/2); H(2,3,:,:) = (1/4).*(A+B); H(2,4,:,:) = (1/4).*((-1).*A+B);
    H(3,1,:,:) = (1/4).*(conj(A)+(-1).*conj(B)); H(3,2,:,:) = (1/4).*(conj(A)+conj(B)); H(3,3,:,:) = (Jz/2); H(3,4,:,:) = 0;
    H(4,1,:,:) = (1/4).*((-1).*conj(A)+(-1).*conj(B)); H(4,2,:,:) = (1/4).*((-1).*conj(A)+conj(B)); H(4,3,:,:) = 0; H(4,4,:,:) = (-Jz/2);

    H2 = num2cell( H , [1,2] );
    
    V2 = H2; D2 = H2;
   % V3 = H; D3 = A;
    V4 = H2; D4 = H2;
    
    U = H;
toc
    
%This one was really slow
%tic
%[V,D] = eigs( blkdiag(sparse([]),A{:}) );
%toc

%loop
tic
for i = 1:N
    for j = 1:N
        [V4{1,1,i,j},D4{1,1,i,j}] = eig(H2{1,1,i,j});
    end
end
toc

%cell function
tic
[V2,D2] = cellfun(@eig,H2,'UniformOutput', false);
toc

%overloaded cell function
%tic
%[V3,D3] = eig(H2);
%toc

%Explicit expression
tic

    a=A.*conj(A)+B.*conj(B);
    b=1+A.*B.*conj(A).*conj(B)+2.*real(B.*conj(A));
    c=(2+a).^2+(-4).*b;
    ep = (1/2).*2.^(-1/2).*(2+a+c.^(1/2)).^(1/2);
    em = (1/2).*2.^(-1/2).*(2+a+(-1).*c.^(1/2)).^(1/2);
    dp = c.^(1/2).*(conj(A)+conj(B))+(conj(A)+(-1).* ...
        conj(B)).*(A.*conj(A)+(-1).*B.*conj(B));
    dm = c.^(1/2).*(conj(A)+conj(B))+(-1).*(conj(A)+(-1).* ...
        conj(B)).*(A.*conj(A)+(-1).*B.*conj(B));
     
U(1,1,:,:) = c.^(1/2).*((-1/2)+em)+(1/2).*(A+(-1).*B).*(conj(A)+(-1).*conj(B));
U(1,2,:,:) = (1/2).*(B.*conj(A)+(-1).*A.*conj(B))+(-1).*em.*(A.*conj(A)+(-1).*B.*conj(B));
U(1,3,:,:) = (1/2).*(1+(-2).*em).^2.*(conj(A)+(-1).*conj(B))+(1/2).*(A+(-1).*B).*conj(A).*conj(B);
U(1,4,:,:) = (1/4).*dm;

U(2,1,:,:) = c.^(1/2).*((-1/2)+(-1).*em)+(1/2).*(A+(-1).*B).*(conj(A)+(-1).*conj(B));
U(2,2,:,:) = (1/2).*(B.*conj(A)+(-1).*A.*conj(B))+em.*(A.*conj(A)+(-1).*B.*conj(B));
U(2,3,:,:) = (1/2).*(1+2.*em).^2.*(conj(A)+(-1).*conj(B))+(1/2).*(A+(-1).*B).*conj(A).*conj(B);
U(2,4,:,:) = (1/4).*dm;

U(3,1,:,:) = c.^(1/2).*((1/2)+(-1).*ep)+(1/2).*(A+(-1).*B).*(conj(A)+(-1).*conj(B));
U(3,2,:,:) = (1/2).*(B.*conj(A)+(-1).*A.*conj(B))+(-1).*ep.*(A.*conj(A)+(-1).*B.*conj(B));
U(3,3,:,:) = (1/2).*(1+(-2).*ep).^2.*(conj(A)+(-1).*conj(B))+(1/2).*(A+(-1).*B).*conj(A).*conj(B);
U(3,4,:,:) = (-1/4).*dp;

U(4,1,:,:) = c.^(1/2).*((1/2)+ep)+(1/2).*(A+(-1).*B).*(conj(A)+(-1).*conj(B));
U(4,2,:,:) = (1/2).*(B.*conj(A)+(-1).*A.*conj(B))+ep.*(A.*conj(A)+(-1).*B.*conj(B));
U(4,3,:,:) = (1/2).*(1+2.*ep).^2.*(conj(A)+(-1).*conj(B))+(1/2).*(A+(-1).*B).*conj(A).*conj(B);
U(4,4,:,:) = (-1/4).*dp;

toc





