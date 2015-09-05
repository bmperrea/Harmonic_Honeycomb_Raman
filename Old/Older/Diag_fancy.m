%Initialization

tic

M = 4;
%N = 100;

A1 = rand(M,M,N,N);
%B1 = rand(M,M,N,N);

D1 = zeros(M,M,N,N);

A = num2cell( A1 , [1,2] );


A2 = rand(N,N,M,M);
%B1 = rand(M,M,N,N);

D2 = zeros(N,N,M,M);

A3 = num2cell( A2 , [3,4] );

%B = num2cell( B1 , [1,2] );

V2=A;
V3=V2; D3=V3;

V=sparse([]);
D=sparse([]);

toc

%This one was really slow
%tic
%[V,D] = eigs( blkdiag(sparse([]),A{:}) );
%toc

tic
[V2,D2] = cellfun(@eig,A,'UniformOutput', false);
toc

tic
[V2,D2] = eig(A);
toc

tic
for i = 1:N
    for j = 1:N
        [V3{1,1,i,j},D3{1,1,i,j}] = eig(A{1,1,i,j});
    end
end
toc






