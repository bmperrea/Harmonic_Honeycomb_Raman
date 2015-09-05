%Initialization

tic

M = 4;
%N = 100;

A1 = rand(M,M,N,N);
B1 = rand(M,M,N,N);

D1 = zeros(M,M,N,N);
D2 = D1;

A = num2cell( A1 , [1,2] );
B = num2cell( B1 , [1,2] );

C=A;

toc

%It is hard for me to tell which of the following approaches is best for
%large N

tic
C = cellfun(@(x,y) x*y, A,B, 'UniformOutput', false);
toc


tic

for n = 1:M
    for m = 1:M
        for k=1:M
           D1(m,n,:,:) = D1(m,n,:,:) + A1(m,k,:,:).*B1(k,n,:,:);
        end
    end
end

toc


tic

for i=1:N
    for j=1:N
        D2(:,:,i,j) = A1(:,:,i,j)*B1(:,:,i,j);
    end
end

toc


