function B = sum(A,n,type)
%The type should be cell 'c' or array 'a' and the function will be applied
%to the nth cell or array dimension (for c and a resp.)
if type=='a'
    B = cellfun(@(x) sum(x,n), A, 'UniformOutput', false);
elseif type == 'c' && n==1
    s = size(A);
    AA = zeros([size(A{1,1}),s(2)]);
    for j=1:(s(2))
        AA(:,:,j)=A{j};
    end
    B = sum(AA,3);
else
    error('type should be c or a for cell or array dimension n, currently it is only coded for n=1 in the cell case');
end
end