function C = mult4(A,B)

%This function does matrix multiplication on the last two indices of a 4d
%array
if isnumeric(A) && isnumeric(B)

C = zeros(size(B));

    if isequal(size(B),size(A))
        
        for n=1:4
            for m=1:4
                for k=1:4
                    C(:,:,n,m) = C(:,:,n,m) + A(:,:,n,k).*B(:,:,k,m);
                end
            end
        end 
    
    else

        for n=1:4
            for m=1:4
                C(:,:,n) = C(:,:,n) + A(:,:,n,m).*B(:,:,m);
            end
        end         
        
    end

    %The other cases are for cells, in which case we recursively call the
    %the function for doubles (arrays)
elseif iscell(A) && iscell(B)
    C = cellfun(@mult4, A,B , 'UniformOutput', false);    
elseif iscell(B) && isnumeric(A)
    C = cellfun(@(y) mult4(A,y) ,B , 'UniformOutput', false);
elseif iscell(A) && isnumeric(B)
    C = cellfun(@(x) mult4(x,B) ,A , 'UniformOutput', false);    
else
    error('inputs need to be arrays or cells of arrays');
end

end