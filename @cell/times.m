function C = times(A,B)

if iscell(A) && iscell(B)
    C = cellfun(@(x,y) x.*y, A,B, 'UniformOutput', false);
elseif iscell(A)
    C = cellfun(@(x) x.*B, A, 'UniformOutput', false);
else %B is a cell array
    C = cellfun(@(y) A.*y, B, 'UniformOutput', false);
end