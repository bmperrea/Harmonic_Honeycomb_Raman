function C = uminus(A)

C = cellfun(@uminus, A, 'UniformOutput', false);