function C = real(A)

C = cellfun(@real, A, 'UniformOutput', false);