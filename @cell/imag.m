function C = imag(A)

C = cellfun(@imag, A, 'UniformOutput', false);