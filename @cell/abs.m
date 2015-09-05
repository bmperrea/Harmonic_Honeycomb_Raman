function B = abs(A)
B = cellfun(@abs, A, 'UniformOutput', false);
end