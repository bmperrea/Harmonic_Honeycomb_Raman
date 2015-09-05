function [V,D] = eig(A)

[V,D] = cellfun(@eig,A,'UniformOutput', false);