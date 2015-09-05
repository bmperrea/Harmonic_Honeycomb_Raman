N=10^5; 
numer = 1:N; 
celler = num2cell(numer); 

tic 
a=sqrt(numer); 
toc 

%The same thing in cell language is
tic 
b = cellfun(@sqrt,celler,'UniformOutput', false); 
toc