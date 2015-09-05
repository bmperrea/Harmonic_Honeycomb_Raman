tic

%type = 'xx';

%nn is number of times to average over
%N  is size of grid (2N)^3 points sampled
%N=50;
%nn=8;
%bins = 2*N;

%initialization
out = zeros(nn,bins,3);
errs = zeros(bins,3);
errm = zeros(1,3); errx = errm;

if strcmp(type,'xx-xy')
    out(1,:,:)= Raman2D3(N,bins,type,flag,nn,Jx);
else
    out(1,:,:)= Raman2D(N,bins,type,flag,nn,Jx);
end
    
if nn>1
for ind = 2:nn
    
    if strcmp(type,'xx-xy')
    out(ind,:,:)= Raman2D3(N,bins,type,0,nn,Jx);
    else
    out(ind,:,:)= Raman2D(N,bins,type,0,nn,Jx);
    end
    
end
end

%Calculate some errors
for ind=1:2
    errs(:,ind) = std(out(:,:,ind+1),0,1);
    errm(ind)   = mean(errs(:,ind));
    errx(ind)   = max (errs(:,ind));
end

disp(errm)
disp(errx)

%mean(errs)
%max(errs)
Ev = out(1,:,1);

I = mean( out(:,:,3), 1);
DE = mean( out(:,:,2), 1);


toc