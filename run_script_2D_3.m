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
    out(1,:,:)= Raman2D_2(N,bins,type,flag,nn,Jx,Jz);
end
    
if nn>1
for ind = 2:nn
    
    if strcmp(type,'xx-xy')
    out(ind,:,:)= Raman2D3(N,bins,type,0,nn,Jx);
    else
    out(ind,:,:)= Raman2D_2(N,bins,type,0,nn,Jx,Jz);
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

%Plot
%position = [pixs(3)/2   20       pixs(3)/2-114   pixs(4)/2-50];
hh=figure;%('Position',position); hold on;
plot(Ev,DD);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
%title(['DOS for HyperHoneycomb Kitaev spinons: J_x=',num2str(Jx),', h=',num2str(h)])
xlabel('\omega/J^z');
ylabel('DOS');
%set(gca,'XTick',-3:3); 
%set(gca,'YTick',2*(0:5));
%hold off;
filename = ['3D_DOS_2_Jx_',num2str(round(100*Jx)),'_h_',num2str(100*h),'_28'];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);



toc