tic

%type = 'xx';

%nn is number of times to average over
%N  is size of grid (2N)^3 points sampled
%N=50;
%nn=8;
%bins = 2*N;

%initialization
out = zeros(8,bins,7);
errs = zeros(bins,7);
errm = zeros(1,7); errx = errm;

if strcmp(type,'xx-xz') || strcmp(type,'+-') || strcmp(type,'-+')
    out(1,:,:)= Raman3D8(N,bins,type,flag,nn,Jx);
else
    out(1,:,:)= Raman3D7(N,bins,type,flag,nn,Jx);
end

if nn>1
for ind = 2:nn
    if strcmp(type,'xx-xz') || strcmp(type,'+-') || strcmp(type,'-+')
        out(ind,:,:)= Raman3D8(N,bins,type,0,nn,Jx);
    else
        out(ind,:,:)= Raman3D7(N,bins,type,0,nn,Jx);
    end
end
end

%Calculate some errors
for ind=1:6
    errs(:,ind) = std(out(:,:,ind+1),0,1);
    errm(ind)   = mean(errs(:,ind));
    errx(ind)   = max (errs(:,ind));
end
errs(:,7) = std(out(:,:,5)+out(:,:,6)+out(:,:,7));
errm(7) = mean(errs(:,7));
errx(7) = max (errs(:,7));

disp(errm)
disp(errx)

%mean(errs)
%max(errs)
Ev = out(1,:,1);

Ipp = mean( out(:,:,5), 1);
Imm = mean( out(:,:,6), 1);
Ipm = mean( out(:,:,7), 1);

%choose the placement of the plot
pixs = get(0,'screensize');
switch type
    case '++'
        position = [0           50+pixs(4)/2   pixs(3)/2   pixs(4)/2-50];
    case '--'
        position = [pixs(3)/2   50+pixs(4)/2   pixs(3)/2   pixs(4)/2-50];
    case '+-'
        position = [0           20             pixs(3)/2   pixs(4)/2-50];
end

%Plot the Raman spectrum just calculated
h=figure('Position',position);
hold on;
plot(Ev,Ipp,Ev,Imm,Ev,Ipm,Ev,Ipp+Imm+Ipm);
%errorbar(Ev,Ipp+Imm+Ipm,errs(:,4)+errs(:,5)+errs(:,6));
title(['Raman for HyperHoneycomb Kitaev for ',type,': bins=',num2str(bins),', N=',num2str(N),' Jx=',num2str(Jx)])
xlabel('E');
ylabel('I(E)');
hold off;
filename = ['3D_Raman_',type,'_10^7_Jx_',num2str(round(100*Jx)),'over100','_2'];
savefig(filename)
print(h, '-dpng', filename);

toc