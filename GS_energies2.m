N=50;
n=3;
bins=200;

EH1 = zeros(1,n);
EH1pi = zeros(1,n);

tic
for nn = 1:n
    EH1(nn) = KitaevRaman_H1_E(N,bins,1,nn,1,1,0);
    %EH1(nn) = (1/2)*sum(H1I{1}.*H1I{3})/4*(-3/200);
    format long e
    disp(EH1(nn)); 
    format short
end
toc
tic
    format long e
disp(mean(EH1)); 
    format short
disp(std(EH1));

%round(N/2^(1/3))
for nn = 1:n
    EH1pi(nn) = KitaevRaman_H1_pi_E(round(N/2^(1/3)),bins,1,nn,1,1,0);
    %EH1pi(nn) = (1/2)*sum(H1Ipi{1}.*H1Ipi{3})/8*(-3/200);
    format long e
    disp(EH1pi(nn)); 
    format short
end

    format long e
disp(mean(EH1pi));
    format short
disp(std(EH1pi));
toc