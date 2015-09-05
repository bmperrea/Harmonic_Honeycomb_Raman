II64p5 = cell(5,5);
slope5 = zeros(5);
intercept5 = zeros(5);

for n = 1:5;
    for m = 1:5;


out = KitaevRaman_a1_64_p_BZcuts5(30,200,1,1,1,1,0,n,m);
II64p5{n,m} = out;
%figure; Ev0=1:20; plot(log(out{1}(Ev0)),log(out{2}(Ev0,1)))
pp = polyfit(log(out{1}(Ev0)),log(out{2}(Ev0,1)), 1);
slope5(n,m) = pp(1);
intercept5(n,m) = pp(2);

    end 
end

slope5

% I find that it is gapped very near the center of the BZ (which is good).