II64p = cell(5,5);
slope = zeros(5);
intercept = zeros(5);

for n = 1:5;
    for m = 1:5;


out = KitaevRaman_a1_64_p_BZcuts4(30,200,1,1,1,1,0,n,m);
II64p{n,m} = out;
%figure; Ev0=1:20; plot(log(out{1}(Ev0)),log(out{2}(Ev0,1)))
pp = polyfit(log(out{1}(Ev0)),log(out{2}(Ev0,1)), 1);
slope(n,m) = pp(1);
intercept(n,m) = pp(2);

    end 
end

% I find that it is gapped very near the center of the BZ (which is good).