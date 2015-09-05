function answ = dxy(x, y)
% finds the derviative  d log y / d log x
%  using  = d y / d x  (x/y)

    l = length(y);
    if ~isequal(l,length(x))
        error('length of Y should be length of X ')
    end
    
    x2 = x(1:(l-1));
    y2 = y(1:(l-1));
    
    answ = diff(y).*x2 ./ ( diff(x).*y2 );
    answ = smooth(answ,5);

end

