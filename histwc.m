function [histw, vint] = histwc(vv, ww, minV, maxV, bins)
    delta = (maxV-minV)/bins;
    vint = linspace(minV,maxV,bins)-delta/2.0;
    histw = zeros(bins,1);
    for i = 1:length(vv)
        ind = find(vint<vv(i),1,'last');
        if ~isempty(ind)
            histw(ind) = histw(ind) + ww(i);
        end
    end
end