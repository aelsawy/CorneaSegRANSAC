function [I,limit] = topArtifact(I)

% vertical projection

colMean = mean(I,2);

limit = 1;

rowMean = mean(I(limit,:));

while rowMean <= mean(colMean)

    limit = limit +1;
    rowMean = mean(I(limit,:));

end

while rowMean > mean(colMean)

    limit = limit +1;
    rowMean = mean(I(limit,:));

end


% precaution condition

if limit > size(I,1)/4 
        
    mu = mean(colMean);
    sd = std(colMean);
    idx = find(colMean > (mu+2*sd));
            
    limit = idx(end);

end


I = I(limit+1:end, :);