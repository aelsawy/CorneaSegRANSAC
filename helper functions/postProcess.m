function [I,limit] = postProcess(I)

% horizontal

colMean = mean(I,2);

limit = 1;

rowMean = mean(I(limit,:));


while rowMean < mean(colMean)

    limit = limit +1;
    rowMean = mean(I(limit,:));

end

% precaution condition

if limit > size(I,1)/5 
    limit = 1;
end


I(1:limit,:) = 0;
