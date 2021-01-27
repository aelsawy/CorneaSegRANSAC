function [mu, sd] = nzStats(I)


for i=1:size(I,1)
    
   row = double(I(i,:));
   row = row(row>0);
   mu(i) = mean(row);
   sd(i) = std(row);

end