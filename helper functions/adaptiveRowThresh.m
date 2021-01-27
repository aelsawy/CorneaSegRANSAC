function I = adaptiveRowThresh(I)


% vertical projection

[mu_v, ~] = nzStats(I);


for i=1:size(I,1)
    
    row=I(i,:); 
    
    row(row < mu_v(i)) = 0;
    
    I(i,:) = row;
      
end
