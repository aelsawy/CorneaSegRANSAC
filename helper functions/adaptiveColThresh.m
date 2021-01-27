function I = adaptiveColThresh(I)


% horizontal projection

mu_h = mean(I,1);

sd_h = std(double(I),0,1);


mu = mean(mu_h);

sd = std(mu_h);

alpha = abs((mu_h-mu)/sd);

alpha(alpha< mean(alpha)) = mean(alpha);



for i=1:size(I,2)
    
    col=I(:,i);  
    
    col(col<mu_h(i)+alpha(i)*sd_h(i)) = 0;
    
    I(:,i) = col;
      
end
