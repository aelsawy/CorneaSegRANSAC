function [best, points] = ransac(data, nparam, niter, threshDist, inlierRatio)
 
% data: a Nx2 dataset with N data points
 % nparam: the minimum number of the points of the model.
 % For line fitting problem, nparam=2
 % niter: the number of iterations
 % threshDist: the threshold of the distances between points and the fitting line
 % inlierRatio: the threshold of the number of inliers 
 
 
 N = size(data,1); % Total number of points
 
 bestInNum = 0; % Best fitting QP with largest number of inliers
 
 % parameters for best fitting QP
 best = [0 0 0];
 points = [];
 
 x = data(:,1);
 y = data(:,2);
 
 for i=1:niter
 
     % Randomly select 3 points
     % y = ax+bx+c
     % y-ax-bx-c = 0
     
     idx = randperm(N,nparam);
     
     sample = data(idx,:);   
 
     
     pol2 = polyfit(sample(:,1),sample(:,2),2);
     yestimate = polyval(pol2,x);

     % Compute the distances between all points with the fitting line      
     
     distance = abs(y-yestimate);
      
     % Compute the inliers with distances smaller than the threshold
     inlierIdx = find( distance <= threshDist);
     inlierNum = length(inlierIdx);
 
     
     % Update the number of inliers and fitting model if better model is found     
     if inlierNum >= round(inlierRatio*N) && inlierNum > bestInNum
         
         bestInNum = inlierNum;
         best = pol2;
         points = [x(inlierIdx) y(inlierIdx)];
         
     end
     
 end