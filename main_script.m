clc;

close all;

addpath(genpath('helper functions'))


img = imread('elsawy_od.png');

img = img(:,:,1);

layers = segment_oct_img(img);

figure, imshow(img)

hold on

for i = 1:size(layers, 1)
    
    plot(layers(i,:), '-', 'linewidth', 2)
    
end