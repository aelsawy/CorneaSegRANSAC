function layers = segment_oct_img(I)

    
[img, limit] = topArtifact(I);

imgC = adaptiveColThresh(img);


imgCR = adaptiveRowThresh(imgC);


imgmed = medfilt2(imgCR,[3 3]);

imgmed = medfilt2(imgmed,[3 3]);

imgCR(imgmed==0) = 0;

imgCR = postProcess(imgCR);

[layers, ~, ~] = flatten_layers(imgCR,img);

layers = layers + limit;

