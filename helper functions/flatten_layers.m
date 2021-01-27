function [layers, points, distances] = flatten_layers(G, I)

    
prom = 5;
nparam = 15;
niter = 50;
thresh = 5;
frac = 0.6;


x=1:size(I,2);


sig = smoothSig(mean(G,2));

[a1, a2] = bwsequence(sig);

G(1:a1,:) = 0;

G(a2+5:end,:) = 0;



% make circular shift after fitting to be using the good estimate

firstlocs=[];

for c=1:size(G,2)
    
    pos = find(G(:,c),1,'first');
    
    if isempty(pos)        
        
        firstlocs=[firstlocs 0];
    
    else
        
        firstlocs=[firstlocs pos];
    
    end
    
end



[epi, epiPnts] = ransac([x(firstlocs~=0); firstlocs(firstlocs~=0)]', nparam, niter, thresh, frac);


if ~any(epi)
    
    epi = polyfit(x(firstlocs~=0), firstlocs(firstlocs~=0),2);
        
end

yepi = polyval(epi,x);


for c=1:size(G,2)            
        
    J(:,c) = circshift(G(:,c),-ceil(yepi(c))+50); % rounding will introduce error, I will replace it with ceil
    
end

sig = smoothSig(mean(J,2));


% remove the noise after flatting by using bwsequence

[a1, a2] = bwsequence(sig);
J(1:a1,:) = 0;
J(a2+5:end,:) = 0;

imgwidth = size(J,1);
L1 = round(imgwidth/4);
L2 = 3*round(imgwidth/4);


imgcenter = J(:,L1:L2);

muCenter = mean(imgcenter, 1);

mu = mean(muCenter);

sd = std(muCenter);

idx = muCenter< mu+2*sd;

sigCenter = smoothSig(mean(J(:,idx),2));






% search for the first and last peaks
% by inspection, first peak is the epithelium & it is the largest 
% last one is the endothelium


try

[pks,locs] = findpeaks(sigCenter,'MinPeakProminence',prom); % ,'SortStr','descend'

thelast = locs(end);


catch
    
        [~, a2] = bwsequence(sig);
        thelast = a2; % boudnary point
end


% finding the last bright points which represent the endothelium

lastlocs=[];

for c=1:size(J,2)
    
    pos = find(J(:,c),1,'last');
    
    if isempty(pos)|| pos < thelast        
        
        lastlocs = [lastlocs 0];
    
    else
        
        lastlocs = [lastlocs pos];
    
    end
    
end


% remove incomplete columns
% see if first and last distance is near the width

widths = [];

for c=1:size(J,2)
    
    colwidth = lastlocs(c)-yepi(c);
    
    widths = [widths colwidth];    
    
end

[~,order] = sort(widths,'ascend');
removed = order(1:floor(size(J,2)/5));
J(:,removed)=0;


% added to remove invalid indices
lastlocs(removed)=0;


% sig=smoothSig(mean(J,2));
sigCenter = smoothSig(mean(J(:,L1:L2),2));


% search for the first & last again


try

[pks,locs] = findpeaks(sigCenter,'MinPeakProminence',prom); % ,'SortStr','descend'

Ep = locs(1);
En = locs(end);

catch
    
    [a1, a2] = bwsequence(sig);    
    Ep = a1;
    En = a2;

end


% from here I will process the original image to have more details


for c=1:size(I,2)        
    
    flatEpi(:,c) = circshift(I(:,c),-ceil(yepi(c))+50);
    
end


sigCenter = smoothSig(mean(flatEpi(:,L1:L2),2));


% search for bowmans 1,2
% bowmans are within 150 micro (approx. 64/70 px)


sigment = sigCenter(1:Ep+70);

% search for epi again
% epi the highest peak so it is the first peak in the sorted peaks

[epipks,epilocs] = findpeaks(sigment,'SortStr','descend','MinPeakDistance',10);
Ep = epilocs(1);
EpPk = epipks(1);


%  put min distance 10 but actually it is within 21 xp (50 micro)
% update segment


sigment = sigCenter(Ep:Ep+70);
[bwpks,bwlocs] = findpeaks(sigment,'MinPeakProminence',5); % ,'SortStr','descend','MinPeakDistance',5



% I compensate the clipped part (Ep-1)

if isempty(bwlocs)
    
    Bw1 = Ep + 30;
    Bw2 = Bw1 + 10;

elseif numel(bwlocs) == 1
    
    Bw1 = bwlocs(1)+Ep-1;
    Bw2 = Bw1 + 10;

else
    
    Bw1 = bwlocs(1)+Ep-1;
    Bw2 = bwlocs(2)+Ep-1;
    
end

% search for black region between epithelium & bowmans
% I search for the local minimum between them

if ~isnan(Bw1)
    
    sig1 = sigCenter(Ep:Bw1);
    [pk,loc] = min(sig1);
    
elseif ~isnan(Bw2)
    
    sig1 = sigCenter(Ep:Bw2);
    [pk,loc] = min(sig1);
    
else
    
    sig1 = nan;
    pk = nan;
    loc = nan;
    
end




x=1:size(I,2);


% adjusting layers

yepi2 = yepi + loc;
ybo1  = yepi + Bw1-Ep;
ybo2  = yepi + Bw2-Ep;

dst_basal = loc;
dst_bw1 = Bw1-Ep;
dst_bw2 = Bw2-Ep;


% accurate estimate of decemets

idx = lastlocs~=0;
[endo, endoPnts] = ransac([x(idx); lastlocs(idx)+yepi(idx)-50]', nparam, niter, thresh, frac);


if ~any(endo) % all zeros
    
    endo = polyfit(x(idx), lastlocs(idx)+yepi(idx)-50,2);
    
end

yendo = polyval(endo,x);



for c=1:size(I,2)
        
   flatEn(:,c) = circshift(I(:,c),-round(yendo(c))+50+En-Ep);               
   
end

flatEn = flatEn(1:En-Ep+100,:);


sigCenter = mean(flatEn(:,L1:L2),2);


% search for the decemets


[pks,locs] = findpeaks(sigCenter,'MinPeakProminence',10); % ,'SortStr','descend'

En=locs(end);
ind = En-1;
Dm=En;

% make sure that the decemet's are within 64 px (150 micro)
% maybe you try to search for local minima instead of stucking into local
% minimum

while sigCenter(ind)<=sigCenter(Dm)
    
    Dm=ind;
    ind = ind -1;

end



idx = lastlocs~=0;

% endo = polyfit(x(idx),lstlocs(idx)+fstlocs(idx)-50,2);
% yendo = polyval(endo,x);

% dcm = polyfit(x(idx),lstlocs(idx)+fstlocs(idx)-(En-Dm)-50,2);

% dcm = QPransac([x(idx); lastlocs(idx)+firstlocs(idx)-(En-Dm)-50]', nparam, niter, thresh, frac);
% ydcm = polyval(dcm,x);

ydcm = yendo - (En-Dm);


dst_dm = En-Dm;
% dst_en = En-Ep;



layers = [yepi; yepi2; ybo1; ybo2; ydcm; yendo]; % -1 for using average filter 
distances = [dst_basal, dst_bw1, dst_bw2, dst_dm];

points{1} = epiPnts;
points{2} = endoPnts;
