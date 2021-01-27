function [l1,l2] = bwsequence(sequence)

idx = find(sequence);

l1 = idx(1);
l2 = l1;

L = [];

for i=2:numel(idx)

    if idx(i) == idx(i-1)+1
        l2 = idx(i);
        
    else
        L = [L; l1 l2];
        l1 = idx(i);
        l2 = l1;
    end

end

L = [L; l1 l2];

D = L(:,2)-L(:,1);
[~,id] = max(D);

l1 = L(id,1);
l2 = L(id,2);
