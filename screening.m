function [indK] = estimating_lf(Js,dimK,vertices,faces,indana)
% Js = abs(diag(S_K));
Js = Js/max(Js);
thresh_mask = 0.9;
indices_vert = 0;
while (length(indices_vert) < dimK) && (thresh_mask > 0)
    indices_vert = find(Js > thresh_mask);
    indices_vert = intersect(indices_vert,indana);
    thresh_mask = thresh_mask-(1e-3);
end
if isempty(indices_vert)
    indK = 0;
else    
    indK = indices_vert;
    count_m = 2;
    count   = 1;
    act_th  = 6;
    d0      = 20;
    while count < count_m
        mask       = zeros(length(vertices),1);
        mask(indK) = 1;
        mask_new   = mask;
        for ii = 1:length(vertices)
            if mask(ii) == 0
                index = ii;
                [row,col] = find(faces==ii);
                findex_old = faces(row,:);
                findex_old = findex_old(:);
                findex_old = setdiff(findex_old,ii);
                vect = vertices(findex_old,:) - repmat(vertices(ii,:),[size(findex_old),1]);
                dold = sum(abs(vect).^2,2).^(1/2);
                findex_old(dold>d0) = [];
                index = [index;findex_old];
                act_near = sum(mask(index));
                if act_near > 2
                    mask_new(ii) = 1;
                end
            end
        end
        count = count+1;
        indK  = find(mask_new > 0);
    end
    count_m = 2;
    count   = 1;
    while count < count_m
        mask       = zeros(length(vertices),1);
        mask(indK) = 1;
        mask_new   = mask;
        for ii = 1:length(vertices)
            if mask(ii) == 1
                index = ii;
                [row,col] = find(faces==ii);
                findex_old = faces(row,:);
                findex_old = findex_old(:);
                findex_old = setdiff(findex_old,ii);
                vect = vertices(findex_old,:) - repmat(vertices(ii,:),[size(findex_old),1]);
                dold = sum(abs(vect).^2,2).^(1/2);
                findex_old(dold>d0) = [];
                index = [index;findex_old];
                act_near = sum(mask(index));
                if act_near <= 3
                    mask_new(ii) = 0;
                end
            end
        end
        count = count+1;
        indK  = find(mask_new > 0);
    end
end
end