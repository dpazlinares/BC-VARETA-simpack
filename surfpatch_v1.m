function [newindex]=surfpatch_v1(Source,vertices,faces,d0,neibormax)

index      = Source;
[row,col]  = find(faces==Source);
findex_old = faces(row,:);
findex_old = findex_old(:);
findex_old = setdiff(findex_old,Source,'stable');
vect       = vertices(findex_old,:) - repmat(vertices(Source,:),[size(findex_old),1]);
dold       = sum(abs(vect).^2,2).^(1/2);
findex_old(dold>d0) = [];
index         = [index;findex_old];


while length(index) < neibormax
    findex = [];
    for i=1:length(findex_old)
        [row,col] = find(faces==findex_old(i));
        index1    = faces(row,:);
        index1    = index1(:);
        index1    = setdiff(index1,index,'stable');
        vect      = vertices(index1,:) - repmat(vertices(findex_old(i),:),[size(index1),1]);
        dold      = sum(abs(vect).^2,2).^(1/2);
        index1(dold>d0) = [];
        findex    = [findex;index1];
    end
    findex_old = findex;
    index      = [index;findex];
    index      = unique(index,'stable');
end

newindex = index;

end