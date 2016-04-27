function [Centroids,refMap] = mapCentroids(Centroids, Maps, FileIndex)


[offsets, refMap] = mapFoVs(Maps);
for findex = unique(FileIndex)'
    Centroids(FileIndex==findex,:) = bsxfun(@plus, Centroids(FileIndex==findex,:), offsets(findex,[1,2]));
end