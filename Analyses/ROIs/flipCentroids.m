function Centroids = flipCentroids(Centroids,type)

switch type
    case 'updown'
        Centroids(:,2,:) = -Centroids(:,2,:);
%         data = Centroids(:,2,:);
%         Lims = [min(data(:)),max(data(:))];
    case 'leftright'
        Centroids(:,1,:) = -Centroids(:,1,:);
end