function whiskPoly(TifFile,ROI,summary)

data = readTiff(TifFile);

figure;
imagesc(data);

x= impoly(gca);