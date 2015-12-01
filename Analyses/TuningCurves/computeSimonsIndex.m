function SI = computeSimonsIndex(Curves)


% compute index
SI = max(Curves, [], 2)./mean(Curves, 2);