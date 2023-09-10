function sMap = smoothing(sMap,h0,h1)
if ~isstruct(sMap)
    error('sMap is not a struct');
end
units = sMap.codebook;
[u1,~] = sort(h0,'descend');
[u2,~] = sort(h1,'descend');
units(h0>h1,:) = units(h0>h1,:)+0.001*(u1(1,:)-units(h0>h1,:));
units(h0<h1,:) = units(h0<h1,:)+0.001*(u2(1,:)-units(h0<h1,:));
sMap.codebook = units;
