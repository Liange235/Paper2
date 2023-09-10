function Dism = test_som(sMap,Dte,loadings,center) 
if ~isstruct(sMap)
    error('sMap is not a struct');
end
sDiris = som_data_struct(Dte,'name','Te (test)');
bmus = som_bmus(sMap,sDiris.data,1);
for i = 1:size(Dte,1)
      sub = Dte(i,:)*loadings{bmus(i)}-center;
      Dism(i) = norm(sub);
end
