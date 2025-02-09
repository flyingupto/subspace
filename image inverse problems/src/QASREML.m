function img = QASREML(coeffs, L,types)
% Reconstruciton L level using quai-affine system
%%%%%%%%%%%%20-07-2013 by Shen ZhengWei%%%%%%%%%%%%%%%%%%%%%%%

for k=L:-1:2
  coeffs{k-1}{1,1}=QASRE(coeffs{k},1,types);
end

img=QASRE(coeffs{1},1,types);