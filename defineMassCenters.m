function [MX,MY] = defineMassCenters(nodes,elements)

numElements = size(elements,1);
MX = zeros(numElements,1);
MY = zeros(numElements,1);
for k = 1:numElements
  index = elements(k,:);
  P1 = nodes(index(1),:);
  P2 = nodes(index(2),:);
  P3 = nodes(index(3),:);
  a = sqrt((P1(1)-P2(1))^2+(P1(2)-P2(2))^2);
  b = sqrt((P2(1)-P3(1))^2+(P2(2)-P3(2))^2);
  c = sqrt((P1(1)-P3(1))^2+(P1(2)-P3(2))^2);
  P = a+b+c;
  MX(k)= 1/P*(a*P3(1)+b*P1(1)+c*P2(1));
  MY(k)= 1/P*(a*P3(2)+b*P1(2)+c*P2(2));
end

end