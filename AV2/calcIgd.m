function igd = calcIgd( fit , PF )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[popSize, M]= size(fit);
 pfSize =  size(PF,1);

for i=1:pfSize
  for j=1:popSize
      D(i,j) = norm(PF(i,:)-fit(j,:),2); 
  end
  
  Dmin(i,:) = min(D(i,:));    
     
end

igd = sqrt(sum(Dmin.^2))/pfSize;

end

