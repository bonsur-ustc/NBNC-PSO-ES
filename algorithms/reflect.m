function [x,v] = reflect(x,xmin,xmax,popsize,DIM,v)
%lower than xmin
relectionAmount = repmat(xmin,popsize,1) - x;
oneForNeedReflectLower = relectionAmount > zeros(popsize,DIM);
relectionAmount = mod(repmat(xmin,popsize,1) - x,repmat(xmax,popsize,1)-repmat(xmin,popsize,1));
relectionAmount = (1-oneForNeedReflectLower).*zeros(popsize,DIM) + oneForNeedReflectLower.*relectionAmount;
% clampfirst
x = (1-oneForNeedReflectLower).*x + oneForNeedReflectLower.*repmat(xmin,popsize,1); 
% then reflect
x = x+ relectionAmount;
% large than xmax
relectionAmount = repmat(xmax,popsize,1) - x;
oneForNeedReflectlarger = relectionAmount < zeros(popsize,DIM);
relectionAmount = mod(x - repmat(xmax,popsize,1),repmat(xmax,popsize,1)-repmat(xmin,popsize,1));
relectionAmount = (1-oneForNeedReflectlarger).*zeros(popsize,DIM) + oneForNeedReflectlarger.*relectionAmount;
% clampfirst
x = (1-oneForNeedReflectlarger).*x + oneForNeedReflectlarger.*repmat(xmax,popsize,1); 
% then reflect
x = x - relectionAmount;
%clamp velocity
if nargin == 6 %mode 1:set velocity for reflected particles to zero
    v = (1-oneForNeedReflectLower).*v + oneForNeedReflectLower.*zeros(popsize,DIM);
    v = (1-oneForNeedReflectlarger).*v + oneForNeedReflectlarger.*zeros(popsize,DIM);
else %mode 2:set velocity to zero
    v = zeros(popsize,DIM);
end
