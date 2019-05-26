function [species,guide,sub_nums,meandis] = NBNC(costX,x,popsize,g,max_g)
%The NBNC algorithm is used to formulate species
alpha = 5.0;  %  scale factor for species
slbest = zeros(popsize,3);
slbest(:,1) = 1:popsize;

matdis = pdist2(x,x);
%Look for the nearest particle
matdis(logical(eye(popsize))) = inf;
[slbest(:, 3), slbest(:, 2)] = min(matdis, [], 2);

%calculate the meandis
meandis = alpha*mean(slbest(:,3));

% get the minimum distance between the current one and its better neighbor one
%find the better individual within meandis radius
matdis(matdis >= meandis) = inf;
[~, rank] = sort(costX, 'descend');
matdis = matdis(rank, rank);
matdis(triu(true(popsize, popsize))) = inf;
matdis(rank, rank) = matdis;
[slbest(:, 3), slbest(:, 2)] = min(matdis, [], 2);
slbest(slbest(:, 3) == inf,2) = slbest(slbest(:, 3) == inf,1);

sgbest = zeros(popsize,2);
sgbest(:,1) = 1:popsize;
for i = 1:popsize
    j = slbest(i,2);
    k = slbest(i,1);
    while j~= k
          k = j;
          j = slbest(k,2);
    end
    sgbest(i,2) = k;
end


%% Construct the raw species
seed_index = unique(sgbest(:,2)); 
seed_len = length(seed_index);
species  = struct();
cost_seed = zeros(1,seed_len);
for i = 1:seed_len
    species(i).seed = seed_index(i);
    species(i).idx= find(sgbest(:,2)==seed_index(i));
    species(i).num = length(species(i).idx);
    species(i).cost = costX(seed_index(i));
    species(i).subcost = costX(species(i).idx);
    cost_seed(i) = costX(seed_index(i));
end

% Sort clusters according to their fitness
%[~,index] = sort(cost_seed);
[~, index] = sort(cat(2, species.cost));
species = species(index);
%seed_x = x(seed_index(index),:);
seed_x = x(cat(2, species.seed), :);
seed_dis = pdist2(seed_x,seed_x);
seed_dis(logical(eye(seed_len))) = inf;
%Sets the flag of the cluster to be deleted

%% the meachism of merging 
mark = zeros(seed_len,2);  % mark all raw species
mark(:,1) = 1:seed_len;
mark(:,2) = mark(:,1);

for i = 1:seed_len
    [~,midx] = min(seed_dis(i,:));
    if species(i).cost < min(species(midx).subcost) 
        mark(i,2) = midx;
    end
end

for i = 1 : seed_len 
    j = mark(i,2);
    k = mark(i,1);
    while j~= k 
         k = j;
         j = mark(k,2); 
    end
    mark(i,2) = k;
end
flag = zeros(1,seed_len);
for i = 1:seed_len
    if mark(i,1)~=mark(i,2)
        flag(i) = 1;
        sgbest(species(i).idx,2) = species(mark(i,2)).seed;
        species(mark(i,2)).idx =[species(i).idx;species(mark(i,2)).idx];
        species(mark(i,2)).num = length(species(mark(i,2)).idx);
        species(mark(i,2)).subcost = costX(species(mark(i,2)).idx); 
    end   
end
species(flag == 1) = [];
sub_nums = length(species);
% sub_nums = size(species,2);
% seed = zeros(1,sub_nums);
% for i = 1 : sub_nums
%     seed(i) = species(i).seed;
% end
per = 0.5;  % the guide swiching  time, set per to be 0.5
if g > per*max_g
    guideIdx = sgbest;
else
    guideIdx = slbest(:,1:2);
end
guide = x(guideIdx(:,2),:);
end