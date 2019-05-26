function [species,sub_nums] = balancePop(pbest,pbestCosts,species,sub_nums)
m_subsize = 2;
% seed_index = zeros(1,sub_nums);
% for i = 1 : sub_nums
%     seed_index(i) = species(i).seed;
% end
seed_index = [species.seed];
%% merge the isolated particle to the nearest species
flag = zeros(1,sub_nums);
if sub_nums > 1
   t = 1;
   while t <= sub_nums
        if species(t).num >= m_subsize
            t = t + 1;
        else
            dist_seed = pdist2(pbest(seed_index,:),pbest(seed_index,:));
            Idx = flag==1;
            dist_seed(t,Idx) = inf;
			dist_seed(t,t) = inf;
			[~,g2] = min(dist_seed(t,:)); %find the nearest species in other species    
			species(g2).idx = [species(g2).idx;species(t).idx];
			species(g2).num = species(g2).num + species(t).num;
			species(g2).subcost = pbestCosts(species(g2).idx);
			[species(g2).cost,max_index] = max(pbestCosts(species(g2).idx));
			species(g2).seed =species(g2).idx(max_index);
			seed_index(g2) = species(g2).seed;
			flag(t) = 1;  %mark the species to be merged
			t = t + 1;
         end
    end
end
species(flag==1) = []; % delete the species with very small individuals
seed_index(flag==1) = []; 
sub_nums = length(seed_index); % update the size of speice
%% sort species according to seed adaptation values
cost_seed = pbestCosts(seed_index);
[~,index] = sort(cost_seed,'descend');
species = species(index);
