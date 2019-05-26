function [species,x,v,pbest,pbestCosts,guide] = balance_species(pbestCosts,pbest,x,v,func_num,species,sub_nums,meandis,guide)
%This is source code for balance_species
dim = size(pbest,2);
MaxSpeciesSize = 20; %the max species size threshold. 
num = cat(2, species.num);
[~,num_index] = sort(num,'descend');
species = species(num_index);
%% Balance species in descending order of seed fitness 
if sub_nums > 1
   Overload_index = []; %Overload_index record the index of redundancy particles
   for i = 1 : sub_nums
       if species(i).num > MaxSpeciesSize
          [~,subindex] = sort(pbestCosts(species(i).idx),'descend');
          worst = subindex(MaxSpeciesSize+1:end);  
          over_index = species(i).idx(worst);
          Overload_index = [Overload_index;over_index];
          species(i).idx(worst) = [];
          species(i).num = MaxSpeciesSize;
          species(i).subcost = pbestCosts(species(i).idx);
       else
           if isempty(Overload_index)
               break;
           else
               if i < sub_nums
                    add_num = min(MaxSpeciesSize - species(i).num,length(Overload_index));
               else
                    add_num = length(Overload_index);
               end
               add_index =  Overload_index(1:add_num);
               Overload_index(1:add_num) = [];
               seed = species(i).seed;
               % update the infomation of particles 
               lb = repmat(-meandis,1,dim)+pbest(seed,:);  
               ub = repmat(meandis,1,dim)+pbest(seed,:);
               max_v = (ub - lb)/2;
               min_v = -max_v;
               %randomly generate particles around seed
               x(add_index,:)= rand(add_num, dim).*(ub - lb)+lb;
               v(add_index,:) = rand(add_num, dim).*(max_v - min_v)+min_v;
               lbounds = get_lb(func_num);
               ubounds = get_ub(func_num);
               % update the current location and the velocity 
               [x(add_index,:),v(add_index,:)] = reflect(x(add_index,:),lbounds,ubounds,add_num,dim); 
               pbest(add_index,:) = repmat(pbest(seed,:),add_num,1);
               pbestCosts(add_index) = repmat(pbestCosts(seed),add_num,1);
               guide(add_index,:) =  repmat(pbest(seed,:),add_num,1);
               species(i).idx = [species(i).idx;add_index];
               species(i).num = length(species(i).idx);
               species(i).subcost = [species(i).subcost;pbestCosts(add_index)];
               [species(i).cost,max_index] = max(species(i).subcost);
               species(i).seed = species(i).idx(max_index);
          end
       end
          
   end
end
end