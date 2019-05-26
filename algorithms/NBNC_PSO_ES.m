function [pbest,pbestCosts] = NBNC_PSO_ES(func_num, x, pbest,v,pbestCosts,xmin, xmax, vmin, vmax)
global FEs;
global maxFEs;
lb = get_lb(func_num);
ub = get_ub(func_num);
[popsize,DIM] = size(x);
generation = 1;   % current generation
eval_cmaes = 0.2*maxFEs; % max evals for CMA-ES
maxGenPso = floor((maxFEs - eval_cmaes - FEs)/popsize); % max generations for NBNC-PSO
%% =========================== Phase1:pso evolution ==========================================================
% PSO parameters -- D. Bratton and J. Kennedy, "Defining a standard for particle swarm optimization," IEEE SIS, 2007, pp. 120?27.
% Search space parameters
w = 0.729;
c1 = 2.05*w;
c2 = 2.05*w;
str = ['The phase 1 of NBNC-PSO-ES for F',num2str(func_num),'......']; 
disp(str);
while generation < maxGenPso
    if generation < maxGenPso
        [species,guide,num,meandis] = NBNC(pbestCosts,pbest, popsize,generation,maxGenPso);
        if generation == floor(0.25*maxGenPso)
            % Strategy of balancing the species
            [species,x,v,pbest,pbestCosts,guide] = balance_species(pbestCosts,pbest,x,v,func_num,species,num,meandis,guide);
        end
        v_tmp = w*v + c1*rand(popsize,DIM).*(pbest-x) + c2*rand(popsize,DIM).*(guide-x);
        % Clamp the veloctiy 
        oneForViolation = v_tmp < repmat(vmin,popsize,1);
        v_tmp = (1-oneForViolation).*v_tmp + oneForViolation.*repmat(vmin,popsize,1); 
        oneForViolation = v_tmp > repmat(vmax,popsize,1); 
        v_tmp = (1-oneForViolation).*v_tmp + oneForViolation.*repmat(vmax,popsize,1);
        % update current location
        x_tmp = x + v_tmp;
        % reflect
        [x_tmp,v_tmp] = reflect(x_tmp,xmin,xmax,popsize,DIM,v_tmp);
        x = x_tmp;
        v = v_tmp;
        % evaluate, Update pbest,pbestCosts
        newCosts = fast_niching_func(x_tmp, func_num);
        FEs = FEs + popsize;
        improved = newCosts > pbestCosts;
        pbest(improved,:) = x(improved,:);
        pbestCosts(improved) = newCosts(improved);
        generation = generation + 1;
    end 
end
%%  =========================== Phase2:CMAES ==========================================================
str = ['The phase 2 of NBNC-PSO-ES for F',num2str(func_num),'......'];
disp(str);
stopval = maxFEs - FEs;   %max fitness evaluations
[species,num] = balancePop(pbest,pbestCosts,species,num); % merged very small species
for i = 1 : num
    index = species(i).idx;
    [~,sort_index] = sort(pbestCosts(index),'descend');
    if DIM < 3
        cut_num = 20;
    else
        cut_num = 30;
    end
    cut = min(species(i).num,cut_num);
    index = index(sort_index(1:cut));
    % Evolution with CMA-ES
    [pbest(index,:),pbestCosts(index),bestcost,bestidx,countval] = CMAES(func_num, pbest(index,:), lb, ub, pbestCosts(index),stopval);
    species(i).seed = index(bestidx);
    species(i).cost = bestcost;
    FEs = FEs + countval;
    stopval = maxFEs - FEs;
    if stopval <= 0
        break;
    end
end
end


