function [popsize] = get_popsize(func_num)
% initial the size of population
popsize_all = [500*ones(1,5),2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000];
popsize = popsize_all(func_num);
end