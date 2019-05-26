function ArrFoundPeaks = run_NBNC_PSO_ES(func_num, run)
global maxFEs;
global FEs;
global initial_flag;
AlgRand = RandStream('mt19937ar','Seed',run); %set random seed
RandStream.setGlobalStream(AlgRand);
FEs = 0; 
initial_flag = 0;
dim = get_dimension(func_num); % get the dimension of function
maxFEs = get_maxfes(func_num);
lbounds = get_lb(func_num);
ubounds = get_ub(func_num);
max_v = (ubounds - lbounds)/2;
min_v = -max_v;
pop_size = get_popsize(func_num); % Set population size
%% Initialize the swarm
pSelf = rand(pop_size, dim).*(ubounds - lbounds)+lbounds;
pV = rand(pop_size, dim).*(max_v - min_v)+min_v;
pFit = fast_niching_func(pSelf,func_num);
pBest = pSelf;
pBestFit = pFit;
FEs = FEs + pop_size;

%% Evolve with NBNC_PSO_ES 
[pBest,pBestFit]=NBNC_PSO_ES(func_num, pSelf, pBest,pV,pBestFit,lbounds, ubounds, min_v, max_v);
%% The results analysis
FinalPop = pBest;
FinalFit = pBestFit;
ArrAccuracy = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5];
ArrFoundPeaks=[];
for accuracy = ArrAccuracy
    [FoundPeaks, ~] = fast_count_goptima(FinalPop, FinalFit, func_num, accuracy);
    ArrFoundPeaks = [ArrFoundPeaks, FoundPeaks];
end
end

