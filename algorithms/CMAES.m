function [pop_best,bestval,bestcosts,bestIdx,countval] = CMAES(func_num, pop_best, lb, ub, bestval,stopval)
% -------------------- Initialization --------------------------------
[popsize, dim] = size(pop_best);
mu = floor(popsize/2);
pop = pop_best(1:mu, :);
xmean = mean(pop)';% objective variables initial point
if dim < 3
    sigma = 0.10;                            % coordinate wise standard deviation (step-size)
else
    sigma = 0.20;
end
% Strategy parameter setting: Selection
weights = log(mu+1/2)-log(1:mu)';       % muXone recombination weights 
mu = floor(mu);                         % number of parents/points for recombination
weights = weights/sum(weights);         % normalize recombination weights array
mueff=sum(weights)^2/sum(weights.^2);   % variance-effective size of mu

% Strategy parameter setting: Adaptation
cc = (4+mueff/dim) / (dim+4 + 2*mueff/dim);     % time constant for cumulation for C
cs = (mueff+2)/(dim+mueff+5);                   % t-const for cumulation for sigma control
c1 = 2 / ((dim+1.3)^2+mueff);                   % learning rate for rank-one update of C
cmu = 2 * (mueff-2+1/mueff) / ((dim+2)^2+2*mueff/2);    % and for rank-mu update
damps = 1 + 2*max(0, sqrt((mueff-1)/(dim+1))-1) + cs;   % damping for sigma

% Initialize dynamic (internal) strategy parameters and constants
pc = zeros(dim,1); ps = zeros(dim,1);           % evolution paths for C and sigma
B = eye(dim);                                   % B defines the coordinate system
D = eye(dim);                                   % diagonal matrix D defines the scaling
C = B*D*(B*D)';                                 % covariance matrix
chiN=dim^0.5*(1-1/(4*dim)+1/(21*dim^2));        % expectation of ||N(0,I)|| == norm(randn(N,1))

% -------------------- Generation Loop --------------------------------
countval = 0;
while countval < stopval
    % Generate and evaluate lambda offspring
    arz = randn(dim,popsize);                          % standard normally distributed vector
    for k=1:popsize
        arx(:,k) = xmean + sigma * (B*D * arz(:,k));
        countval = countval + 1;
    end
    arx = arx';
    for i = 1 : dim
        reflect = find(arx(:,i) > ub(i));
        arx(reflect,i) = ub(i) - mod((arx(reflect,i) - ub(i)), (ub(i) - lb(i)));
        reflect = find(arx(:,i) < lb(i));
        arx(reflect,i) = lb(i) + mod((lb(i) - arx(reflect,i)), (ub(i) - lb(i)));
    end
    x = arx;
    arfitness = fast_niching_func(x, func_num);
    arx = x';

    % Sort by fitness and compute weighted mean into xmean
    [arfitness, arindex] = sort(arfitness,'descend');         % maxmization
    pop_best = pop_best(arindex,:);                    %exchange pbest
    bestval = bestval(arindex);
    xmean = arx(:,arindex(1:mu))*weights;
    zmean = arz(:,arindex(1:mu))*weights;
    
    % Cumulation: Update evolution paths
    ps = (1-cs)*ps + (sqrt(cs*(2-cs)*mueff)) * (B * zmean);
    hsig = norm(ps)/sqrt(1-(1-cs)^(2*countval/popsize))/chiN < 1.4+2/(dim+1);
    pc = (1-cc)*pc + hsig * sqrt(cc*(2-cc)*mueff) * (B*D*zmean);
    
    % Adapt covariance matrix C
    C = (1-c1-cmu) * C ... 
        + c1 * (pc*pc' ... % plus rank one update
        + (1-hsig) * cc*(2-cc) * C) ... % minor correction
        + cmu ... % plus rank mu update
        * (B*D*arz(:,arindex(1:mu))) ...
        * diag(weights) * (B*D*arz(:,arindex(1:mu)))';
    
    % Adapt step-size sigma
    sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1)); 
    
    % Update B and D from C
    C = triu(C) + triu(C,1)';       % enforce symmetry
    [B,D] = eig(C);                 % eigen decomposition, B==normalized eigenvectors
    D = diag(sqrt(diag(abs(D))));   % D contains standard deviations now
%%  update bestval
     improved = arfitness > bestval;
     arx = arx';
     pop_best(improved,:) = arx(improved,:);
     bestval(improved) = arfitness(improved);
     arx = arx';
%%
     % Break, if fitness satisfies stop condition
    if abs(arfitness(1)-arfitness(popsize))<1e-7
       [bestcosts,bestIdx]  = max(bestval);
        break;
    end
end
[bestcosts,bestIdx]  = max(bestval);
end
