% This is the source code of NBNC-PSO-ES.
% Author: Yingying Qiao
% Date: 2019/5/14 
%==========================main program===============================%
function ex() 
    max_run = 50;       % Total numbers of running
    test_func = 1:20;   % Function index
    mkdir(sprintf('./result/ALG'));
    PR = zeros(20, 5); %The peak ratio of algorithm 
    SR = zeros(20, 5); %The success rate of algorithm
    for func = test_func   
        if ~ismember(func, test_func) 
           continue;    
        end
        delete(gcp('nocreate'));
	    parpool('local',25); 
        spmd(25)
            result1 = run_NBNC_PSO_ES(func, labindex);% parallel for 25 workers
        end
        result = cat(1, result1{1:end});
        delete(gcp('nocreate'));
        parpool('local',25);
        spmd(25)
            result1 = run_NBNC_PSO_ES(func, labindex+25); %parallel for 25 workers
        end
        result = [result; cat(1, result1{1:end})];  % result records the number of peaks found.
        pr = mean(result) / get_no_goptima(func);  
        sr = sum(result == get_no_goptima(func)) / max_run;
        result = [pr; sr; result];
        PR(func, :) = pr;
        SR(func, :) = sr;
        dlmwrite(sprintf('./result/ALG/F%d',func), result); 
    end
    dlmwrite(sprintf('./result/ALG/BPR'), PR); 
    dlmwrite(sprintf('./result/ALG/BSR'), SR);