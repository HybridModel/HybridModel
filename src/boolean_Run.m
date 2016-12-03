%% Update of the Boolean network
% input: initial feature vector x0
%        regulatory connection func
function var = boolean_Run(x0,func,max_T)
N = length(x0);  % number of genes
var = zeros(max_T,N);
var(1,:) = x0;
p = ones(N,1)/N;
pro = cumsum(p);
for t = 2:max_T
    var(t,:) = var(t-1,:);
    x = var(t-1,:);
    i = find(pro>rand(1),1);
    cons = func{i};
    num1 = length(cons);
    s = zeros(num1,1);
    flag = 0;
    cnt = 0;
    for k = 1:num1
        mtx = cons{k};   % constraint
        b = sum(x(mtx(1:end-1)))==(size(mtx,2)-1); % whether the constraints are satisfied
        s(k) = (mtx(end)==-1)+(mtx(end)*b);
        if mtx(end)==-1&&b==1
            flag = 1;
            break;
        else if mtx(end)==1
                cnt = cnt + s(k);
            end
        end
    end
    if flag == 1   % valid inhibition
        var(t,i) = 0;
    else   % there is no valid inhibition
        var(t,i) = cnt>0;   % combination of the effects
    end
end


