%% Update of the Boolean network
% Boolean Network model
% input: initial feature vector x0
%        initial continuous feature vector x0_F
%        regulatory connection func
%        indicator of whether it's a continous node or not: Flag
function [var,varF] = boolean_Run1(x0,x0_F,func,para, Flag,max_T)
N = length(x0);  % number of genes
var = zeros(max_T,N);
varF = zeros(max_T,N);  % continuous values
% para = para;      % parameters: generation rate lambda, degradation rate gamma, threshold theta
var(1,:) = x0;
varF(1,:) = x0_F;
p = ones(N,1)/N;
pro = cumsum(p);

for t = 2:max_T
    %if mod(t,10000)==0
    %    fprintf('run:%d\r\n',t);
    %end
    var(t,:) = var(t-1,:);
    varF(t,:) = varF(t-1,:);
    x = var(t-1,:);
    i = find(pro>rand(1),1);
    cons = func{i};
    num1 = length(cons);
    s = zeros(num1,1);
    flag = 0;
    cnt = 0;
    % r_func = 0;  % value of regulator function
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
        r_func = 0;         % value of regulator function
    else   % there is no valid inhibition
        r_func = cnt>0;     % combination of the effects
    end
    if Flag(i)==1  % it is a continuous5 node
        % fprintf('Node %d: %d\r\n',i,r_func);
        varF(t,i) = max(0,varF(t-1,i)+para(i,1)*r_func-para(i,2)*varF(t-1,i));
        var(t,i) = varF(t,i)>para(i,3);
    else
        var(t,i) = r_func;
    end
end



