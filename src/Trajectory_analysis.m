%% Trajectory analysis
load('varF.mat');
num1 = length(track);
state_track = cell(num1,3);
s_id = [];
for k = 1:num1
    if mod(k,1000)==0
        fprintf('%d\n',k);
    end
    trace = trans_mtx{k};
    change = trace(2:end,:)-trace(1:end-1,:);
    b = find(change~=0);
    if isempty(b)
        continue;
    end
    
    v1 = [trace(1);trace(b+1)];
    stay = [b(1);b(2:end)-b(1:end-1);size(trace,1)-b(end)];
    
    num2 = size(v1,1);
    state_v = zeros(num2,1);
    for j = 1:num2
        id1 = find(value2a==v1(j),1);
        if ~isempty(id1)
            state_v(j) = id1;
            if id1>34
                s_id = [s_id;k];
            end
        end
    end
    state_track(k,:) = {v1,state_v,stay};
end

%%
tmp1 = state_track(s_id,:);
n1 = length(s_id);
record = cell(n1,1);
for i = 1:n1
   tmp2 = tmp1{i,2};
   b = find(tmp2>33);
   record{i} = tmp2(b);
end

%%
s1 = de2bi(value2,N);
