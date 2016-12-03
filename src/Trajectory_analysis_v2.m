%% Trajectory analysis
% load('varF.mat');
num1 = length(trans_mtx);
state_track = cell(num1,3);
s_id = [];
s_init = [];
num2a = length(value2a);
for k = 1:num1
    if mod(k,100)==0
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
    
    idx = find(value2a(34:end)==trace(1),1);
    if ~isempty(idx)
        s_init = [s_init;[k,idx,trace(1)]];
    end
    
    % num2 = size(v1,1);
    state_v = zeros(num2a,2);
    % state_v = cell(num2,1);
    for j = 1:num2a;
        id1 = find(v1==value2a(j));
        if ~isempty(id1)
            state_v(j,:) = [id1(1),length(id1)];
        end
    end
    b1 = find(state_v(34:end,2)>0);
    if ~isempty(b1)
        s_id = [s_id;b1(1);length(b1)];
    end
%     for j = 1:num2
%         id1 = find(value2a==v1(j),1);
%         if ~isempty(id1)
%             state_v(j) = id1;
%             if id1>33
%                 s_id = [s_id;[k,id1,v1(j)]];
%             end
%         end
%     end
    state_track(k,:) = {v1,state_v,stay};
end

%%
num2 = length(state_track);
s_id = [];
for k = 1:num2
    if mod(k,100)==0
        fprintf('%d\n',k);
    end
    state_v = state_track{k,2};
    if isempty(state_v)
        continue;
    end
    b1 = find(state_v(34:end,2)>0);
    if ~isempty(b1)
        s_id = [s_id;[k, b1(end),length(b1)]];
    end
end

% save('s_id.mat','s_id','state_track','s_int','-v7.3');

%%
s1 = de2bi(value2,N);
tmp = state_track(s_id,:);
n1 = length(s_id);
record = cell(n1,1);
for i = 1:n1
   tmp = tmp{i,2};
   b = find(tmp>33);
   record{i} = tmp(b);
end

%%
load('s_id.mat'); % load the trajectory results


%%
id_vec = [];
load('state_value_special.mat');
[v,idx] = max(s_id(:,end));
b = find(s_id(:,end)>2); 
serial = s_id(b,1);
state_special = state_track(serial,2);
state_traj_special = state_track(serial,1);

id = 1;
n1 = length(serial);
threshold = 10;
for id = 1:n1
base = 34;
ss1 = state_special{id};
st1 = state_traj_special{id}; 
b1 = find(ss1(base:end,2)>0); %how many cell types occur
state1 = value2a(b1+base-1);
n1 = length(state1);
state_special_occur = cell(n1,1);
% if state1(1)~=2039
%     continue;
% end
for k=1:n1
   b2 = find(st1==state1(k));
   state_special_occur{k} = b2; % where the cell type occur
   if state1(k)~=2039 && b2(end)>threshold
       id_vec = [id_vec;serial(id)];
   end
%    if k>1 && b2(end)>state_special_occur{1}(1)
%        id_vec = [id_vec; serial(id)];
%    end
end
% if state1(k)==2039&&size(state1,1)>2
%    
% end

end

%%
num = length(state_track);
M = 2^11;
trans_vec = cell(num,2);
% serial1 = serial(1);
serial1 = 37985; %(21017,21201,35531,37985)
for i = serial1:serial1
   % trace = trans_mtx{i};
   if mod(i,1000)==0
       fprintf('%d\n',i);
   end
   v1 = state_track{i,1};
   if isempty(v1)
       continue;
   end
   num2 = size(v1,1);
   [state_value,ia,ic] = unique(v1,'stable');
   n1 = size(state_value,1);
   trans2 = zeros(n1);
   trans1 = [ic(1:end-1) ic(2:end)];
   idx = n1*(ic(2:end)-1)+ic(1:end-1);  % n1*(j-1)+i
   [v3,ia1,ic1] = unique(idx,'stable');
   n2 = length(v3);
   cnt = zeros(n2,1);
   for k = 1:n2;
       cnt(k) = sum(idx==v3(k));
   end
   trans2(v3) = cnt;
   trans_vec{i,1} = trans2;
   trans_vec{i,2} = state_value;
end



%% 
%match the state_value to value2a, mark the state_value with index
%in value2a if state_value(i) equals to value2a(j)

tmp = [];
for i = 1: length(state_value);
    for j = 1:43;
        if state_value(i) == value2a(j)
            tmp(i) = j;

        end
    end
end

index_update = [state_value, tmp'];
length(index_update) == length(state_value)
%save('steady_state_21017.mat','trans2','index_update')
save('trans_37985.mat','trans2')

%%
num1 = [21017,21201,35531,37985];
n1 = length(num1);
for i = 1:n1
   num1(i)
   filename1 = sprintf('steady_state_%d.mat',num1(i));
   filename2 = sprintf('transition_%d.txt',num1(i));
   load(filename1);
   dlmwrite(filename2,trans2,'delimiter','\t');
end
