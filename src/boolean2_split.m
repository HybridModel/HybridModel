%% Boolean network
%% Regulatory connection of the network
clc; clear all;
filename = './network_th_1A.net';
[name1,inf_type,name2] = textread(filename,'%s %s %s');
type_vec = {'->','-|'};
[gene_Name,ia,ic] = unique(name2);
num1 = length(name1);
N = length(gene_Name);
func = cell(N,1);
for k = 1:N
    idx1 = find(strcmp(name2,gene_Name{k})==1);
    num2 = length(idx1);
    tmp = cell(num2,1); 
    for i = 1:num2
        str1 = strsplit(name1{idx1(i)},'&');
        idx = mapping(str1,gene_Name,1);
        if strcmp(inf_type{idx1(i)},type_vec{1})==1
            type_id = 1;
        else
            type_id = -1;
        end
        tmp{i} = [idx' type_id];
    end
    func{k} = tmp;
end

%% load the previous simulation results
% load('varF.mat');

%%
flag = zeros(N,1);
sel = [3,4,5,9];
flag(sel) = 1;
rep = 10;
a = 0.1;
max_T = 10000;
[cell_state, para, thr_v] = load_geneExp();
para(:,3) = para(:,3)./thr_v(:,4);
n2 = length(a);
var_rep = cell(n2,2);
track = cell(n2);
Flag = flag;
% sel = find(flag>0);
for a = 0.02:0.02
    k = 1;
    M = 2^4;
    for j = 1:M
    M1 = 2^7;
    xvec = zeros(M1,N);
    var_mtx = zeros(M1,1);
    varF_mtx = zeros(M1,length(sel));
    trans_vec_s = cell(M1,1);
    
    for i = 1:M1
        init = (j-1)*M1 + i-1;
        % init = i-1;
        x0 = de2bi(init,N);
        fprintf('%d, %f\r\n',i,a);
        para(sel,1) = a;
        % para(5,2) = 5*a;
        % para(3,2) = a/10;
        % para(sel,1) = para(sel,2).*para(sel,3)*1.5;
        [var, varF, trans_vec] = hybrid_RunNorm(x0,func,para,thr_v,flag,max_T,rep);
        var_mtx((i-1)*rep+1:i*rep,:) = var;
        varF_mtx((i-1)*rep+1:i*rep,:) = varF;
        trans_vec_s{i} = trans_vec;
    end
    
    var_mtx = de2bi(var_mtx,N);
    b1 = find(var_mtx(:,4)==0&sum(var_mtx,2)>0);
    sel = [1:3,5:8,10:11];
    b2 = find(var_mtx(:,4)==1&var_mtx(:,9)==1);
    b2a = find(var_mtx(:,4)==1&var_mtx(:,9)==1&sum(var_mtx(:,sel),2)==0);
    b3 = find(sum(var_mtx,2)==0);
    
    filename = sprintf('track_%.2f_%db.mat',a,j);
    var_mtx = bi2de(var_mtx);
    save(filename,'var_mtx','varF_mtx','trans_vec_s','b1','b2','b2a','b3','-v7.3');
    clear var_mtx varF_mtx trans_vec_s
    %var_rep(k,:) = {var_mtx, varF_mtx};
    %track{k} = trans_vec_s;
    end
    k = k+1;
end

%%
% var_mtxsub = de2bi(var_mtx((i-1)*rep+1:i*rep,:),N);
% 
% figure(1);
% plot(trans_vec{1,2}(:,4));
% figure(2);
% tmp1 = de2bi(trans_vec{1,1},N);
% plot(tmp1(:,9));
% ylim([0,2])
% 
% figure;
% plot(trans_vec{1,2}(:,1));
% figure;
% tmp1 = de2bi(trans_vec{1,1},N);
% plot(tmp1(:,3));
% ylim([0,2])

% %%
% figure(1);
% data1 = trans_vec_s{1}{1,2};
% plot(data1(:,3));
% hold on;
% plot(data1(:,4),'r');
% plot(data1(:,9),'k');
% hold off;
% 
% %%
% for i = 1:1
%    filename = sprintf('./track_%d.mat',i);
%    load(filename);
%     
% end 

%%
% b1 = find(varF(:,4)==0&sum(varF,2)>0);
% sel = [1:3,5:8,10:11];
% b2 = find(varF(:,4)==1&varF(:,9)==1);
% b2a = find(varF(:,4)==1&varF(:,9)==1&sum(varF(:,sel),2)==0);
% b3 = find(sum(varF,2)==0);

%%
% N = length(gene_Name);
% x0 = zeros(N,1);
% tmp = rand(N,1);
% x0(tmp>0.2) = 1;
% max_T = 100;
% var = boolean_Run(x0,func,max_T);
% var(1,:) = x0';
% p = ones(N,1)/N;
% pro = cumsum(p);
% for t = 2:max_T
%     var(t,:) = var(t-1,:);
%     x = var(t-1,:);
%     i = find(pro>rand(1),1);
%     cons = func{i};
%     num1 = length(cons);
%     s = zeros(num1,1);
%     for k = 1:num1
%         mtx = cons{k};   % constraint
%         b = sum(x(mtx(1:end-1)))==(size(mtx,2)-1); % whether the constraints are satisfied
%         s(k) = (mtx(end)==-1)+(mtx(end)*b);
%     end
%     var(t,i) = sum(s)>num1*0.5;   % combination of the effects
% end


