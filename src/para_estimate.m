%% parameter estimate
load('varF.mat');

%% gene expression
gene_exp1 = [2.8	3.02	3	2.97	2.99	8.15	11.76	12.62	4.11	2.78
7.47	5.65	5.59	5.61	5.63	8.71	7.98	7.28	7.98	6.68
14.25	10.56	13.05	11.16	13.45	10.54	13.12	12.83	13.12	14.04
2.16	2.65	3.16	2.57	2.97	11.56	5.26	2.57	4.12	2.07
3.02	3.27	3.15	3.13	3.1	3.92	4.78	9.69	3.27	3.05
12.48	2.78	7.87	2.58	6.49	7.48	12.8	10.98	12.05	11.69
6.62	4.99	4.91	4.52	5.53	4.59	10.03	6.36	10.21	4.63
11.23	8.63	11.39	5.06	9.93	10.2	9.92	9.55	9.92	11.19
2.49	2.6	2.72	2.49	2.5	13.69	7.78	11.66	5.34	2.49
1.89	1.94	1.92	1.9	1.9	1.91	1.91	1.9	1.92	1.9
6.47	7.93	6.45	8.1	6.77	12.43	5.92	6.32	3.31	7.68]';

%% Cell type states
load('varF.mat');
cell_state = zeros(10,N);
ery_id = [1,2,4,6,8:9,11];
cd8a_id = [3,11]; cd4a_id = [3,8,11];
mono_id = [2:3,6:8];
bcell_id = [2:3,6:8,11];
cd8_id = [3,6:8,11];
cd4_id = [3,6,8,11];
nk_id = [2:3,6,8,11];
granu_id = [1,3,6:9,11];
hspc_id = [1:3,5:11];
cell_id = cell(10,1);
cell_id{1} = ery_id; cell_id{2} = cd8a_id;
cell_id{3} = cd4a_id; cell_id{4} = mono_id;
cell_id{5} = bcell_id; cell_id{6} = cd8_id;
cell_id{7} = cd4_id; cell_id{8} = nk_id;
cell_id{9} = granu_id; cell_id{10} = hspc_id;
for i = 1:10
   cell_state(i,cell_id{i}) = 1; 
end

%% Update of the network
% N = length(gene_Name);
% max_T = 1000;
% M = 2^11;
% rep = 20;
% M1 = M*rep;
% xvec = zeros(M1,N);
% varF = zeros(size(xvec));
% trans_vec = cell(M1,1);
% for i = 1:M
%     xvec((i-1)*rep+1:i*rep,:) = ones(rep,1)*de2bi(i-1,N);
% end
% for i = 1:M1
%     if mod(i,100)==0 
%         fprintf('Run: %d\n',i);
%     end
%     x0 = xvec(i,:);
%     var = boolean_Run(x0,func,max_T);
%     varF(i,:) = var(end,:);
%     trans_vec{i} = var;
% end

%%
n1 = length(b1);
start = zeros(n1,2);
for i = 1:n1
    idx = b1(i);
    trace = de2bi(trans_mtx{idx},N);
    trace1 = trace(:,9);
    b = find(trace1>0,1);
    start(i,:) = [b,find(trace1(b:end)<1,1)+b-1];
end

%%
b = find(start(:,2)~=start(:,1));
diff_idx = b;
id1 = b1(b);
sub_trace1 = start(b,:);
n1 = length(b);
figure;
for k = 1:10
id = ceil(n1*rand(1));
subplot(5,2,k);
trace = de2bi(trans_mtx{id1(k)},N);
plot(trace(:,9));
hold on;
plot(trace(:,5),'r');
plot(trace(:,3),'k');
ylim([0,2]);
hold off;
end

%%
cell_type = {'Erythroid','CD8-activated','CD4-activated','Monocyte','B-Cell',...
             'CD8','CD4','NK','Granulocyte','HSPC'};
cell_type1 = {'B-Cell','CD4-activated','CD4','CD8-activated','CD8',...
              'Erythroid','Granulocyte','HSPC','Monocyte','NK'};
n1 = length(cell_type1);
sel = zeros(n1,1);
for i = 1:n1
    b = find(strcmp(cell_type,cell_type1{i})==1,1);    
    sel(i) = b;
end
gene_exp(sel,:) = gene_exp1;

%%
para = zeros(N,3);
flag = zeros(N,1);
thr_v = zeros(N,5);
for id2 = 1:N
s1 = cell_state(:,id2);
c1 = find(s1==0);
c2 = find(s1==1);
g1 = gene_exp(c1,id2);
g2 = gene_exp(c2,id2);
thr1 = min(g1); thr1a = max(g1);
thr2 = max(g2); thr2a = min(g2);
thresh = (max(g1)+min(g2))*0.5;
para(id2,3) = thresh;
thr_v(id2,:) = [thr1,thr1a,thr2a,thr2,thresh];
end

%%
id1 = b1(diff_idx);
dis = sub_trace1(:,2)-sub_trace1(:,1); % time from off to on
sel = [3,9];
flag(sel)=1;
para(sel,2) = 1;
para(sel,1) = para(sel,2).*para(sel,3)*1.5;






