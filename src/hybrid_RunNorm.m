%% boolean network: continuous
function [var, varF, trans_vec] = hybrid_RunNorm(x0,func,para,thr_v,Flag,max_T,rep)
N = length(Flag);
M1 = rep;
xvec = zeros(M1,N);
flag = Flag;
sel = find(flag>0);
var = zeros(M1,1);
varF = zeros(M1,length(sel));
trans_vec = cell(M1,2);

xvec(1:rep,:) = ones(rep,1)*x0;
for i = 1:rep
    % fprintf('Replicate: %d\r\n',i);
    x0 = xvec(i,:);
    x0_F = x0;
    tmp = zeros(length(sel),1);
    for k = 1:length(sel)
      if x0(sel(k))==1
          tmp(k) = 1;
      else
          tmp(k) = 0;
      end
    end
    dis = tmp - para(sel,3);
    x0_F(sel) = tmp - dis*0.2.*rand(length(sel),1);
    
    [var_mtx,varF_mtx] = boolean_RunNorm(x0,x0_F,func,para,Flag,max_T);
    varF(i,:) = varF_mtx(end,sel);
    var(i,1) = bi2de(var_mtx(end,:));
    trans_vec(i,:) = {bi2de(var_mtx),varF_mtx(:,sel)};
end
