%% mapping a vector to the reference
function idx = mapping(vec,ref,type_id)
num1 = length(vec);
num2 = length(ref);
idx = zeros(num1,1);
if type_id==1   % string
   for k = 1:num1
       b = find(strcmp(ref,vec{k})==1,1);
       if ~isempty(b)
           idx(k) = b;
       end
   end
else   % number
   for k = 1:num1
       b = find(ref==vec(k),1);
       if ~isempty(b)
           idx(k) = b;
       end
   end    
end
end

