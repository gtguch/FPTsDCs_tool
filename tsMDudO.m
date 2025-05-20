function [D, leng] = tsMDudO(tp,n,l,Du,Dd,res,maxid)
% Function to generate sample paths of the two-state Markov process
%   Input variables
%       tp: transition probability between the two states
%        n: the number of sample paths
%        l: time length for calculation
%       Du: diffusion coefficient at + state
%       Dd: diffusion coefficient at - state
%      res: time resolution
%    maxid: time length for output
%   Output variables
%        D: diffusion coefficient (two-state Markov process)
%     leng: the number of elements in the time direction for each sample
%           path

tau = exprnd(1/tp,n,2*l)/res;

D = cell(1,n);
tauid = round(cumsum(tau,2));
a = zeros(1,n);
for i = 1:n
    b = find(tauid(i,:)>maxid/res,1);
    if isempty(b)
        a(i) = 2*l;
    else
        a(i) = b;
        tauid(i,b) = maxid/res;
    end
end    

leng = zeros(1,n);
for i = 1:n
    disp(i)
    Dm = Du*ones(1,tauid(i,a(i)));
    if i <= n/2
        for j = 2:2:a(i)
            Dm(tauid(i,j-1)+1:tauid(i,j)) = Dd;
        end
    else
        Dm(1:tauid(i,1)) = Dd;
        for j = 3:2:a(i)
            Dm(tauid(i,j-1)+1:tauid(i,j)) = Dd;
        end
    end
    D{1,i} = single(Dm);
    leng(i) = length(Dm);
    clear Dm
end 
end