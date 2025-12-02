%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A script to evaluate the cumulative probability of early arriving particles
% ppp: the cumulative probability of early arriving particles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all parameters and close all figures
clear;
close all;

% set parameter values
tp = .5;        % transition probability between the two states
n = 8000;       % the number of sample paths
l = 2000;       % time length for calculation
Du = 1;         % diffusion coefficient at + state
Dd = 0.2;       % diffusion coefficient at - state
res = 0.01;     % time resolution
maxid = 1000;   % time length

% calculate the cumulative probability of early arriving particles
[Dc, leng] = tsMDudO(tp,n,l,Du,Dd,res,maxid);

tt = res:res:maxid;
Dm = (Du + Dd)/2;
hp = Du/Dm;
hm = Dd/Dm;

D1 = cell2mat(Dc);
D2 = reshape(D1,maxid/res,[]);
D = D2';
clear D1 D2;

nD = D/Dm; 

H = taveragDO(nD, leng);
erfcH = erfc((2)./sqrt(H));
ppp = mean(erfcH);





