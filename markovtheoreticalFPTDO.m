%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A script to evaluate first passae time distribution (FPTD) based on the
% theoretical equation
%  fptd: FPTD for diffusion with the stochastic diffusion coefficient
% dfptd: FPTD for the corresponding diffusion with ensemble-averaged diffusivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all parameters and close all figures
clear;
close all;

% set parameter values
tp = .5;        % transition probability between the two states
n = 8000;       % the number of sample paths   
l = 200;        % time length for calculation
Du = 0.6;       % diffusion coefficient at + state
Dd = 0.4;       % diffusion coefficient at - state
res = .01;      % time resolution
maxid = 100;    % time length

% simulate stochastic diffusion coefficient
[Dc, leng] = tsMDudO(tp,n,l,Du,Dd,res,maxid);

tt = res:res:maxid;
ltt = length(tt);
Dm = 0.5;

D1 = cell2mat(Dc);
D2 = reshape(D1,maxid/res,[]);
D = D2';
clear D1 D2;

nD = D/Dm;
H = taveragDO(nD, leng);

% estimate first passage time distributions (FPTDs)
ts = 0.1;
Dmd = 1/(16*ts);    

dfptd = ((1./sqrt(4 * pi * Dmd * tt.^3)) .* exp(-1./(4 * Dmd * tt)) * res)';    

fptdw = zeros(ltt,n);
for j = 1:n
    fptdw(:,j) = (nD(j,:)./sqrt(4 * pi * Dmd * H(j,:).^3 .* tt.^3)) .* exp(-1./(4 * Dmd * H(j,:) .* tt)) * res;
end
fptd = mean(fptdw,2);   

























