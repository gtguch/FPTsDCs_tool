%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A script to perform simulation for diffusion with a stochastic diffusion
% coefficient and to estimate first passage time
% Solve the Langevin equation numerically using the Euler method
%   D: diffusion coefficient
%   x: position of a diffusing particle
% fpt: first passage time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all parameters and close all figures
clear;
close all;

% set parameter values
tp = .5;        % transition probability between the two states
n = 10000;      % the number of sample paths 
l = 5;          % time length for calculation
Du = 1;         % diffusion coefficient at + state
Dd = 0.2;       % diffusion coefficient at - state
res = 0.001;    % time resolution
maxid = 1;      % time length

% simulate stochastic diffusion coefficient
[Dc, leng] = tsMDudO(tp,n,l,Du,Dd,res,maxid);

D1 = cell2mat(Dc);
D2 = reshape(D1,maxid/res,[]);
D = D2';
clear D1 D2;

% simulate diffusion
dt = .001;            
D = double(D);
[n1,n2] = size(D);
t = dt:dt:size(D,2)*dt;

rv = randn(n1,n2);

xx = sqrt(2*D*dt).*rv; 
x = cumsum(xx,2);

% estimate first passage time
simn = size(x,1);

Dm = (Du + Dd)/2;
ts = [0.1 1];
Dmd = 1./(16*ts);
lDmd = length(Dmd);

nx = x/sqrt(Dm);
clear x;
fpt = zeros(simn,lDmd);
for i = 1:lDmd
    disp(i)
    x = nx * sqrt(Dmd(i));
    bx = imbinarize(x,1);
    dbx = diff(bx,1,2);
    for j = 1:simn
        a = find(dbx(j,:));
        if isempty(a)
            a = NaN;
            fpt(j,i) = a;
        else    
            fpt(j,i) = t(a(1));
        end    
    end
    clear x bx dbx;
end


