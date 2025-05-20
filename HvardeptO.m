%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A script to evaluate the time dependence of the variance of H
% Diffusion coefficient at + state is fixed at 1. 
% Hvar: the variance of H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all parameters and close all figures
clear;
close all;

% set parameter values
ap = .5;        % transition probability to + state
am = .5;        % transition probability to - state
Dr = 0.2;       % diffusion coefficient at - state
res = 0.01;     % time resolution
maxid = 1000;   % time length

% calculate the time dependence of the variance of H
tt = res:res:maxid;
Dm = (1 + Dr)/2;
hp = 1/Dm;
hm = Dr/Dm;

Hvar = (hp - hm)^2 * ((ap + am) * tt + exp(-(ap + am) * tt) - 1)./(2 * (ap + am)^2 * (tt.*tt));
