function S = taveragDO(D,leng)
% Function to evaluate the time average of a stochastic diffusion coefficient
%   Input variables
%       D: stochastic diffusion coefficient
%    leng: time length
%   Output variable
%       S: the time average of the diffusion coefficient

n = size(D,1);
SS = cumsum(D,2);
T = 1:leng;
iT = 1./T;
iTm = repmat(iT,n,1);
S = SS.*iTm;

end

