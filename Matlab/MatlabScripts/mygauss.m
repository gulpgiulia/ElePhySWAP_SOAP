function Y = mygauss(Params, X)
%
%  Y = mygauss(Params, X)
%
%
%   Copyright 2008 Maurizio Mattia @ Ist. Super. Sanità, Rome - Italy
%   Version: 1.0 - Jun. 24, 2008
%

Mu = Params(1);
Sigma = Params(2);
Ampl = Params(3);

Y = Ampl * exp(-(X-Mu).^2/(2*Sigma^2));