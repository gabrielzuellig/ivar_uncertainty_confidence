function X = VARmakelags(DATA,lags)
% =======================================================================
% Builds a matrix with lagged values of DATA, i.e. if DATA = [x y],
% VARmakelags(DATA,1) yields X = [x y x(-1) y(-1)]
% =======================================================================
% X = VARmakelags(DATA,lags)
% -----------------------------------------------------------------------
% INPUT
%   DATA : matrix containing the original data
%   lags : lag order
%
% OUPUT
%   X    : matrix of lagged values
% =======================================================================
% Ambrogio Cesa Bianchi, February 2012
% ambrogio.cesabianchi@gmail.com


[T]= size(DATA,1);

% Create the lagged matrix
X = [];
for jj=0:lags-1
    X = [DATA(jj+1:T-lags+jj,:), X];
end

% Finally, save the non-lagged values...
aux = DATA(lags+1:end,:);

%... and append to the lagged matrix
X = [aux X];
