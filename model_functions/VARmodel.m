function results = VARmodel(DATA,nlag,c_case,extexoshock,DATA_EX,nlag_ex)

% =======================================================================
% Performs vector autogressive estimation
% =======================================================================
% results = VARmodel(DATA,nlag,c_case,DATA_EX)
% -----------------------------------------------------------------------
% INPUTS 
%   DATA = an (nobs x neqs) matrix of y-vectors
%	nlag = the lag length
%
% OPTIONAL INPUTS
%   c_case  = 0 no constant; 1 with constant ; 2 constant + time trend (default = 1)
%	DATA_EX = optional matrix of variables (nobs x nvar_ex)
%
% OUTPUT
%    results.X         = independent variable (ordered as in VARmakexy)
%    results.Y         = dependent variable
%    results.X_EX      = matrix of exogenous variables
%    results.beta      = matrix of estimated coefficients (ordered as in VARmakexy)
%    results.sigma     = VCV matrix of residuals (adjusted for dof)
%    results.nobs      = nobs, # of observations adjusted (i.e., removing the lags)
%    results.neqs      = neqs, # of equations
%    results.nlag      = nlag, # of lags
%    results.nvar      = nlag*neqs, # of endogenous variables per equation
%    results.nvar_ex   = nvar_ex, # of exogenous variables;
%    results.ntotcoeff = nlag*neqs+nvar_ex+c_case , # coefficients to be estimated per per equation
%--------------------------------------------------------------------------
%    results(eq).beta  = bhat for equation eq
%    results(eq).tstat = t-statistics 
%    results(eq).tprob = t-probabilities
%    results(eq).resid = residuals 
%    results(eq).yhat  = predicted values 
%    results(eq).y     = actual values 
%    results(eq).sige  = e'e/(n-nvar)
%    results(eq).rsqr  = r-squared
%    results(eq).rbar  = r-squared adjusted
%    results(eq).boxq  = Box Q-statistics
%    results(eq).ftest = Granger F-tests
%    results(eq).fprob = Granger marginal probabilities
%
% Note: compared to Eviews, there is a difference in the estimation of the 
% constant when lag is > 2. This is because Eviews initialize the trend
% with the number of lags (i.e., when lag=2, the trend is [2 3 ...T]), 
% while VARmakexy.m initialize the trend always with 1.
% =======================================================================
% Ambrogio Cesa Bianchi, May 2012
% ambrogio.cesabianchi@gmail.com



%% Check inputs
%==============
[nobs, neqs] = size(DATA);

% Check if ther are constant, trend, both, or none
if ~exist('extexoshock','var')
    extexoshock = 0;
end

if ~exist('c_case','var')
    c_case = 1;
end

% Check if there are exogenous variables
if exist('DATA_EX','var')
    [nobs2, num_ex] = size(DATA_EX);
    % Check that DATA and DATA_EX are conformable
    if (nobs2 ~= nobs)
        error('var: nobs in DATA_EX-matrix not the same as y-matrix');
    end
    clear nobs2
else
    num_ex = 0;
end

% Check if there is lag orderd of DATA_EX, otherwise set it to 0
if ~exist('nlag_ex','var')
    nlag_ex = 0;
end


%% Save some parameters and create data for VAR estimation
%=========================================================
nobse                 = nobs - max(nlag,nlag_ex);
    results.nobs      = nobse;
    results.neqs      = neqs;
    results.nlag      = nlag;
    results.nlag_ex   = nlag_ex;
nvar                  = neqs*nlag; 
    results.nvar      = nvar;
nvar_ex               = num_ex*(nlag_ex+1);
    results.nvar_ex   = nvar_ex;
ntotcoeff             = nvar + nvar_ex + c_case;
    results.ntotcoeff = ntotcoeff;
    results.c_case    = c_case;


% Create independent vector and lagged dependent matrix
[Y, X] = VARmakexy(DATA,nlag,c_case);

% Create (lagged) exogeanous matrix
if nvar_ex
    X_EX  = VARmakelags(DATA_EX,nlag_ex);
    if nlag == nlag_ex
        disp('hello')
        X = [X X_EX];
    elseif nlag > nlag_ex
        diff = nlag - nlag_ex;
        X_EX = X_EX(diff+1:end,:);
        X = [X X_EX];
    elseif nlag < nlag_ex
        diff = nlag_ex - nlag;
        Y = Y(diff+1:end,:);
        X = [X(diff+1:end,:) X_EX];
    end
end


%% OLS estimation equation by equation
%=====================================

% pull out each y-vector and run regressions
for j=1:neqs;
    aux = ['eq' num2str(j)];
    Yvec = Y(:,j);
    if sum(j==extexoshock)==1
        %{
        % define zero restr.: 0 for all except for own lags
        if c_case == 0
            zerorest(1) = 1;
        elseif c_case == 1
            zerorest(2) = 1;  % unc. series has non-zero level, so need that const. coeff.
        elseif c_case == 2
            zerorest(3) = 1;
        end  
        svtfp is demeaned. If it was not use the following:
        %}
        zerorest = zeros(1, size(X,2));
        if c_case == 1
            zerorest(1) = 1;  % 
        elseif c_case == 2
            zerorest(1:2) = 1;
        end
        for ll = 1:nlag
            zerorest(c_case + extexoshock+((ll-1)*neqs)) = 1;
        end
        
        ols_struct = ols(Yvec,X(:, logical(zerorest)));
        eval( ['results.' aux '.beta = zeros(size(X,2),1) ;'] );
        eval( ['results.' aux '.beta(logical(zerorest), 1)  =  ols_struct.beta ;'] );        % bhats
        eval( ['results.' aux '.tstat = zeros(size(X,2),1) ;'] );
        eval( ['results.' aux '.tstat(logical(zerorest), 1) = ols_struct.tstat ;'] );  % t-stats
        tstat = zeros(sum(zerorest),1);
        tstat = ols_struct.tstat;
        tout = tdis_prb(tstat,nobse-sum(zerorest));
    else
        ols_struct = ols(Yvec,X);
        eval( ['results.' aux '.beta  = ols_struct.beta;'] );        % bhats
        eval( ['results.' aux '.tstat = ols_struct.tstat;'] );       % t-stats
        % compute t-probs
        tstat = zeros(nvar,1);
        tstat = ols_struct.tstat;
        tout = tdis_prb(tstat,nobse-nvar);
    end
    eval( ['results.' aux '.tprob = tout;'] );                   % t-probs
    eval( ['results.' aux '.resid = ols_struct.resid;'] );       % resids 
    eval( ['results.' aux '.yhat = ols_struct.yhat;'] );         % yhats
    eval( ['results.' aux '.y    = Yvec;'] );                    % actual y
    eval( ['results.' aux '.rsqr = ols_struct.rsqr;'] );         % r-squared
    eval( ['results.' aux '.rbar = ols_struct.rbar;'] );         % r-adjusted
    eval( ['results.' aux '.sige = ols_struct.sige;'] );

    % do the Q-statistics
    % use residuals to do Box-Pierce Q-stats
    % use lags = nlag in the VAR
    % NOTE: a rule of thumb is to use (1/6)*nobs but this seems excessive to me
    elag = mlag(ols_struct.resid,nlag);
    % feed the lags
    etrunc = elag(nlag+1:nobse,:);
    rtrunc = ols_struct.resid(nlag+1:nobse,1);
    qres   = ols(rtrunc,etrunc);
    if nlag ~= 1
    	boxq = (qres.rsqr/(nlag-1))/((1-qres.rsqr)/(nobse-nlag));
    else
        boxq = (qres.rsqr/(nlag))/((1-qres.rsqr)/(nobse-nlag));
    end

    eval( ['results.' aux '.boxq = boxq;'] );

    % TO DO: This part, til line 206 should be adjusted for zero-restrictions, but not important!!
    
    % form matrices for joint F-tests (exclude each variable sequentially)
    for r=1:neqs;
        xtmp = [];
        for s=1:neqs
            if s ~= r
                xlag = mlag(DATA(:,s),nlag);
                if nlag == nlag_ex
                    xtmp = [xtmp trimr(xlag,nlag,0)];
                elseif nlag > nlag_ex
                    xtmp = [xtmp trimr(xlag,nlag,0)];
                elseif nlag < nlag_ex
                    xtmp = [xtmp trimr(xlag,nlag+diff,0)];
                end
            end
        end
        % we have an xtmp matrix that excludes 1 variable
        % add deterministic variables (if any) and constant term
        if nvar_ex > 0
            [xtmp X_EX ones(nobse,1)];
        else
            xtmp = [xtmp ones(nobse,1)];
        end
        % get ols residual vector
        b = xtmp\Yvec; % using Cholesky solution
        etmp = Yvec-xtmp*b;
        sigr = etmp'*etmp;
        % joint F-test for variables r
        sigu = ols_struct.resid'*ols_struct.resid;
        ftest(r,1) = ((sigr - sigu)/nlag)/(sigu/(nobse-nvar)); 
    end

    eval( ['results.' aux '.sige = ols_struct.sige;'] );
%     eval( ['results.' aux '.ftest = ftest;' ]);     
%     eval( ['results.' aux '.fprob = fdis_prb(ftest,nlag,nobse-nvar);' ]); 

end % end of loop over equations 

%% Compute the matrix of coefficients & VCV
%==========================================
BETA = (X'*X)\(X'*Y); % ordered by block lags (i.e., all coefficients with first lag, then second,...)
if extexoshock
    for j = extexoshock
        eval( ['BETA(:,' num2str(j) ') = results.eq' num2str(j) '.beta;'] );    % overwrite beta matrix with zero restrictions
    end
end
results.beta = BETA;
SIGMA = (1/(nobse))*(Y-X*BETA)'*(Y-X*BETA); %(nobse-ntotcoeff) adjusted for # of estimated coeff per equation
results.sigma = SIGMA;
results.residuals = Y - X*BETA;

results.X = X;
results.Y = Y;
if nvar_ex > 0
    results.X_EX = X_EX;
end

if exist('DATA_EX','var')
    results.nlagsNL=size(DATA_EX,2);
end

results.extexoshock = extexoshock;
  
