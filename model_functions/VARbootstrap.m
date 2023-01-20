function [INFavgL,SUPavgL,MEDavgL,INFavgH,SUPavgH,MEDavgH,IRFdraws_histL,IRFdraws_histH] = VARirbandtrue_endint1shock_historiesEnd2_futureN(VAR,OIRF_optavg,pick,horplus,val,ndraws,pctg,typeImpulse,fast)
% =======================================================================
%% Retrieve some parameters from the structure VAR and IRF_opt
% =============================================================

nsteps = OIRF_optavg.nsteps;
ordend= OIRF_optavg.ordend;
ordend2= OIRF_optavg.ordend2;

beta     = VAR.beta;  % rows are coefficients, columns are equations
nvars    = VAR.neqs;
nvar_ex  = VAR.nvar_ex;
nlags    = VAR.nlag;
c_case   = VAR.c_case;
nobs     = VAR.nobs;
Y        = VAR.Y;
resid    = VAR.residuals;
extexoshock = VAR.extexoshock;
defstate = VAR.defstate;

if pick~=ordend
    error('Error: pick must be equal to the variable used as ordend in the VAR.')
end

numExpl=0;

mode=OIRF_optavg.mode;
draws=OIRF_optavg.draws;

IRFdraws_histL=NaN(nsteps,nvars,nvars,ndraws); 
IRFdraws_histH=NaN(nsteps,nvars,nvars,ndraws); 


%% Create the matrices for the loop
%==================================

%% Loop over the number of draws, generate data, estimate the var and then
%% calculate impulse response functions
%==========================================================================
tt = 1; %cycle over draws
while tt<=ndraws
    
    %% STEP 1: use the Bootstrap method to generate the residuals
    
    y_artificial = zeros(nobs,nvars);
    
          % Use the residuals to bootstrap: generate a random number bounded
        % between 0 and # of residuals, then use the ceil function to select
        % that row of the residuals (this is equivalent to sampling with replacement)
        u = resid(ceil(size(resid,1)*rand(nobs,1)),:);
        

    
    %%
    
    %%IMPORTANT:
    % to have the right results when there are interactions between
    % endogenous variables one has to construct properly also the
    % series for the interaction variables from each draw. It is important because the
    % interaction variables have to respect the link that they have
    % with our chosen endogenous variables.
    exog_draw=zeros(nobs,nvar_ex);
    
    %% STEP 2: generate the artifcial data
    
    %% STEP 2.1: initial values for the artificial data
    % Intialize the first nlags observations with real data + plus artificial
    % residuals. Nonetheless, in the estimation of the VAR on the simulated data,
    % I through away the first nobs observations so it should not matter.
    
    LAG=[];
    
        ii=1;
    
    for jj = 1:nlags
        y_artificial(jj,:) = Y(ii+jj-1,:) + u(jj,:);
        LAG = [y_artificial(jj,:) LAG];
        % Initialize the artificial series and the LAGplus vector
        if c_case==0
            LAGplus = LAG;
        elseif c_case==1
            LAGplus = [1 LAG];
        elseif c_case==2
            T = (1:nobs)';
            LAGplus = [1 T(jj) LAG] ;
        end
    end
    
    
    if nvar_ex~=0 % to contruct the interaction variables
        if nlags==1
            LAGplus = [LAGplus LAG(1,ordend2)*LAG(1,ordend)];
            
        elseif nlags==2
            LAGplus = [LAGplus LAG(1,ordend2)*LAG(1,ordend) LAG(1,ordend2+nvars)*LAG(1,ordend+nvars)];
        elseif nlags==3
            LAGplus = [LAGplus LAG(1,ordend2)*LAG(1,ordend) LAG(1,ordend2+nvars)*LAG(1,ordend+nvars) LAG(1,ordend2+2*nvars)*LAG(1,ordend+2*nvars)];
        elseif nlags==4
            LAGplus = [LAGplus LAG(1,ordend2)*LAG(1,ordend) LAG(1,ordend2+nvars)*LAG(1,ordend+nvars) LAG(1,ordend2+2*nvars)*LAG(1,ordend+2*nvars) LAG(1,ordend2+3*nvars)*LAG(1,ordend+3*nvars)];
        elseif nlags==5
            LAGplus = [LAGplus LAG(1,ordend2)*LAG(1,ordend) LAG(1,ordend2+nvars)*LAG(1,ordend+nvars) LAG(1,ordend2+2*nvars)*LAG(1,ordend+2*nvars) LAG(1,ordend2+3*nvars)*LAG(1,ordend+3*nvars) LAG(1,ordend2+4*nvars)*LAG(1,ordend+4*nvars)];
        elseif nlags==6
            LAGplus = [LAGplus LAG(1,ordend2)*LAG(1,ordend) LAG(1,ordend2+nvars)*LAG(1,ordend+nvars) LAG(1,ordend2+2*nvars)*LAG(1,ordend+2*nvars) LAG(1,ordend2+3*nvars)*LAG(1,ordend+3*nvars) LAG(1,ordend2+4*nvars)*LAG(1,ordend+4*nvars) LAG(1,ordend2+5*nvars)*LAG(1,ordend+5*nvars)];
        else disp('Error: for now it is possible only to have nlags up to 4.')
        end
    end
    
    
    
    
    %% STEP 2.2: generate artificial series
    % From observation nlags+1 to nobs, compute the artificial data
    for jj = nlags+1:nobs
        
        for mm = 1:nvars
            % Compute the value for time=jj
            y_artificial(jj,mm) = LAGplus * beta(1:end,mm) + u(jj,mm);
        end
        % now update the LAG matrix
        LAG = [y_artificial(jj,:) LAG(1,1:(nlags-1)*nvars)];
        if c_case==1
            LAGplus = [1 LAG];
        elseif c_case==2
            LAGplus = [1 T(jj) LAG];
        end
        
        if nvar_ex~=0 % to construct properly also the interaction variables
            if nlags==1
                LAGplus = [LAGplus LAG(1,ordend2)*LAG(1,ordend)];
            elseif nlags==2
                LAGplus = [LAGplus LAG(1,ordend2)*LAG(1,ordend) LAG(1,ordend2+nvars)*LAG(1,ordend+nvars)];
            elseif nlags==3
                LAGplus = [LAGplus LAG(1,ordend2)*LAG(1,ordend) LAG(1,ordend2+nvars)*LAG(1,ordend+nvars) LAG(1,ordend2+2*nvars)*LAG(1,ordend+2*nvars)];
            elseif nlags==4
                LAGplus = [LAGplus LAG(1,ordend2)*LAG(1,ordend) LAG(1,ordend2+nvars)*LAG(1,ordend+nvars) LAG(1,ordend2+2*nvars)*LAG(1,ordend+2*nvars) LAG(1,ordend2+3*nvars)*LAG(1,ordend+3*nvars)];
            elseif nlags==5
                LAGplus = [LAGplus LAG(1,ordend2)*LAG(1,ordend) LAG(1,ordend2+nvars)*LAG(1,ordend+nvars) LAG(1,ordend2+2*nvars)*LAG(1,ordend+2*nvars) LAG(1,ordend2+3*nvars)*LAG(1,ordend+3*nvars) LAG(1,ordend2+4*nvars)*LAG(1,ordend+4*nvars)];
            elseif nlags==6
                LAGplus = [LAGplus LAG(1,ordend2)*LAG(1,ordend) LAG(1,ordend2+nvars)*LAG(1,ordend+nvars) LAG(1,ordend2+2*nvars)*LAG(1,ordend+2*nvars) LAG(1,ordend2+3*nvars)*LAG(1,ordend+3*nvars) LAG(1,ordend2+4*nvars)*LAG(1,ordend+4*nvars) LAG(1,ordend2+5*nvars)*LAG(1,ordend+5*nvars)];
            else
                disp('Error: for now it is possible only to set nlags up to 2.')
            end
            
        end
    end
    

    numb=isnan(y_artificial(nobs,ordend));
    if numb~=1    % if one bootstrap draw delivers an explosive series for the endogenous
    % variable the draw is repeated with a new sequence of drawn residuals 
        % sometimes, however, the artificial sample although divergirg,
        % still doesn't take NaN values for the number of iterations considered. In
        % these cases the final simulated observations take extreme values with
        % respect to initial values, so that to give problems in the
        % VAR estimation step. Hence the bootstrap draw necessitates to be repeated. 
        % Instead than using an arbitrary criterion, we discard just the
        % draws giving problems in the VARmodel estimation function. In
        % particular: i) we repeat the draw if errors or specific warning
        % and ii) repeat the draw if all histories are discarded.
        
        warning('error', 'MATLAB:nearlySingularMatrix'); % turn this specific warning into an error. 
        try 
            %% STEP 3: estimate VAR on artificial data.
            
            % here, after having simulated a series for the endogneous
            % variables, we also obtain the corresponding simulated series
            % for the interaction variables
            if nvar_ex~=0
                for ii=nlags+1:nobs
                    if nlags==1
                        exog_draw(ii,:)=y_artificial(ii-1,ordend2)*y_artificial(ii-1,ordend);
                    elseif nlags==2
                        exog_draw(ii,:)=[y_artificial(ii-1,ordend2)*y_artificial(ii-1,ordend) y_artificial(ii-2,ordend2)*y_artificial(ii-2,ordend)];
                    elseif nlags==3
                        exog_draw(ii,:)=[y_artificial(ii-1,ordend2)*y_artificial(ii-1,ordend) y_artificial(ii-2,ordend2)*y_artificial(ii-2,ordend) y_artificial(ii-3,ordend2)*y_artificial(ii-3,ordend)];
                    elseif nlags==4
                        exog_draw(ii,:)=[y_artificial(ii-1,ordend2)*y_artificial(ii-1,ordend) y_artificial(ii-2,ordend2)*y_artificial(ii-2,ordend) y_artificial(ii-3,ordend2)*y_artificial(ii-3,ordend) y_artificial(ii-4,ordend2)*y_artificial(ii-4,ordend)];
                    elseif nlags==5
                        exog_draw(ii,:)=[y_artificial(ii-1,ordend2)*y_artificial(ii-1,ordend) y_artificial(ii-2,ordend2)*y_artificial(ii-2,ordend) y_artificial(ii-3,ordend2)*y_artificial(ii-3,ordend) y_artificial(ii-4,ordend2)*y_artificial(ii-4,ordend) y_artificial(ii-5,ordend2)*y_artificial(ii-5,ordend)];
                    elseif nlags==6
                        exog_draw(ii,:)=[y_artificial(ii-1,ordend2)*y_artificial(ii-1,ordend) y_artificial(ii-2,ordend2)*y_artificial(ii-2,ordend) y_artificial(ii-3,ordend2)*y_artificial(ii-3,ordend) y_artificial(ii-4,ordend2)*y_artificial(ii-4,ordend) y_artificial(ii-5,ordend2)*y_artificial(ii-5,ordend) y_artificial(ii-6,ordend2)*y_artificial(ii-6,ordend)];
                    end
                end
                
                VAR_draw = VARmodel(y_artificial(nlags+1:end,:),nlags,c_case,extexoshock,exog_draw(nlags+1:end,:));% the first observaation are not considered
                 
            else
                VAR_draw = VARmodel(y_artificial(1:end,:),nlags,c_case,extexoshock);
            end
            
            
            %% STEP 4: calculate "ndraws" impulse responses and store them
            
            if nvar_ex~=0
                %algoritm to select all histories for the Low and High state:
                
                clear historiesL historiesH
                
                state = VAR_draw.X(1:(end-1),c_case+ordend2);
                
                if defstate.growth == 1
                    state = [NaN; state(2:end) - state(1:end-1)];
                end
                if defstate.cumul > 0
                    helpstate = NaN*ones(size(state,1),defstate.cumul);
                    helpstate(:,1) = state;
                    for i = 2:defstate.cumul
                        helpstate(:,i) = [repmat(NaN,i-1,1); state(1:(end-i+1))];
                    end
                    state = sum(helpstate,2);
                end
                
                stateH = logical(state > val);
                mean(stateH)
                stateL = logical(state <= val);
                
                historiesH = VAR_draw.X(stateH,:);                
                historiesL = VAR_draw.X(stateL,:);
                
                [irf_drawL, irf_draw_optL] = VARgirf(VAR_draw,nsteps,pick,historiesL,ordend,ordend2,draws,mode,typeImpulse,fast);
                [irf_drawH, irf_draw_optH] = VARgirf(VAR_draw,nsteps,pick,historiesH,ordend,ordend2,draws,mode,typeImpulse,fast);
            

            else
                historiesL=zeros(1,VAR_draw.ntotcoeff);
                historiesH=zeros(1,VAR_draw.ntotcoeff);
                history = VAR_draw.X;
                [irf_drawL, irf_draw_optL] = VARgirf(VAR_draw,nsteps,pick,history,ordend,0,draws,mode,typeImpulse,fast);
                irf_drawH = irf_drawL;
                irf_draw_optH = irf_draw_optL;
                
            end
            
            disp(['Loop ' num2str(tt) ' / ' num2str(ndraws) ' draws'])
            
            if irf_draw_optL.nexplos~=length(historiesL) && irf_draw_optH.nexplos~=length(historiesH)
                IRFdraws_histL(:,:,:,tt)=irf_drawL;
                IRFdraws_histH(:,:,:,tt)=irf_drawH;
                tt = tt+1;
            else
                disp(['because of too large values of y_artificial repeat draw #: ' num2str(tt) ' (since all histories were discarded)' ]);
                numExpl=numExpl+1;
            end
            
        catch 
            numExpl=numExpl+1;
        end
        
        warning ('on','all'); % restore warning state
        
    else 
        numExpl=numExpl+1;
    end
    
end


%% Compute the error bands
%=========================
%  using boostrap, use percentile (upper and lower bounds) bands type
% 1) CASE LOW PERCENTILE

INFavgL = zeros(nsteps,nvars,nvars);
SUPavgL = zeros(nsteps,nvars,nvars);
MEDavgL = zeros(nsteps,nvars,nvars);

pctg_inf = (100-pctg)/2;
pctg_sup = 100 - (100-pctg)/2;

IRFdraws_histL=IRFdraws_histL(1:nsteps-horplus,:,:,:);

INFavgL = prctile(IRFdraws_histL(:,:,:,:),pctg_inf,4);
SUPavgL = prctile(IRFdraws_histL(:,:,:,:),pctg_sup,4);
MEDavgL = prctile(IRFdraws_histL(:,:,:,:),50,4);

%%% 2) CASE HIGH PERCENTILE

INFavgH = zeros(nsteps,nvars,nvars);
SUPavgH = zeros(nsteps,nvars,nvars);
MEDavgH = zeros(nsteps,nvars,nvars);

%do a cycle here over histories TO GET BOOSTRAPPED GIRFS FOR EACH HISTORY

pctg_inf = (100-pctg)/2;
pctg_sup = 100 - (100-pctg)/2;

IRFdraws_histH=IRFdraws_histH(1:nsteps-horplus,:,:,:);

INFavgH = prctile(IRFdraws_histH(:,:,:,:),pctg_inf,4);
SUPavgH = prctile(IRFdraws_histH(:,:,:,:),pctg_sup,4);
MEDavgH = prctile(IRFdraws_histH(:,:,:,:),50,4);


disp('the percentage of draws repeated is of:')
percExpl=numExpl/ndraws
if percExpl>=0.25
    disp('Many discarded draws in the bootstrap. SUGGESTION: Start the sample in another date might help')
end
