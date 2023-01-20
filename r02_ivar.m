
%% Initialize
tic
rand('seed',1234)
randn('seed',1234)


%% split data into endogenous and exogenous
data_end=data(:,1:numend); % endogenous data
data_ex=data(:,numend+1:numend+nlags); % exogenous data + interaction terms (remember to put interaction variables are last ordered in data_ex)
data_exFull=data(:,numend+1:numend+12);

[nobs] = size(data_end,1);
labels_ex=labels(1,numend+1:numend+nlags);
labels = labels(1,1:numend);


%% OLS estimation of I-VAR (interactions as exogenous variables)
VAR = VARmodel(data_end,nlags,c_case,extexoshock,data_ex);
VARprint(VAR,labels,labels_ex);


%% Nested linear VAR
if dolinear==1

    % estimation
    VARlin= VARmodel(data_end,nlags,c_case,extexoshock);

    % model fit with IVAR vs. linear VAR
    LMstat=(nobs-nlags)*(log(det(VARlin.sigma))-log(det(VAR.sigma))); % see Canova ch. 4
    dof=nlags*VAR.neqs ; % number of restrictions (additional interacted terms)
    pvalue=1-chi2cdf(LMstat,dof);

    % sample responses:
    history=zeros(1,VARlin.ntotcoeff);
    if proxyvar == 1
        VARlin.z = instr(nlags+1:end);
        VARlin.select = ones(size(VARlin.z));
    end
    [IRFs,IRFs_optavg]=VARgirf(VARlin,nsteps+horplus,pick,history,ordend,ordend2_lin,ndraws,mode,shocksize,1);
    IRFs=IRFs(1:nsteps,:,:);

    % CI with bootstrap
    if dobootstrap == 1

        VARlin.defstate = defstate;
        [INFavg2, SUPavg2, MEDavg2, ~,~,~,INFlindraws,~] = VARbootstrap(VARlin,IRFs_optavg,pick,horplus,0,ndraws,confLev2,shocksize,1);

        if irfmedian==1
            IRFs = MEDavg2;
        end
        if irfmeanfromdraws==1
            IRFs = mean(INFlindraws, 4);
        end

        % Figure
        fig=figure;
        VARirplot(IRFs,pick,labels_print,INFavg2,SUPavg2)
        tightfig;
        set(gcf,'paperpositionmode','auto')
        print(gcf,'-depsc2','-loose',strcat('./',spec,'/irf_lin.eps'));

    end

end


%% Split histories into two different states

state = data(nlags:end-1, ordend2);
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
stateH = state > prctile(state, perc);
stateL = state <= prctile(state, perc);

historiesH = VAR.X(stateH, :);
nhistH=size(historiesH,1);

historiesL = VAR.X(stateL, :);
nhistL=size(historiesL,1);

val = prctile(state, perc); % Save cutoff percentile and absolute value


%% State-dependent GIRFs, point estimates

% for sample state-conditional GIRFs
if proxyvar == 1
    VAR.z = instr(nlags+1:end);
    VAR.select = ones(size(VAR.z));
end
[OIRFavg,OIRF_optavg]=VARgirf(VAR,nsteps+horplus,pick,historiesL,ordend,ordend2,ndraws,mode,shocksize,0);
[OIRF2avg,OIRF_opt2avg]=VARgirf(VAR,nsteps+horplus,pick,historiesH,ordend,ordend2,ndraws,mode,shocksize,0);
OIRFavg=OIRFavg(1:nsteps,:,:);
OIRF2avg=OIRF2avg(1:nsteps,:,:);

% Diagnostic explosiveness/instability:
disp('Diagnostics for GIRFS point estimates:')
disp(['- Number of explosive initial histories discarded: ' num2str(OIRF_optavg.nexplos) ' for pessimistic state, ', num2str(OIRF_opt2avg.nexplos) ' for normal state.'])
if max(OIRF_optavg.countRep)>0
    disp('- It was needed to repeat some particularly extreme residual extractions for at least one initial history. For details see OIRF_optavg.countRep and OIRF_opt2avg.countRep')
end


%% Generalized Forecast error variance decomposition
if dofevd == 1
    clear OIRFhist
    for j = 1:numend
        [~,~,OIRFhist_temp]=VARgirf(VAR,nsteps+horplus,j,historiesL,ordend,ordend2,ndraws,mode,shocksize,0);
        OIRFhist(:,:,j,:)= squeeze(OIRFhist_temp(:,:,j,:));
        [~,~,OIRFhist_temp]=VARgirf(VAR,nsteps+horplus,j,historiesH,ordend,ordend2,ndraws,mode,shocksize,0);
        OIRFhist2(:,:,j,:)= squeeze(OIRFhist_temp(:,:,j,:));
    end
    clear OIRFhist_temp

    allfevd = NaN*ones(numend,numend,2);  % rows: shock variable, column: response variable, 3rd: state
    for j = 1:numend
        % low state:
        CSshock=cumsum(OIRFhist(:,:,j,:).^2,1) ; % cumsum of girfs for h=1:nstep (numerator formula 9 in Lanne et al, 2016)
        CSall= sum(cumsum(OIRFhist(:,:,:,:).^2,1),3) ;
        fevd_hist_pick=squeeze(CSshock)./squeeze(CSall); % contains the gfevd for each step ahead (dim1), each variable (dim2) and for each history (dim3), for the picked shock.
        fevd=nanmean(fevd_hist_pick,3);
        fevd = fevd(1:nsteps, :);  % nsteps x numend
        allfevd(j,:,1) = fevd(25,:);
        % high state
        CSshock=cumsum(OIRFhist2(:,:,j,:).^2,1) ; % cumsum of girfs for h=1:nstep (numerator formula 9 in Lanne et al, 2016)
        CSall= sum(cumsum(OIRFhist2(:,:,:,:).^2,1),3) ;
        fevd_hist_pick=squeeze(CSshock)./squeeze(CSall); % contains the gfevd for each step ahead (dim1), each variable (dim2) and for each history (dim3), for the picked shock.
        fevd2=nanmean(fevd_hist_pick,3);
        fevd2 = fevd2(1:nsteps, :);  % nsteps x numend
        allfevd(j,:,2) = fevd2(25,:);
        if j == pick
            close all
            VARplotSDirf(IRFs,fevd2,fevd,reverseHL,1, labels_print)
            tightfig;
            set(gcf,'paperpositionmode','auto')
            print(gcf,'-depsc2','-loose',strcat('./',spec,'/fevd.eps'));
        end
        clear fevd2 fevd CSshock CSall fevd_hist_pick
    end
    save(strcat(spec,'/fevd.mat'), 'allfevd')
    clear OIRFhist OIRFhist2
end


%% Bootstrapped confidence intervals
if dobootstrap == 1

    OIRF_optavg.nsteps=nsteps+horplus;
    VAR.defstate = defstate;
    [INFavgL2, SUPavgL2, MEDavgL2, INFavgH2, SUPavgH2, MEDavgH2, IRFdrawsL, IRFdrawsH] = VARbootstrap(VAR,OIRF_optavg,pick,horplus,val,ndraws,confLev2,shocksize,fast);

    if irfmedian==1
        OIRF2avg = MEDavgH2;
        OIRFavg = MEDavgL2;
    end
    if irfmeanfromdraws==1
        OIRF2avg = nanmean(IRFdrawsH, 4);
        OIRFavg = nanmean(IRFdrawsL, 4);
    end

    % Figure
    fig=figure;
    VARplotSDirf(IRFs,OIRF2avg,OIRFavg,reverseHL,pick, labels_print,INFavgH2,SUPavgH2,INFavgL2,SUPavgL2)
    tightfig;
    set(gcf,'paperpositionmode','auto')
    print(gcf,'-depsc2','-loose',strcat('./',spec,'/irf_nonlin.eps'));

    % Difference between bootstrapped confidence intervals
    dOIRFavg=OIRFavg(:,:,pick)-OIRF2avg(:,:,pick);
    [dIRFinf,dIRFsup,dIRFmed,dIRFdraws]=VAR_intstat(IRFdrawsH,IRFdrawsL,pick,confLev);
    [dIRFinf2,dIRFsup2,dIRFmed2,dIRFdraws2]=VAR_intstat(IRFdrawsH,IRFdrawsL,pick,confLev2);

    if irfmedian==1
       dOIRFavg = dIRFmed;
    end

    % Figure
    fig=figure;   %
    VARplotdiff(dOIRFavg,dOIRFavg,pick,labels_print,dIRFinf,dIRFsup,dIRFinf2,dIRFsup2)
    tightfig;
    set(gcf,'paperpositionmode','auto')
    print(gcf,'-depsc2','-loose',strcat('./',spec,'/irf_diff.eps'));

end


%% Save
if dobootstrap == 1
    save(strcat(spec,'/out.mat'), 'VAR', 'OIRF2avg', 'OIRFavg', 'INFavgL2', ...
        'SUPavgL2', 'MEDavgL2', 'INFavgH2', 'SUPavgH2', 'MEDavgH2', 'IRFdrawsL',...
        'IRFdrawsH', 'IRFs', 'INFavg2', 'SUPavg2', 'MEDavg2', 'INFlindraws',...
        'dOIRFavg','dIRFinf','dIRFinf2','dIRFsup','dIRFsup2')
else
    save(strcat(spec,'/out_nobtstrp.mat'), 'VAR', 'OIRF2avg', 'OIRFavg', 'IRFs')
end

close all;

toc

