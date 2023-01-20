
nhist_L=size(historiesL,1);
nhist_H=size(historiesH,1);


% find the G-FEVD for each single history in the ZLB regime
CSshock=cumsum(OIRFhist(:,:,pick,:).^2,1) ; % cumsum of girfs for h=1:nstep (numerator formula 9 in Lanne et al, 2016)
CSall= sum(cumsum(OIRFhist(:,:,:,:).^2,1),3) ; %  h-step forecast error variance (for h=1:nstep), for each variable (denomimator formula 9 in Lanne et al, 2016)

    %----
    % h-step ahead forecast error standard error
    fese=sqrt(CSall);% this, for the linear model, coincides with the output of GRETL when typeImpulse=0 (i.e., when 1 st.dev. shock)
    %----

fevd_hist_pick=squeeze(CSshock)./squeeze(CSall); % contains the gfevd for each step ahead (dim1), each variable (dim2) and for each history (dim3) , for the picked shock. 
 % Notice: It is a number between 0 and 1 giving the contributions to each variable of uncertainty shocks 

% ---- proof that gfevd for each shock sum to 1
% it can just be seen from computation above that holds by construction. 
% CSshock_all=cumsum(OIRFhist(:,:,:,:).^2,1) ;
% sum_CSshock_all=sum(CSshock_all,3); % sum the contribution of each shock (think to the sum of the denominator of formula 9 in Lanne et al. for each shock)
% sum_fevd_hist=sum_CSshock_all./CSall; % --> this is a matrix full of ones!
% ---- 

% find the average G-FEVD for the ZLB regime, and plot it
fevd_ZLB_pick=mean(fevd_hist_pick,3);


%%%%%%%
% find the G-FEVD for each single history in the Normal times regime
CSshock_2=cumsum(OIRFhist2(:,:,pick,:).^2,1) ; 
CSall_2= sum(cumsum(OIRFhist2(:,:,:,:).^2,1),3) ; 

fevd_hist_2_pick=squeeze(CSshock_2)./squeeze(CSall_2);  

% find the average G-FEVD for the Normal times regime, and plot it
fevd_Normal_pick=mean(fevd_hist_2_pick,3);



disp('Contribution uncertainty shocks:')
disp('---- 3 years ahead (normal times, zlb) ----')
labels'
[fevd_Normal_pick(12,:)',fevd_ZLB_pick(12,:)']
