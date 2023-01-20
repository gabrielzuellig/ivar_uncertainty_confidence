function [dIRFinf,dIRFsup,dIRFmed,dIRF]=VAR_intstat(IRFdraws1,IRFdraws2,select,pctg)
% Obtain the difference between state-conditional GIRFs and their confidence band.
% by G. Pellegrino

% check inputs
[a1,a2,a3,a4]=size(IRFdraws1);
[b1,b2,b3,b4]=size(IRFdraws2);
if a1~=b1 
    disp('error: the 2 IRFs must be comparable')
end

if a2~=b2 
    disp('error: the 2 IRFs must be comparable')
end

if a4~=b4 
    disp('error: the 2 IRFs must be comparable')
end

dIRF=zeros(a1,a2,a4); %matrix that will contain the difference among the two IRFs for the shoch selected
%compute the difference between the two IRFs for each draw
for ii=1:a4
    dIRF(:,:,ii)=IRFdraws2(:,:,select,ii)-IRFdraws1(:,:,select,ii);
end
%obtaining the prabability band with percentiles
    pctg_inf = (100-pctg)/2; 
    pctg_sup = 100 - (100-pctg)/2;
    dIRFinf(:,:,:) = prctile(dIRF(:,:,:),pctg_inf,3);
    dIRFsup(:,:,:) = prctile(dIRF(:,:,:),pctg_sup,3);
    dIRFmed(:,:,:) = prctile(dIRF(:,:,:),50,3);

