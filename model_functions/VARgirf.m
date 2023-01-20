function [OIRFavg,OIRF_optavg,OIRFhist]=VARirtrue_sim1shock_histories_futureN(VAR,nsteps,pick,histories,ordend,ordend2,draws,mode,typeImpulse,fast)
% =======================================================================
%

%% Check inputs

if ~exist('fast','var')
    fast = 0; 
end 
%% Retrieve parameters and preallocate variables
c_case   = VAR.c_case;
nvars    = VAR.neqs;
nlags    = VAR.nlag;
beta     = VAR.beta;
sigma    = VAR.sigma;
if isfield(VAR, 'z')
    z = VAR.z(logical(VAR.select));
    proxyvar = 1;
    disp('Instrumental variable identification')
else
    proxyvar = 0;
    s=NaN;
end

ntotcoeff=VAR.ntotcoeff ; % # coefficients to be estimated per equation

 [Cho, chol_flag] = chol(sigma);
    S = Cho';
    if chol_flag~=0
        disp('-------------------------------------');
        disp('Error: VCV is not positive definite');
    end

    if fast==1
      draws=2;% in this case no extraction future shocks, and only 2 draws are sufficient (not one because it should be an even n°)
    end
    
resid=VAR.residuals;
if proxyvar == 1  %instrumental variable identification

    s = NaN*ones(1,nvars);
    b = logical(1:nvars == pick);
    
    eps_p=resid(logical(VAR.select),b);  % residuals of uncertainty variable
    eps_q=resid(logical(VAR.select),~b);  % all other residuals
    
    % 1st stage
    % beta1 = X'*X \ X'*eps_p;  
    stage1 = fitlm(z, eps_p);
    beta1 = table2array(stage1.Coefficients(1:2,1));  % elast.: +165
    X = [ones(size(z)), z];
    eps_p_hat = X * beta1;  % fitted values of residuals
    regsummary = anova(stage1,'summary');
    Fstat = table2array(regsummary(2,4));
    
    % 2nd stage
    s(~b) = (eps_p_hat'*eps_p_hat)\eps_p_hat'*eps_q;
    s(b) = 1;
    
    % so far, s gives us the vector of relative responses. We want absolute values for a 1sd
    % shock. To get that, follow Piffer (2020):
    sigma11 = sigma(b,b);
    sigma12 = sigma(b,~b);
    sigma21 = sigma(~b,b);
    sigma22 = sigma(~b,~b);
    mu = s(~b)';
    
    Gamma = sigma22 + mu*sigma11*mu' - sigma21*mu' - mu*sigma21';
    
    b11b11prime = sigma11 - (sigma21-mu*sigma11)'*Gamma^(-1)*(sigma21-mu*sigma11);
    b11 = chol(b11b11prime)';
    
    s = s*b11;
  
end


nexplos=0; % to count for the number of explosive histories
timeExp=[];
%% Compute the impulse responses

jj=pick;   
nhist=size(histories,1);
countRep=zeros(1,nhist);
OIRFhist=zeros(nsteps,nvars,nvars,nhist);



parfor ii=1:nhist

    %1) construct (non-orthogonalized) responses
   
    IRF      = zeros(nsteps,nvars,nvars,draws);
    init_shock=[];
    
    % recovering the time which history influences
    time=histories(ii,2);
    
    %give some initial values to the residuals       
        epsilon=zeros(1,nvars);
        Imp_eps=zeros(1,nvars);       
        u= zeros(nsteps,nvars);
        
        ic=1;
        while ic<=draws
                
            if fast==0
                if mode==0
                    tempp=rand(nsteps,1);
                    u = resid(ceil(size(resid,1)*tempp),:);% in this case future and current residuals/shocks are extracted from the empirical distribution of residuals, with replacement. Notice that the distribution of residuals is assumed joint (so that to respect the contemporaneous  relationship among them)
                else
                    u=randn(nsteps,nvars)*Cho  ; % to extract from a multivariate normal distribution respecting the VCV matrix (seen from the help for the function "randn")
                end
            end
        
            X=histories(ii,:);     %predefinite matrix for the starting values for the RHS variables
            W=histories(ii,:);
            
            epsilon(1,:)=(u(1,:)');% to consider reduced form residuals

            if proxyvar == 0

                Imp_eps(1,:)= (S\u(1,:)')' ; % to consider the orthogonalized shock
                
                if typeImpulse==0
                    init_shock=1 ; %initial shock (in terms of st.dev.)
                    Imp_eps(1,pick)=Imp_eps(1,pick)+init_shock;% to ADD a shock of 1 st.dev.
                elseif typeImpulse==1
                    init_shock=(1/S(pick,pick)) ; %initial shock (in terms of st.dev.)
                    Imp_eps(1,pick)=Imp_eps(1,pick)+init_shock;
                else
                    init_shock=typeImpulse;
                    Imp_eps(1,pick)=Imp_eps(1,pick)+init_shock;
                end
                
                Imp_eps(1,:)=S*Imp_eps(1,:)' ; % to re-get again the reduced form residuals to agevolate computations
            
            else
                
                Imp_eps(1,:) = epsilon(1,:);
                
                if typeImpulse==0   % 1 standard deviation
                    init_shock = 1;
                    Imp_eps(1,:) = Imp_eps(1,:) + init_shock*s;
                elseif typeImpulse == 1  % unit increase
                    init_shock = 1/s(pick);
                    Imp_eps(1,:) = Imp_eps(1,:) + init_shock*s;
                else % x standard deviations
                    init_shock = typeImpulse;
                    Imp_eps(1,:) = Imp_eps(1,:) + init_shock*s;
                end
                
            end
            
            Y=zeros(nsteps,nvars);%to store the  responses to the only shock considered by the cycle
            Z=zeros(nsteps,nvars); %to store the  responses when no shock is considered
            
            Y(1,:)=X*beta+Imp_eps; %to get the response of Y when the shock is imposed to an underlying sequence of shocks
            Z(1,:)=W*beta+epsilon; % to get the response of Z for a particular sequence of shocks estracted
            IRF(1,:,jj,ic)=Y(1,:)-Z(1,:) ;%response in the period 0
            
          
            
            %NOW we need to compute the evolution of the system from here on
            %for the case of the shock and for the case of no shock. Then we will get the empirical response by taking the
            % difference between these two.
            
            
            for tt=1:nsteps-1
                
                NewX=zeros(1,ntotcoeff);
                NewW=zeros(1,ntotcoeff);
                
                %New lines to take into account also of the constant and time
                %trend to construct IRFs
                
                if c_case==2
                NewX(1,1:2)=[1,time+tt];
                NewW(1,1:2)=[1,time+tt];
                elseif c_case==1
                NewX(1,1)=1;
                NewW(1,1)=1;                   
                end
                
                if nlags==6
                    NewX(1,1+c_case+5*nvars:c_case+6*nvars)=X(1,1+c_case+4*nvars:c_case+5*nvars);
                    NewW(1,1+c_case+5*nvars:c_case+6*nvars)=W(1,1+c_case+4*nvars:c_case+5*nvars);                    
                    NewX(1,1+c_case+4*nvars:c_case+5*nvars)=X(1,1+c_case+3*nvars:c_case+4*nvars);
                    NewW(1,1+c_case+4*nvars:c_case+5*nvars)=W(1,1+c_case+3*nvars:c_case+4*nvars);
                    NewX(1,1+c_case+3*nvars:c_case+4*nvars)=X(1,1+c_case+2*nvars:c_case+3*nvars);
                    NewW(1,1+c_case+3*nvars:c_case+4*nvars)=W(1,1+c_case+2*nvars:c_case+3*nvars);
                    NewX(1,1+c_case+2*nvars:c_case+3*nvars)=X(1,1+c_case+nvars:c_case+2*nvars);
                    NewW(1,1+c_case+2*nvars:c_case+3*nvars)=W(1,1+c_case+nvars:c_case+2*nvars);
                    NewX(1,1+c_case+nvars:c_case+2*nvars)=X(1,1+c_case:c_case+nvars);
                    NewW(1,1+c_case+nvars:c_case+2*nvars)=W(1,1+c_case:c_case+nvars);
                    NewX(1,1+c_case:c_case+nvars)=Y(tt,:) ;
                    NewW(1,1+c_case:c_case+nvars)=Z(tt,:);                      
                elseif nlags==5
                    NewX(1,1+c_case+4*nvars:c_case+5*nvars)=X(1,1+c_case+3*nvars:c_case+4*nvars);
                    NewW(1,1+c_case+4*nvars:c_case+5*nvars)=W(1,1+c_case+3*nvars:c_case+4*nvars);
                    NewX(1,1+c_case+3*nvars:c_case+4*nvars)=X(1,1+c_case+2*nvars:c_case+3*nvars);
                    NewW(1,1+c_case+3*nvars:c_case+4*nvars)=W(1,1+c_case+2*nvars:c_case+3*nvars);
                    NewX(1,1+c_case+2*nvars:c_case+3*nvars)=X(1,1+c_case+nvars:c_case+2*nvars);
                    NewW(1,1+c_case+2*nvars:c_case+3*nvars)=W(1,1+c_case+nvars:c_case+2*nvars);
                    NewX(1,1+c_case+nvars:c_case+2*nvars)=X(1,1+c_case:c_case+nvars);
                    NewW(1,1+c_case+nvars:c_case+2*nvars)=W(1,1+c_case:c_case+nvars);
                    NewX(1,1+c_case:c_case+nvars)=Y(tt,:) ;
                    NewW(1,1+c_case:c_case+nvars)=Z(tt,:);             
                elseif nlags==4
                    NewX(1,1+c_case+3*nvars:c_case+4*nvars)=X(1,1+c_case+2*nvars:c_case+3*nvars);
                    NewW(1,1+c_case+3*nvars:c_case+4*nvars)=W(1,1+c_case+2*nvars:c_case+3*nvars);
                    NewX(1,1+c_case+2*nvars:c_case+3*nvars)=X(1,1+c_case+nvars:c_case+2*nvars);
                    NewW(1,1+c_case+2*nvars:c_case+3*nvars)=W(1,1+c_case+nvars:c_case+2*nvars);
                    NewX(1,1+c_case+nvars:c_case+2*nvars)=X(1,1+c_case:c_case+nvars);
                    NewW(1,1+c_case+nvars:c_case+2*nvars)=W(1,1+c_case:c_case+nvars);
                    NewX(1,1+c_case:c_case+nvars)=Y(tt,:) ;
                    NewW(1,1+c_case:c_case+nvars)=Z(tt,:);
                elseif nlags==3
                    NewX(1,1+c_case+2*nvars:c_case+3*nvars)=X(1,1+c_case+nvars:c_case+2*nvars);
                    NewW(1,1+c_case+2*nvars:c_case+3*nvars)=W(1,1+c_case+nvars:c_case+2*nvars);
                    NewX(1,1+c_case+nvars:c_case+2*nvars)=X(1,1+c_case:c_case+nvars);
                    NewW(1,1+c_case+nvars:c_case+2*nvars)=W(1,1+c_case:c_case+nvars);
                    NewX(1,1+c_case:c_case+nvars)=Y(tt,:) ;
                    NewW(1,1+c_case:c_case+nvars)=Z(tt,:);
                elseif nlags==2
                    NewX(1,1+c_case+nvars:c_case+2*nvars)=X(1,1+c_case:c_case+nvars);
                    NewW(1,1+c_case+nvars:c_case+2*nvars)=W(1,1+c_case:c_case+nvars);
                    NewX(1,1+c_case:c_case+nvars)=Y(tt,:) ;
                    NewW(1,1+c_case:c_case+nvars)=Z(tt,:);
                elseif nlags==1
                    NewX(1,1+c_case:c_case+nvars)=Y(tt,:) ;
                    NewW(1,1+c_case:c_case+nvars)=Z(tt,:);
                end
                
                % the following code -- activating if ordend2 is different from 0 -- is fundamental to compute
                %GIRF for the case of endogenous interaction variables, given that in this
                %case, the value of the interaction variables change over time
                %with the evolution of the endogenous variables. 
                
                if ordend2~=0
                    
                    if nlags==1
                        ordint=c_case+ordend;
                        ordint2=c_case+ordend2;
                        NewX(1,ntotcoeff)=NewX(1,ordint2)*NewX(1,ordint);
                        NewW(1,ntotcoeff)=NewW(1,ordint2)*NewW(1,ordint);
                    elseif nlags==2
                        ordint(1,:)=[c_case+ordend,c_case+nvars+ordend]; %order of the betas associated with the endogenous variables that I want to interact with the exogenous one
                        ordint2(1,:)=[c_case+ordend2,c_case+nvars+ordend2];
                NewX(1,ntotcoeff-1)=NewX(1,ordint2(1,1))*NewX(1,ordint(1,1));
                NewX(1,ntotcoeff)=NewX(1,ordint2(1,2))*NewX(1,ordint(1,2));
                NewW(1,ntotcoeff-1)=NewW(1,ordint2(1,1))*NewW(1,ordint(1,1));
                NewW(1,ntotcoeff)=NewW(1,ordint2(1,2))*NewW(1,ordint(1,2));
                    elseif nlags==3
                        ordint(1,:)=[c_case+ordend,c_case+nvars+ordend,c_case+2*nvars+ordend]; %order of the betas associated with the endogenous variables that I want to interact with the exogenous one
                        ordint2(1,:)=[c_case+ordend2,c_case+nvars+ordend2,c_case+2*nvars+ordend2];
                NewX(1,ntotcoeff-2)=NewX(1,ordint2(1,1))*NewX(1,ordint(1,1));
                NewX(1,ntotcoeff-1)=NewX(1,ordint2(1,2))*NewX(1,ordint(1,2));
                NewX(1,ntotcoeff)=NewX(1,ordint2(1,3))*NewX(1,ordint(1,3));
                NewW(1,ntotcoeff-2)=NewW(1,ordint2(1,1))*NewW(1,ordint(1,1));
                NewW(1,ntotcoeff-1)=NewW(1,ordint2(1,2))*NewW(1,ordint(1,2));
                NewW(1,ntotcoeff)=NewW(1,ordint2(1,3))*NewW(1,ordint(1,3));
                    elseif nlags==4
                ordint(1,:)=[c_case+ordend,c_case+nvars+ordend,c_case+2*nvars+ordend,c_case+3*nvars+ordend]; %order of the betas associated with the endogenous variables that I want to interact with the exogenous one
                ordint2(1,:)=[c_case+ordend2,c_case+nvars+ordend2,c_case+2*nvars+ordend2,c_case+3*nvars+ordend2];
                 
                NewX(1,ntotcoeff-3)=NewX(1,ordint2(1,1))*NewX(1,ordint(1,1));
                NewX(1,ntotcoeff-2)=NewX(1,ordint2(1,2))*NewX(1,ordint(1,2));
                NewX(1,ntotcoeff-1)=NewX(1,ordint2(1,3))*NewX(1,ordint(1,3));
                NewX(1,ntotcoeff)=NewX(1,ordint2(1,4))*NewX(1,ordint(1,4));
                NewW(1,ntotcoeff-3)=NewW(1,ordint2(1,1))*NewW(1,ordint(1,1));
                NewW(1,ntotcoeff-2)=NewW(1,ordint2(1,2))*NewW(1,ordint(1,2));
                NewW(1,ntotcoeff-1)=NewW(1,ordint2(1,3))*NewW(1,ordint(1,3));
                NewW(1,ntotcoeff)=NewW(1,ordint2(1,4))*NewW(1,ordint(1,4));
                
                    elseif nlags==5
                ordint(1,:)=[c_case+ordend,c_case+nvars+ordend,c_case+2*nvars+ordend,c_case+3*nvars+ordend,c_case+4*nvars+ordend]; %order of the betas associated with the endogenous variables that I want to interact with the exogenous one
                ordint2(1,:)=[c_case+ordend2,c_case+nvars+ordend2,c_case+2*nvars+ordend2,c_case+3*nvars+ordend2,c_case+4*nvars+ordend2];

                NewX(1,ntotcoeff-4)=NewX(1,ordint2(1,1))*NewX(1,ordint(1,1));                
                NewX(1,ntotcoeff-3)=NewX(1,ordint2(1,2))*NewX(1,ordint(1,2));
                NewX(1,ntotcoeff-2)=NewX(1,ordint2(1,3))*NewX(1,ordint(1,3));
                NewX(1,ntotcoeff-1)=NewX(1,ordint2(1,4))*NewX(1,ordint(1,4));
                NewX(1,ntotcoeff)=NewX(1,ordint2(1,5))*NewX(1,ordint(1,5));
                
                NewW(1,ntotcoeff-4)=NewW(1,ordint2(1,1))*NewW(1,ordint(1,1));                
                NewW(1,ntotcoeff-3)=NewW(1,ordint2(1,2))*NewW(1,ordint(1,2));
                NewW(1,ntotcoeff-2)=NewW(1,ordint2(1,3))*NewW(1,ordint(1,3));
                NewW(1,ntotcoeff-1)=NewW(1,ordint2(1,4))*NewW(1,ordint(1,4));
                NewW(1,ntotcoeff)=NewW(1,ordint2(1,5))*NewW(1,ordint(1,5));
                
                    elseif nlags==6
                ordint(1,:)=[c_case+ordend,c_case+nvars+ordend,c_case+2*nvars+ordend,c_case+3*nvars+ordend,c_case+4*nvars+ordend,c_case+5*nvars+ordend]; %order of the betas associated with the endogenous variables that I want to interact with the exogenous one
                ordint2(1,:)=[c_case+ordend2,c_case+nvars+ordend2,c_case+2*nvars+ordend2,c_case+3*nvars+ordend2,c_case+4*nvars+ordend2,c_case+5*nvars+ordend2];

                NewX(1,ntotcoeff-5)=NewX(1,ordint2(1,1))*NewX(1,ordint(1,1));   
                NewX(1,ntotcoeff-4)=NewX(1,ordint2(1,2))*NewX(1,ordint(1,2));                
                NewX(1,ntotcoeff-3)=NewX(1,ordint2(1,3))*NewX(1,ordint(1,3));
                NewX(1,ntotcoeff-2)=NewX(1,ordint2(1,4))*NewX(1,ordint(1,4));
                NewX(1,ntotcoeff-1)=NewX(1,ordint2(1,5))*NewX(1,ordint(1,5));
                NewX(1,ntotcoeff)=NewX(1,ordint2(1,6))*NewX(1,ordint(1,6));
                
                NewW(1,ntotcoeff-5)=NewW(1,ordint2(1,1))*NewW(1,ordint(1,1)); 
                NewW(1,ntotcoeff-4)=NewW(1,ordint2(1,2))*NewW(1,ordint(1,2));                
                NewW(1,ntotcoeff-3)=NewW(1,ordint2(1,3))*NewW(1,ordint(1,3));
                NewW(1,ntotcoeff-2)=NewW(1,ordint2(1,4))*NewW(1,ordint(1,4));
                NewW(1,ntotcoeff-1)=NewW(1,ordint2(1,5))*NewW(1,ordint(1,5));
                NewW(1,ntotcoeff)=NewW(1,ordint2(1,6))*NewW(1,ordint(1,6));
                        
                    else disp('Error: for now it is possible only to set nlags up to 6. But if you need it is very simple to modify for that');
                    end
                    
                end
                
                X=NewX;
                W=NewW;
               
                epsilon(1,:)=(u(1+tt,:)')' ; % to consider the reduced form residuals
                
           
                Y(tt+1,:)=((X*beta+epsilon)')'; %the extracted residuals are used for both cases
                Z(tt+1,:)=((W*beta+epsilon)')';
                
                IRF(tt+1,:,jj,ic)=Y(tt+1,:)-Z(tt+1,:); %store the empirical response fot the period tt
                
                
            end
            
            %New part: to repeat draw if residuals estracted make the GIRF explosive or NaN.
                %The condition is imposed to the last step ahead (looking to the shocked variable)        
            if IRF(nsteps,jj,jj,ic)<=S(pick,pick)*10*abs(init_shock) && IRF(nsteps,jj,jj,ic)>=-S(pick,pick)*10*abs(init_shock); % the GIRF is consided non-explosive if at the last step ahead considered the shocked variable takes values smaller than extreme value, such as 10 standard deviations of the structural shock considered
                % If this is not the case the GIRFs is considered explosive
                % and hence the draw is repeated since it might depend on
                % the particular sequence of future shocks considered. In
                % case of explosiveness (according to the definition
                % given), the GIRF will have no economic interpretation.
                 % Admittedly, the criterion is somehow arbitrary, but there is no easy
                 % analytical condition to establish the stability of a
                 % polinomially nonlinear VAR model. The user has though to make sure
                 % that reasonable changes in the criterion do not affect
                 % results. The code allows to take account of discarded
                 % responses.
                   % [Notice: the DSGE model literature has developped the concept of Pruning to ensure stable GIRFs. 
                   % We plan to work in that direction for future versions]
                
                 ic=ic+1;
                 
            else %the draw is repeated since ic remain the same...
                countRep(1,ii)=countRep(1,ii)+1;
                
                    if countRep(1,ii)==draws/2% ...but if the same extraction has been repeated for much of the draws considered it means that is not only due to a particular sequence of future shocks, but rather that the starting history considered is explosive in itself (maybe because is an extreme history or because it is a history near the final date of the sample) 
                        nexplos=nexplos+1;
                        IRF=zeros(nsteps,nvars,nvars,draws); % to not count at all when computing the average (discarding such a explosive histories is important since we are using the average (rather than the median for example))
                        timeExp=[timeExp, time]; 
                        break % to leave the while loop and consider the next history
                    end
            end
           
            %% here we have to average across (future) shocks the GIRFs
        end
       
       %IRF(:,2:5,pick,:)=cumsum(IRF(:,2:5,pick,:)); % to use when there are variables in growth rates in positions 2:5
        
   IIRF=mean(IRF,4);    

            
  OIRFhist(:,:,:,ii)=IIRF;% here we store the single-history conditional GIRF      

    
end


%% NOW WE NEED TO AVERAGE ACROSS HISTORIES
% to get a consistent estimate of our state-conditional GIRF

OIRFavg=zeros(nsteps,nvars,nvars);
temp=zeros(nsteps,nvars,nvars);
for yy=1:nhist
        temp=temp+OIRFhist(:,:,:,yy);
end
OIRFavg=temp/(nhist-nexplos); %take the sample mean across the non explosive IRF

if nexplos~=0
    if c_case>1
        disp(['to note that girf computed starting from some histories are esplosive -- msg from VARirtrue_sim1shock_histories...:' num2str(nexplos) '  / ' num2str(nhist) ' that correspond to the histories at time: ' num2str(timeExp) ] )
    else
        disp(['to note that girf computed starting from some histories are esplosive -- msg from VARirtrue_sim1shock_histories...:' num2str(nexplos) '  / ' num2str(nhist) ] )
    end
end

OIRF_optavg.S=S;
if proxyvar == 1
    OIRF_optavg.s = s;
    OIRF_optavg.beta1 = beta1;
    OIRF_optavg.Fstat = Fstat;
end
OIRF_optavg.nsteps    = nsteps;
OIRF_optavg.ordend=ordend; 
OIRF_optavg.ordend2=ordend2;
OIRF_optavg.pick=pick;
OIRF_optavg.chol_flag = chol_flag;
OIRF_optavg.nexplos=nexplos;
OIRF_optavg.countRep=countRep;
OIRF_optavg.mode=mode;
OIRF_optavg.draws=draws;
