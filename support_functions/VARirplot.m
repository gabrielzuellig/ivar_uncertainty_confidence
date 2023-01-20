function VARirplot(IRF,pick,labels,INF,SUP)
% =======================================================================
% Plot the IRFs computed with VARir
% =======================================================================
% VARirplot(IRF,pick,labels,filename,INF,SUP)
% -----------------------------------------------------------------------
% INPUT
%   IRF(:,:,:) : matrix with periods, variable, shock
%
% OPTIONAL INPUT
%   pick     : a vector containing the shocks to be plotted [default 0 = all]
%   labels   : name of the variables, ordered as in the Cholesky
%   filename : name for file saving
%   INF      : lower error band
%   SUP      : upper error band
% =======================================================================
% Ambrogio Cesa Bianchi, May 2012
% ambrogio.cesabianchi@gmail.com


%% Check optional inputs & Define some parameters
%================================================
if ~exist('filename','var') 
    filename = 'shock_';
end

if ~exist('FontSize','var') 
    FontSize = 16;
end

% Initialize IRF matrix
[nsteps, nvars, nshocks] = size(IRF);

% If one shock is chosen, set the right value for nshocks
if ~exist('pick','var') 
    pick = 1;
else
    if pick<0 || pick>nvars
        error('The selected shock is non valid')
    else
        if pick==0
            pick=1;
        else
            nshocks = pick;
        end
    end
end

% Define the rows and columns for the subplots
row = round(sqrt(nvars));
col = ceil(sqrt(nvars));

% Define a timeline
steps = 1:1:nsteps;
x_axis = zeros(1,nsteps);


%% Plot
%=========
%FigSize(24,8)
for jj = pick:nshocks                
    for ii=1:nvars
        subplot(row,col,ii);
        plot(steps,IRF(:,ii,jj),'LineStyle','-','Color','k','LineWidth',2);
        hold on
        plot(x_axis,'k','LineWidth',1)
        if exist('INF','var') && exist('SUP','var')
            plot(steps,INF(:,ii,jj),'LineStyle','--','Color','k','LineWidth',1);
            hold on
            plot(steps,SUP(:,ii,jj),'LineStyle','--','Color','k','LineWidth',1);
        end
        xlim([1 nsteps]);
        box off 
        grid on
        if exist('labels','var') 
            title(labels(ii),'FontSize',12); 
        end
%         xlabel('Place your label here');
%         ylabel('Place your label here');
        %FigFont(FontSize);
    end

%     % Save
%     set(gcf, 'Color', 'w');
%     FigName = [filename num2str(jj)];
%     export_fig(FigName,'-pdf','-png','-painters')
%     clf('reset');
end

% close all
