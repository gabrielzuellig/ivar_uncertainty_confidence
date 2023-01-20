function VARplotdiff(IRFmed1,IRFmed2,pick,labels,INF1,SUP1,INF2,SUP2,leg,plPE)
% =======================================================================
% Plot the difference between two GIRFs
% =======================================================================

%% Check optional inputs & Define some parameters
%================================================
if ~exist('filename','var')
    filename = 'shock_';
end

if ~exist('leg','var')
    leg = 0;
end

if ~exist('plPE','var')
    plPE = 1;
end


% Initialize IRF matrix
[nsteps, nvars, nshocks] = size(IRFmed1);

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

for ii=1:nvars
    subplot(row,col,ii);
    if exist('INF1','var') && exist('SUP1','var')

        p2=shadedplot(steps,INF2(:,ii)',SUP2(:,ii)',[0.8 0.8 0.8],[0.8 0.8 0.8] );
        hold on
        p1=shadedplot(steps,INF1(:,ii)',SUP1(:,ii)',[0.6 0.6 0.6],[0.6 0.6 0.6] );
        hold on
    end
    box off
    grid on
    set(gca,'GridLineStyle','-','Layer','bottom')

    if plPE==1

        plot(steps,IRFmed1(:,ii),'LineStyle','-','Color',rgb('black'),'LineWidth',2);
        hold on
        plot(steps,IRFmed2(:,ii),'LineStyle','-','Color',rgb('black'),'LineWidth',2);
    end
    plot(x_axis,'k','LineWidth',1)
    xlim([1 nsteps]);
    if exist('labels','var')
        title(labels(ii),'FontSize',12);
    end
    if nsteps == 20
        set(gca,'XTickLabelMode', 'manual','XTickLabel',{'5','10','15','20'},'XTick',[5 10 15 20])
    elseif nsteps==48
        set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36','48'},'XTick',[12 24 36 48])
    end
    if ii==1 && leg==1
        legend([p1(2) p2(2)],'68% conf. bands','90% conf. bands','Location','NorthWest')
    end
end

