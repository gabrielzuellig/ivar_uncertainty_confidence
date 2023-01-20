function VARplotSDirf(IRFlin,IRFmed1,IRFmed2,reverseHL,pick,labels,INF1,SUP1,INF2,SUP2)
% =======================================================================
% plot two different GIRFS
% =======================================================================

if reverseHL==1

    % central
    tmp = IRFmed1;
    IRFmed1 = IRFmed2;
    IRFmed2 = tmp;

    %confLev
    if exist('INF1')
        tmp = INF1;
        INF1 = INF2;
        INF2 = tmp;

        %confLev2
        tmp = SUP1;
        SUP1 = SUP2;
        SUP2 = tmp;
    end

end


%% Check optional inputs & Define some parameters
%================================================
if ~exist('filename','var')
    filename = 'shock_';
end

if ~exist('FontSize','var')
    FontSize = 16;
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
for jj = pick:nshocks
    for ii=1:nvars
        subplot(row,col,ii);
        if exist('INF1','var') && exist('SUP1','var')
            plot1=shadedplot(steps,INF1(:,ii,jj)',SUP1(:,ii,jj)',[0.8 0.8 0.8],[0.8 0.8 0.8] );
            box off
            grid on
            set(gca,'GridLineStyle','-','Layer','bottom')
            if nsteps == 20
                set(gca,'XTickLabelMode', 'manual','XTickLabel',{'5','10','15','20'},'XTick',[5 10 15 20])
            elseif nsteps == 48
                set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36','48'},'XTick',[12 24 36 48])
            end
            hold on
            plot1=plot(steps,INF1(:,ii,jj),'LineStyle','-','Color',rgb('dark blue'),'LineWidth',1);
            hAnnotation = get(plot1,'Annotation');
            hLegendEntry = get(hAnnotation','LegendInformation');
            set(hLegendEntry,'IconDisplayStyle','off')
            hold on
            plot1=plot(steps,SUP1(:,ii,jj),'LineStyle','-','Color',rgb('dark blue'),'LineWidth',1);
            hAnnotation = get(plot1,'Annotation');
            hLegendEntry = get(hAnnotation','LegendInformation');
            set(hLegendEntry,'IconDisplayStyle','off')
            hold on
            plot1=plot(steps,INF2(:,ii,jj),'LineStyle','-','Color','k','LineWidth',1);
            hAnnotation = get(plot1,'Annotation');
            hLegendEntry = get(hAnnotation','LegendInformation');
            set(hLegendEntry,'IconDisplayStyle','off')
            hold on
            plot1=plot(steps,SUP2(:,ii,jj),'LineStyle','-','Color','k','LineWidth',1);
            hAnnotation = get(plot1,'Annotation');
            hLegendEntry = get(hAnnotation','LegendInformation');
            set(hLegendEntry,'IconDisplayStyle','off')
        end
        hold on
        plot1=plot(x_axis,'k','LineWidth',1);
        hAnnotation = get(plot1,'Annotation');
        hLegendEntry = get(hAnnotation','LegendInformation');
        set(hLegendEntry,'IconDisplayStyle','off')
        hold on
        p2=plot(steps,IRFmed1(:,ii,jj),'LineStyle','--','Color',rgb('dark blue'),'LineWidth',2);
        hold on
        p3=plot(steps,IRFmed2(:,ii,jj),'-k','LineWidth',3);
        xlim([1 nsteps]);
        if exist('labels','var')
            title(labels(ii),'FontSize',12);
        end
    end

end

