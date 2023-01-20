
%% Load full data
mkdir(spec)
load data/data_m_full.mat


%% column subsetting
numend = length(labels_sel)-1;
data = ones(size(data_lib, 1), length(labels_sel)+12-1)*NaN;
labels_print = labels_sel(1:numend);

for i = 1:length(labels_sel)
    if i <= numend
        col = find(strcmp(labels_lib, labels_sel{i}));
        data(:,i) = data_lib(:,col);
        labels_print{i} = labels_print_lib{col};
    else
        col = find(strcmp(labels_lib, strcat(labels_sel{i},'(-1)')));
        data(:,i:end) = data_lib(:,col:col+(12-1));
        intname = labels_sel{i};
        for l = 1:12
            if l ==1
                labels_sel{i} = strcat(intname,'(-',num2str(l),')');
            else
                labels_sel{length(labels_sel) + 1} = strcat(intname,'(-',num2str(l),')');
            end
        end      
    end
end

labels = labels_sel;

% instrument
if proxyvar == 1
    col = find(strcmp(labels_lib, instr_label));
    instr = data_lib(:,col);
end

%% shifter (not all variables are available from start)
% make sure data ends 2020m1
data = data(1:481,:);
if proxyvar == 1
    instr = instr(1:481);
end
% exclude NAs
allcomplete = logical(sum(~isnan(data(:,1:numend)),2) == numend);
data = data(allcomplete, :);
if proxyvar == 1
    instr = instr(allcomplete);
end
xt00 = xt00 - min(find(allcomplete == 1)) + 1;
clear allcomplete
if xt00 > 12
    data = data(xt00 - 12:end,:);  % make sure starting data is 1999m1 and end 2020m1
    if proxyvar == 1
        instr = instr(xt00 - 12:end);
    end
    xt00 = 12;
end


%% plot state
x = 1:size(data,1);
xx = [x;x];
col = find(strcmp(labels_lib, labels{ordend2}));
y = data(:, ordend2)';
yorig = y;
if defstate.growth == 1
    y = [NaN; data(2:end, ordend2) - data(1:end-1, ordend2)]';
end
if defstate.cumul > 0
    helpstate = NaN*ones(size(y,2),defstate.cumul);
    helpstate(:,1) = y;
    for i = 2:defstate.cumul
        helpstate(:,i) = [repmat(NaN,i-1,1); y(1:(end-i+1))'];
    end
    y = sum(helpstate,2)';
end
s = y <= prctile(y, perc);
yy = [1.1*max([y, yorig])*s;
      1.1*max([y, yorig])*s];
yyneg = [min([y, yorig]-0.1*abs(min([y, yorig])))*s;
         min([y, yorig]-0.1*abs(min([y, yorig])))*s];


figure()
h1 = area(xx([2:end end]), yy(1:end), 'LineStyle', 'none', 'FaceColor', [0.9 0.9 0.9]);
hold on
h2 = area(xx([2:end end]), yyneg(1:end), 'LineStyle', 'none', 'FaceColor', [0.9 0.9 0.9]);
h3 = plot(x, y, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2);
h4 = plot(x, yorig, 'Color', 'k', 'LineWidth', 1.5);
plot(x, zeros(size(x,2)), 'k')
axis('tight')
grid on
xt = [flip(xt00:-48:1) (xt00+48):48:length(yy)];
set(gca, 'XTick',xt, 'XTickLabel',(xt-xt00)/12 + 2000)
title(labels_print{ordend2})
legend(h1, 'Pessimistic times','Location','SouthEast')
set(gca,'FontSize',16)
set(gcf, 'position', [0 0 800 400]);
set(gcf,'paperpositionmode','auto')
print(gcf,'-depsc2','-loose',strcat('./',spec,'/state'));
if strcmp(spec, 'm01_baseline')
    mkdir('paper_results')
    print(gcf,'-depsc2','-loose',strcat('./paper_results/fig1.eps'));
end

clear x y ysmooth s helpstate xx yy yyneg h1 h2 h3 xt

xt00 = xt00 - nlags;

