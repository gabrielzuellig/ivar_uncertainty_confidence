
%% Import all series
[data labels] = xlsread('data/all_series.xlsx','m');
data = data(:,2:end);
labels = labels(1,:);


%% Data treatments

% 100*log of variables in levels
var = {'Output','Manu','ECBAssets'};
for vv = 1:length(var)
    col = find(strcmp(labels, var{vv}));
    if ~isempty(col)
        data(:, col) = 100*log(data(:, col));
    end
end

% level-100
var = {'EconSent'};
for vv = 1:length(var)
    col = find(strcmp(labels, var{vv}));
    data(:, col) = data(:, col)-100;
end

% level+10
var = {'ConsConfEC','IndConfEC'};
for vv = 1:length(var)
    col = find(strcmp(labels, var{vv}));
    if strcmp(var{vv}, 'ConsConfEC')
        data(:, col) = data(:, col)+10;
    elseif strcmp(var{vv}, 'IndConfEC')
        data(:, col) = data(:, col)+5;
    end
end

% y/y growth rates 
var = {'CPI','Manu','ECBAssets'};
for vv = 1:length(var)
    col = find(strcmp(labels, var{vv}));
    if ~isempty(col)
        data(:, end+1) = NaN;
        data(13:end, end) = data(13:end, col) - data(1:end-12, col);
        labels{length(labels) + 1} = strcat(var{vv},'Gr12');
    end
end

% m/m growth rates
var = {'ConsConfEC','IndConfEC'};
for vv = 1:length(var)
    col = find(strcmp(labels, var{vv}));
    if ~isempty(col)
        data(:, end+1) = NaN;
        data(2:end, end) = data(2:end, col) - data(1:end-1, col);
        labels{length(labels) + 1} = strcat(var{vv},'Gr');
    end
end

% Orthogonalize ConsConf for VSTOXX
xx = data(:,find(strcmp(labels, 'VSTOXX')));
xx = [ones(size(xx)), xx];
yy = data(:,find(strcmp(labels, 'ConsConfEC')));
choose = logical(sum(isnan([xx, yy]), 2) == 0);
xx = xx(choose,:);
yy = yy(choose,:);
beta = inv(xx'*xx)*(xx'*yy);
yyfit = xx*beta;
figure()
plot(yy)
hold on
plot(yyfit)
yyres = yy-yyfit;
plot(yyres)
legend('ConcConfEC','Fitted values (pred. with VSTOXX)','Residual')
data(:,end+1) = NaN;
data(choose, end) = yyres;
labels{length(labels) + 1} = strcat('ConsConfECres');
clear xx yy choose beta yyfit yyres
% orthogonalize VSTOXX for ConsConf
xx = data(:,find(strcmp(labels, 'ConsConfEC')));
xx = [ones(size(xx)), xx];
yy = data(:,find(strcmp(labels, 'VSTOXX')));
choose = logical(sum(isnan([xx, yy]), 2) == 0);
xx = xx(choose,:);
yy = yy(choose,:);
beta = inv(xx'*xx)*(xx'*yy);
yyfit = xx*beta;
figure()
plot(yy)
hold on
plot(yyfit)
yyres = yy-yyfit;
plot(yyres)
legend('VSTOXX','Fitted values (pred. with ConsConf)','Residual')
data(:,end+1) = NaN;
data(choose, end) = yyres;
labels{length(labels) + 1} = strcat('VSTOXXres');
clear xx yy choose beta yyfit yyres


%% Cosmetics
% "anchoring" time at 2000m1
xt00 = 241;
% get rid of unnecessary variables (year, month indices)
% 'nice' labels for printing
labels_print = labels;
labels_print(find(strcmp(labels,'Output'))) = {'Output'};
labels_print(find(strcmp(labels,'Manu'))) = {'Manuf. output'};
labels_print(find(strcmp(labels,'ManuGr12'))) = {'Manuf. output (y/y)'};
labels_print(find(strcmp(labels,'VSTOXX'))) = {'Uncertainty (VSTOXX)'};
labels_print(find(strcmp(labels,'VIX'))) = {'Uncertainty (VIX)'};
labels_print(find(strcmp(labels,'VSTOXXres'))) = {'Uncertainty (VSTOXX)'};
labels_print(find(strcmp(labels,'ConsConfEC'))) = {'Consumer confidence'};
labels_print(find(strcmp(labels,'ConsConfECres'))) = {'Consumer confidence'};
labels_print(find(strcmp(labels,'IndConfEC'))) = {'Business confidence'};
labels_print(find(strcmp(labels,'EconSent'))) = {'Economic sentiment'};
labels_print(find(strcmp(labels,'CPI'))) = {'Consumer prices'};
labels_print(find(strcmp(labels,'CPIGr12'))) = {'Inflation'};
labels_print(find(strcmp(labels,'NFCSprd'))) = {'Bond spread'};
labels_print(find(strcmp(labels,'EONIA'))) = {'Money market rate'};
labels_print(find(strcmp(labels,'GoldProxy'))) = {'Gold price instrument (Piffer & Podstawski)'};
labels_print(find(strcmp(labels,'JLNFin'))) = {'Uncertainty (JLN)'};
labels_print(find(strcmp(labels,'ConsDisp'))) = {'Uncertainty (Consumer disagreement)'};
labels_print(find(strcmp(labels,'IndDisp'))) = {'Uncertainty (Industry disagreement)'};
labels_print(find(strcmp(labels,'ECBAssets'))) = {'ECB balance sheet'};
labels_print(find(strcmp(labels,'CommCarGr12'))) = {'Commercial car registrations'};
labels_print(find(strcmp(labels,'ShadowWX'))) = {'Shadow rate (Wu & Xia)'};
labels_print(find(strcmp(labels,'ShadowDR'))) = {'Shadow rate (De Rezende)'};


%% Interactions and lags of interactions
for i = 1:length(uctyvars)
    for j = 1:length(statevars)
        col1 = find(strcmp(labels, uctyvars{i}));
        int = ones(size(data,1), 12+1)*NaN;
        % interaction with state variable
        col2 = find(strcmp(labels, statevars{j}));
        if ~isempty(col1) && ~isempty(col2) 
            int(:, 1) = data(:,col1) .* data(:,col2);
            % lags
            for l = 1:12
                int(l+1:end, l+1) = int(1:end-l, 1);
                labels{length(labels) + 1} = strcat(uctyvars{i}, 'x',statevars{j},'(-',num2str(l),')');
            end
            % put together
            data = [data, int(:,2:end)];
        end
    end
end


%% Export, housekeeping
data_lib = data;
labels_lib = labels;
labels_print_lib = labels_print;
clear col col1 col2 int vv i j l var data labels labels_print m x y
save('data/data_m_full.mat', 'data_lib', 'labels_lib', 'labels_print_lib','xt00')
close all
