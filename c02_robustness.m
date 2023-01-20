


%% Model 2: VIX
ordend = 1;
ordend2 = 4;
spec = strcat('m02_vix');
labels_sel =  {'VIX','ManuGr12','CPIGr12','ConsConfEC','EONIA','VIXxConsConfEC'};

r01_subset_data;

r02_ivar;


%% Model 3: LMN financial uncertainty
ordend = 1;
ordend2 = 4;
spec = strcat('m03_lmn');
labels_sel =  {'JLNFin','ManuGr12','CPIGr12','ConsConfEC','EONIA','JLNFinxConsConfEC'};

r01_subset_data;

r02_ivar;


%% Model 4: Consumer outlook disagreement
ordend = 1;
ordend2 = 4;
spec = strcat('m04_consdisp');
labels_sel =  {'ConsDisp','ManuGr12','CPIGr12','ConsConfEC','EONIA','ConsDispxConsConfEC'};

r01_subset_data;

r02_ivar;


%% Model 5: Industry outlook disagreement
ordend = 1;
ordend2 = 4;
spec = strcat('m05_inddisp');
labels_sel =  {'IndDisp','ManuGr12','CPIGr12','ConsConfEC','EONIA','IndDispxConsConfEC'};

r01_subset_data;

r02_ivar;


%% Model 6: Economic sentiment
ordend = 1;
ordend2 = 4;
spec = strcat('m06_econsent');
labels_sel =  {'VSTOXX','ManuGr12','CPIGr12','EconSent','EONIA','VSTOXXxEconSent'};

r01_subset_data;

r02_ivar;


%% Model 7: Industry confidence
ordend = 1;
ordend2 = 4;

spec = strcat('m07_indconf');
labels_sel =  {'VSTOXX','ManuGr12','CPIGr12','IndConfEC','EONIA','VSTOXXxIndConfEC'};

r01_subset_data;

r02_ivar;


%% Model 8: Commercial cars with consumer confidence
nlags = 3;  
ordend = 1;
ordend2 = 4;
spec = strcat('m08_commcar_consconf');
labels_sel =  {'VSTOXX','CommCarGr12','CPIGr12','ConsConfEC','EONIA','VSTOXXxConsConfEC'};

r01_subset_data;

r02_ivar;



%% Model 9: Commercial cars with industry confidence
nlags = 3;  
ordend = 1;
ordend2 = 4;
spec = strcat('m09_commcar_indconf');
labels_sel =  {'VSTOXX','CommCarGr12','CPIGr12','IndConfEC','EONIA','VSTOXXxIndConfEC'};

r01_subset_data;

r02_ivar;


%% Model 10: Mood swings
nlags = 4;  
ordend = 1;
ordend2 = 4;
spec = strcat('m10_moodswings');
labels_sel =  {'VSTOXX','ManuGr12','CPIGr12','ConsConfEC','EONIA','VSTOXXxConsConfEC'};
defstate.growth = 1; % 1 for growth, 0 for level of state
defstate.cumul = 6;

r01_subset_data;

r02_ivar;


%% Figures 4 + 5: Load data

load('./m01_baseline/out.mat')
ipL.inf = INFavgL2(:,2,1);
ipL.sup = SUPavgL2(:,2,1);
ipH.inf = INFavgH2(:,2,1);
ipH.sup = SUPavgH2(:,2,1);
piL.inf = INFavgL2(:,3,1);
piL.sup = SUPavgL2(:,3,1);
piH.inf = INFavgH2(:,3,1);
piH.sup = SUPavgH2(:,3,1);
ipL.m1 = OIRFavg(:,2,1);
ipH.m1 = OIRF2avg(:,2,1);
piL.m1 = OIRFavg(:,3,1);
piH.m1 = OIRF2avg(:,3,1);
load('./m02_vix/out.mat')
ipL.m2 = OIRFavg(:,2,1);
ipH.m2 = OIRF2avg(:,2,1);
piL.m2 = OIRFavg(:,3,1);
piH.m2 = OIRF2avg(:,3,1);
load('./m03_lmn/out.mat')
ipL.m3 = OIRFavg(:,2,1);
ipH.m3 = OIRF2avg(:,2,1);
piL.m3 = OIRFavg(:,3,1);
piH.m3 = OIRF2avg(:,3,1);
load('./m04_consdisp/out.mat')
ipL.m4 = OIRFavg(:,2,1);
ipH.m4 = OIRF2avg(:,2,1);
piL.m4 = OIRFavg(:,3,1);
piH.m4 = OIRF2avg(:,3,1);
load('./m05_inddisp/out.mat')
ipL.m5 = OIRFavg(:,2,1);
ipH.m5 = OIRF2avg(:,2,1);
piL.m5 = OIRFavg(:,3,1);
piH.m5 = OIRF2avg(:,3,1);
load('./m06_econsent/out.mat')
ipL.m6 = OIRFavg(:,2,1);
ipH.m6 = OIRF2avg(:,2,1);
piL.m6 = OIRFavg(:,3,1);
piH.m6 = OIRF2avg(:,3,1);
load('./m07_indconf/out.mat')
ipL.m7 = OIRFavg(:,2,1);
ipH.m7 = OIRF2avg(:,2,1);
piL.m7 = OIRFavg(:,3,1);
piH.m7 = OIRF2avg(:,3,1);
load('./m10_moodswings/out.mat')
ipL.m8 = OIRFavg(:,2,1);
ipH.m8 = OIRF2avg(:,2,1);
piL.m8 = OIRFavg(:,3,1);
piH.m8 = OIRF2avg(:,3,1);
load('./m08_commcar_consconf/out.mat')
ipL.m9 = OIRFavg(:,2,1);
ipH.m9 = OIRF2avg(:,2,1);
piL.m9 = OIRFavg(:,3,1);
piH.m9 = OIRF2avg(:,3,1);
load('./m09_commcar_indconf/out.mat')
ipL.m10 = OIRFavg(:,2,1);
ipH.m10 = OIRF2avg(:,2,1);
piL.m10 = OIRFavg(:,3,1);
piH.m10 = OIRF2avg(:,3,1);


%% Figure 4: Robustness for industrial production (4x2 panel)

figure()

subplot(4,2,1)
plot1 = shadedplot(1:36,ipH.inf',ipH.sup',[220, 207, 244] / 244,[220, 207, 244] / 244);
hold on
% plot(1:36, ipH.m1, 'k', 'LineWidth', 2)
plot(1:36,ipH.m1,'LineStyle','--','Color',rgb('dark blue'),'LineWidth',2);
f2 = plot(1:36, ipH.m2, '--r', 'LineWidth', 2);
f3 = plot(1:36, ipH.m3, '-b*');
axis('tight')
ylim([-1.8 0.7])
box off
grid on
set(gca,'GridLineStyle','-','Layer','bottom')
set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36'},'XTick',[12 24 36 48])
title({'Industrial production:','Normal times'},'FontSize',14)
legend([f2 f3], {'VIX', 'Fin. uncertainty (LMN)'},'Location','SouthWest','FontSize',10)
legend boxoff

subplot(4,2,2)
plot1=plot(1:36,ipL.inf','LineStyle','-','Color','k','LineWidth',1);hold on
plot1=plot(1:36,ipL.sup','LineStyle','-','Color','k','LineWidth',1);
plot(1:36, ipL.m1, 'k', 'LineWidth', 2)
plot(1:36, ipL.m2, '--r', 'LineWidth', 2);
plot(1:36, ipL.m3, '-b*');
axis('tight')
ylim([-1.8 0.7])
box off
grid on
set(gca,'GridLineStyle','-','Layer','bottom')
set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36'},'XTick',[12 24 36 48])
title({'Industrial production:','Pessimistic times'},'FontSize',14)

subplot(4,2,3)
plot1 = shadedplot(1:36,ipH.inf',ipH.sup',[220, 207, 244] / 244,[220, 207, 244] / 244 );
hold on
plot(1:36,ipH.m1,'LineStyle','--','Color',rgb('dark blue'),'LineWidth',2);
f4 = plot(1:36, ipH.m4, '--r', 'LineWidth', 2);
f5 = plot(1:36, ipH.m5, '-b*');
axis('tight')
ylim([-1.8 0.7])
box off
grid on
set(gca,'GridLineStyle','-','Layer','bottom')
set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36'},'XTick',[12 24 36 48])
legend([f4 f5], {'Consumer outlook disagreement',...
    'Industry outlook disagreement',},'Location','SouthWest','FontSize',10)
legend boxoff

subplot(4,2,4)
plot1=plot(1:36,ipL.inf','LineStyle','-','Color','k','LineWidth',1);hold on
plot1=plot(1:36,ipL.sup','LineStyle','-','Color','k','LineWidth',1);
plot(1:36, ipL.m1, 'k', 'LineWidth', 2)
plot(1:36, ipL.m4, '--r', 'LineWidth', 2);
plot(1:36, ipL.m5, '-b*');
axis('tight')
ylim([-1.8 0.7])
box off
grid on
set(gca,'GridLineStyle','-','Layer','bottom')
set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36'},'XTick',[12 24 36 48])

subplot(4,2,5)
plot1 = shadedplot(1:36,ipH.inf',ipH.sup',[220, 207, 244] / 244,[220, 207, 244] / 244 );
hold on
plot(1:36,ipH.m1,'LineStyle','--','Color',rgb('dark blue'),'LineWidth',2);
f6 = plot(1:36, ipH.m6, '--r', 'LineWidth', 2);
f7 = plot(1:36, ipH.m7, '-b*');
f8 = plot(1:36, ipH.m8, '--mo');
axis('tight')
ylim([-1.8 0.7])
box off
grid on
set(gca,'GridLineStyle','-','Layer','bottom')
set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36'},'XTick',[12 24 36 48])
legend([f6 f7 f8], {'Economic sentiment', 'Industry confidence','Consumer mood swings'},'Location','SouthWest','FontSize',10)
legend boxoff

subplot(4,2,6)
plot1=plot(1:36,ipL.inf','LineStyle','-','Color','k','LineWidth',1);hold on
plot1=plot(1:36,ipL.sup','LineStyle','-','Color','k','LineWidth',1);
plot(1:36, ipL.m1, 'k', 'LineWidth', 2)
plot(1:36, ipL.m6, '--r', 'LineWidth', 2);
plot(1:36, ipL.m7, '-b*');
plot(1:36, ipL.m8, '--mo');
axis('tight')
ylim([-1.8 0.7])
box off
grid on
set(gca,'GridLineStyle','-','Layer','bottom')
set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36'},'XTick',[12 24 36 48])

subplot(4,2,7)
plot1 = shadedplot(1:36,ipH.inf',ipH.sup',[220, 207, 244] / 244,[220, 207, 244] / 244 );
hold on
plot(1:36,ipH.m1,'LineStyle','--','Color',rgb('dark blue'),'LineWidth',2);
f9 = plot(1:36, ipH.m9, '--r', 'LineWidth', 2);
f10 = plot(1:36, ipH.m10, '-b*');
axis('tight')
ylim([-2 1])
box off
grid on
set(gca,'GridLineStyle','-','Layer','bottom')
set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36'},'XTick',[12 24 36 48])
legend([f9 f10], {'Consumer confidence','Industry confidence'},'Location','SouthWest','FontSize',10)
legend boxoff
title({'Commercial car registrations'},'FontSize',14)

subplot(4,2,8)
plot1=plot(1:36,ipL.inf','LineStyle','-','Color','k','LineWidth',1);hold on
plot1=plot(1:36,ipL.sup','LineStyle','-','Color','k','LineWidth',1);
plot(1:36, ipL.m1, 'k', 'LineWidth', 2)
plot(1:36, ipL.m9, '--r', 'LineWidth', 2);
plot(1:36, ipL.m10, '-b*');
axis('tight')
ylim([-2 1])
box off
grid on
set(gca,'GridLineStyle','-','Layer','bottom')
set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36'},'XTick',[12 24 36 48])
title({'Commercial car registrations'},'FontSize',14)

set(gcf, 'position', [0 0 600 1100]);
set(gcf,'paperpositionmode','auto')
print(gcf,'-depsc2','-loose','paper_results/fig4.eps');


%% Figure 5: Robustness for inflation (4x2 panel)

figure()

subplot(4,2,1)
plot1 = shadedplot(1:36,piH.inf',piH.sup',[220, 207, 244] / 244,[220, 207, 244] / 244 );
hold on
plot(1:36,piH.m1,'LineStyle','--','Color',rgb('dark blue'),'LineWidth',2);
f2 = plot(1:36, piH.m2, '--r', 'LineWidth', 2);
f3 = plot(1:36, piH.m3, '-b*');
axis('tight')
ylim([-0.2 0.08])
box off
grid on
set(gca,'GridLineStyle','-','Layer','bottom')
set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36'},'XTick',[12 24 36 48])
title({'Inflation:','Normal times'},'FontSize',14)
legend([f2 f3], {'VIX', 'Fin. uncertainty (LMN)'},'Location','SouthWest','FontSize',10)
legend boxoff

subplot(4,2,2)
plot1=plot(1:36,piL.inf','LineStyle','-','Color','k','LineWidth',1);hold on
plot1=plot(1:36,piL.sup'','LineStyle','-','Color','k','LineWidth',1);
plot(1:36, piL.m1, 'k', 'LineWidth', 2)
plot(1:36, piL.m2, '--r', 'LineWidth', 2)
plot(1:36, piL.m3, '-b*')
axis('tight')
ylim([-0.2 0.08])
box off
grid on
set(gca,'GridLineStyle','-','Layer','bottom')
set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36'},'XTick',[12 24 36 48])
title({'Inflation:','Pessimistic times'},'FontSize',14)

subplot(4,2,3)
plot1 = shadedplot(1:36,piH.inf',piH.sup',[220, 207, 244] / 244,[220, 207, 244] / 244 );
hold on
plot(1:36,piH.m1,'LineStyle','--','Color',rgb('dark blue'),'LineWidth',2);
f4 = plot(1:36, piH.m4, '--r', 'LineWidth', 2);
f5 = plot(1:36, piH.m5, '-b*');
axis('tight')
ylim([-0.2 0.08])
box off
grid on
set(gca,'GridLineStyle','-','Layer','bottom')
set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36'},'XTick',[12 24 36 48])
legend([f4 f5], {'Consumer outlook disagreement',...
    'Industry outlook disagreement'},'Location','SouthWest','FontSize',10)
legend boxoff

subplot(4,2,4)
plot1=plot(1:36,piL.inf','LineStyle','-','Color','k','LineWidth',1);hold on
plot1=plot(1:36,piL.sup'','LineStyle','-','Color','k','LineWidth',1);
hold on
plot(1:36, piL.m1, 'k', 'LineWidth', 2)
plot(1:36, piL.m4, '--r', 'LineWidth', 2)
plot(1:36, piL.m5, '-b*')
axis('tight')
ylim([-0.2 0.08])
box off
grid on
set(gca,'GridLineStyle','-','Layer','bottom')
set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36'},'XTick',[12 24 36 48])

subplot(4,2,5)
plot1 = shadedplot(1:36,piH.inf',piH.sup',[220, 207, 244] / 244,[220, 207, 244] / 244 );
hold on
plot(1:36,piH.m1,'LineStyle','--','Color',rgb('dark blue'),'LineWidth',2);
f6 = plot(1:36, piH.m6, '--r', 'LineWidth', 2);
f7 = plot(1:36, piH.m7, '-b*');
f8 = plot(1:36, piH.m8, '--mo');
axis('tight')
ylim([-0.2 0.08])
box off
grid on
set(gca,'GridLineStyle','-','Layer','bottom')
set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36'},'XTick',[12 24 36 48])
legend([f6 f7 f8], {'Economic sentiment', 'Industry confidence', 'Consumer mood swings'},'Location','SouthWest','FontSize',10)
legend boxoff

subplot(4,2,6)
plot1=plot(1:36,piL.inf','LineStyle','-','Color','k','LineWidth',1);hold on
plot1=plot(1:36,piL.sup'','LineStyle','-','Color','k','LineWidth',1);
plot(1:36, piL.m1, 'k', 'LineWidth', 2)
plot(1:36, piL.m6, '--r', 'LineWidth', 2)
plot(1:36, piL.m7, '-b*')
plot(1:36, piL.m8, '--mo');
axis('tight')
ylim([-0.2 0.08])
box off
grid on
set(gca,'GridLineStyle','-','Layer','bottom')
set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36'},'XTick',[12 24 36 48])

subplot(4,2,7)
plot1 = shadedplot(1:36,piH.inf',piH.sup',[220, 207, 244] / 244,[220, 207, 244] / 244 );
hold on
plot(1:36,piH.m1,'LineStyle','--','Color',rgb('dark blue'),'LineWidth',2);
f9 = plot(1:36, piH.m9, '--r', 'LineWidth', 2);
f10 = plot(1:36, piH.m10, '-b*');
axis('tight')
ylim([-0.2 0.08])
box off
grid on
set(gca,'GridLineStyle','-','Layer','bottom')
set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36'},'XTick',[12 24 36 48])
legend([f9 f10], {'Consumer confidence','Industry confidence'},'Location','SouthWest','FontSize',10)
legend boxoff
title({'using commercial car registrations'},'FontSize',10)

subplot(4,2,8)
plot1=plot(1:36,piL.inf','LineStyle','-','Color','k','LineWidth',1);hold on
plot1=plot(1:36,piL.sup','LineStyle','-','Color','k','LineWidth',1);
plot(1:36, piL.m1, 'k', 'LineWidth', 2)
plot(1:36, piL.m9, '--r', 'LineWidth', 2);
plot(1:36, piL.m10, '-b*');
axis('tight')
ylim([-0.2 0.08])
box off
grid on
set(gca,'GridLineStyle','-','Layer','bottom')
set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36'},'XTick',[12 24 36 48])
title({'using commercial car registrations'},'FontSize',10)

set(gcf, 'position', [0 0 600 1100]);
set(gcf,'paperpositionmode','auto')
print(gcf,'-depsc2','-loose','paper_results/fig5.eps');



%% Model 11: Bond spread 
ordend = 1;
ordend2 = 5;
spec = strcat('m11_bondspreadGilch');
labels_sel =  {'VSTOXX','NFCSprd','ManuGr12','CPIGr12','ConsConfEC','EONIA','VSTOXXxConsConfEC'};
defstate.growth = 0;
defstate.cumul = 0;

r01_subset_data;

r02_ivar;


%% Model 12: Shadow rate by DeRezende
ordend = 1;
ordend2 = 4;
spec = strcat('m12_shadowDR');
labels_sel =  {'VSTOXX','ManuGr12','CPIGr12','ConsConfEC','ShadowDR','VSTOXXxConsConfEC'};

r01_subset_data;

r02_ivar;


%% Model 13: Shadow rate by Wu & Xia
ordend = 1;
ordend2 = 4;
spec = strcat('m13_shadowWX');
labels_sel =  {'VSTOXX','ManuGr12','CPIGr12','ConsConfEC','ShadowWX','VSTOXXxConsConfEC'};

r01_subset_data;

r02_ivar;


%% Model 14: Add CB balance sheet
ordend = 1;
ordend2 = 4;
spec = strcat('m14_balancesheet');
labels_sel =  {'VSTOXX','ManuGr12','CPIGr12','ConsConfEC','EONIA','ECBAssets','VSTOXXxConsConfEC'};

r01_subset_data;

r02_ivar;


%% Model 15: VSTOXX orthogonalized
ordend = 1;
ordend2 = 4;
spec = strcat('m15_vstoxxres');
labels_sel =  {'VSTOXXres','ManuGr12','CPIGr12','ConsConfEC','EONIA','VSTOXXresxConsConfEC'};

r01_subset_data;

r02_ivar;


%% Model 16: Consumer confidence orthogonalized
ordend = 1;
ordend2 = 4;
spec = strcat('M16_consconfres');
labels_sel =  {'VSTOXX','ManuGr12','CPIGr12','ConsConfECres','EONIA','VSTOXXxConsConfECres'};

r01_subset_data;

r02_ivar;


%% Model 17: Identification using gold price
dobootstrap = 0;
ordend = 1;
ordend2 = 4;
spec = strcat('m11_gold');
labels_sel =  {'VIX','ManuGr12','CPIGr12','ConsConfEC','EONIA','VIXxConsConfEC'};
proxyvar = 1;
instr_label = {'GoldProxy'};

r01_subset_data;

r02_ivar;


%% Model 18: Identification using Bloom dates for VSTOXX
ordend = 1;
ordend2 = 4;
spec = strcat('m12_Bloomdummyvstoxx');
labels_sel =  {'VSTOXX','ManuGr12','CPIGr12','ConsConfEC','EONIA','VSTOXXxConsConfEC'};
proxyvar = 1;
instr_label = {'BloomDummyVSTOXX'};

r01_subset_data;

r02_ivar;

proxyvar = 0;
dobootstrap = 1;


%% Figure 6: Load data

load('./m11_bondspreadGilch/out.mat')
ipL.m2 = OIRFavg(:,3,1);
ipH.m2 = OIRF2avg(:,3,1);
piL.m2 = OIRFavg(:,4,1);
piH.m2 = OIRF2avg(:,4,1);
load('./m12_shadowDR/out.mat')
ipL.m3 = OIRFavg(:,2,1);
ipH.m3 = OIRF2avg(:,2,1);
piL.m3 = OIRFavg(:,3,1);
piH.m3 = OIRF2avg(:,3,1);
load('./m13_shadowWX/out.mat')
ipL.m4 = OIRFavg(:,2,1);
ipH.m4 = OIRF2avg(:,2,1);
piL.m4 = OIRFavg(:,3,1);
piH.m4 = OIRF2avg(:,3,1);
load('./m14_balancesheet/out.mat')
ipL.m5 = OIRFavg(:,2,1);
ipH.m5 = OIRF2avg(:,2,1);
piL.m5 = OIRFavg(:,3,1);
piH.m5 = OIRF2avg(:,3,1);
load('./m15_vstoxxres/out.mat')
ipL.m6 = OIRFavg(:,2,1);
ipH.m6 = OIRF2avg(:,2,1);
piL.m6 = OIRFavg(:,3,1);
piH.m6 = OIRF2avg(:,3,1);
load('./M16_consconfres/out.mat')
ipL.m7 = OIRFavg(:,2,1);
ipH.m7 = OIRF2avg(:,2,1);
piL.m7 = OIRFavg(:,3,1);
piH.m7 = OIRF2avg(:,3,1);
load('./M11_Gold/out_nobtstrp.mat')
ipL.m8 = OIRFavg(:,2,1);
ipH.m8 = OIRF2avg(:,2,1);
piL.m8 = OIRFavg(:,3,1);
piH.m8 = OIRF2avg(:,3,1);
load('./m12_Bloomdummyvstoxx/out_nobtstrp.mat')
ipL.m9 = OIRFavg(:,2,1);
ipH.m9 = OIRF2avg(:,2,1);
piL.m9 = OIRFavg(:,3,1);
piH.m9 = OIRF2avg(:,3,1);


%% Figure 6: Robustness II 

figure()

subplot(4,2,1)
plot1 = shadedplot(1:36,ipH.inf',ipH.sup',[220, 207, 244] / 244,[220, 207, 244] / 244 );
hold on
plot(1:36,ipH.m1,'LineStyle','--','Color',rgb('dark blue'),'LineWidth',2);
f2 = plot(1:36, ipH.m2, '--r', 'LineWidth', 2);
f3 = plot(1:36, ipH.m3, '-b*');
f4 = plot(1:36, ipH.m4, '--mo');
f5 = plot(1:36, ipH.m5, '-cd');
axis('tight')
ylim([-1.8 0.7])
box off
grid on
set(gca,'GridLineStyle','-','Layer','bottom')
set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36'},'XTick',[12 24 36 48])
legend([f2 f3 f4 f5], {'Bond spread','Shadow rate (DR&R)','Shadow rate (W&X)','ECB assets'},'Location','SouthWest','FontSize',9)
title({'Industrial production:','Normal times'},'FontSize',14)
legend boxoff

subplot(4,2,2)
plot1=plot(1:36,ipL.inf','LineStyle','-','Color','k','LineWidth',1);hold on
plot1=plot(1:36,ipL.sup'','LineStyle','-','Color','k','LineWidth',1);
plot(1:36, ipL.m1, 'k', 'LineWidth', 2)
plot(1:36, ipL.m2, '--r', 'LineWidth', 2)
plot(1:36, ipL.m3, '-b*')
plot(1:36, ipL.m4, '--mo')
plot(1:36, ipL.m5, '-cd');
axis('tight')
ylim([-1.8 0.7])
box off
grid on
set(gca,'GridLineStyle','-','Layer','bottom')
set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36'},'XTick',[12 24 36 48])
title({'Industrial production:','Pessimistic times'},'FontSize',14)

subplot(4,2,3)
plot1 = shadedplot(1:36,ipH.inf',ipH.sup',[220, 207, 244] / 244,[220, 207, 244] / 244 );
hold on
plot(1:36,ipH.m1,'LineStyle','--','Color',rgb('dark blue'),'LineWidth',2);
f6 = plot(1:36, ipH.m6, '--r', 'LineWidth', 2);
f7 = plot(1:36, ipH.m7, '-b*');
f8 = plot(1:36, ipH.m8, '--mo');
f9 = plot(1:36, ipH.m9, '-cd');
axis('tight')
ylim([-1.8 0.7])
box off
grid on
set(gca,'GridLineStyle','-','Layer','bottom')
set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36'},'XTick',[12 24 36 48])
legend([f6 f7 f8 f9] , {'Orthogonalized uncertainty','Orthogonalized confidence', 'Gold price instrument','Bloom dummy identification'},'Location','SouthWest','FontSize',9)
legend boxoff

subplot(4,2,4)
plot1=plot(1:36,ipL.inf','LineStyle','-','Color','k','LineWidth',1);hold on
plot1=plot(1:36,ipL.sup'','LineStyle','-','Color','k','LineWidth',1);
hold on
plot(1:36, ipL.m1, 'k', 'LineWidth', 2)
plot(1:36, ipL.m6, '--r', 'LineWidth', 2)
plot(1:36, ipL.m7, '-b*')
plot(1:36, ipL.m8, '--mo');
plot(1:36, ipL.m9, '-cd');
axis('tight')
ylim([-1.8 0.7])
box off
grid on
set(gca,'GridLineStyle','-','Layer','bottom')
set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36'},'XTick',[12 24 36 48])

subplot(4,2,5)
plot1 = shadedplot(1:36,piH.inf',piH.sup',[220, 207, 244] / 244,[220, 207, 244] / 244 );
hold on
plot(1:36,piH.m1,'LineStyle','--','Color',rgb('dark blue'),'LineWidth',2);
f2 = plot(1:36, piH.m2, '--r', 'LineWidth', 2);
f3 = plot(1:36, piH.m3, '-b*');
f4 = plot(1:36, piH.m4, '--mo');
f5 = plot(1:36, piH.m5, '-cd');
axis('tight')
ylim([-0.2 0.08])
box off
grid on
set(gca,'GridLineStyle','-','Layer','bottom')
set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36'},'XTick',[12 24 36 48])
title({'Inflation:','Normal times'},'FontSize',14)
legend([f2 f3 f4 f5], {'Bond spread','Shadow rate (DR&R)','Shadow rate (W&X)','ECB assets'},'Location','SouthWest','FontSize',9)
legend boxoff

subplot(4,2,6)
plot1=plot(1:36,piL.inf','LineStyle','-','Color','k','LineWidth',1);hold on
plot1=plot(1:36,piL.sup'','LineStyle','-','Color','k','LineWidth',1);
plot(1:36, piL.m1, 'k', 'LineWidth', 2)
plot(1:36, piL.m2, '--r', 'LineWidth', 2)
plot(1:36, piL.m3, '-b*')
plot(1:36, piL.m4, '--mo')
plot(1:36, piL.m5, '-cd');
axis('tight')
ylim([-0.2 0.08])
box off
grid on
set(gca,'GridLineStyle','-','Layer','bottom')
set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36'},'XTick',[12 24 36 48])
title({'Inflation:','Pessimistic times'},'FontSize',14)

subplot(4,2,7)
plot1 = shadedplot(1:36,piH.inf',piH.sup',[220, 207, 244] / 244,[220, 207, 244] / 244 );
hold on
plot(1:36, piH.m1, 'k', 'LineWidth', 2)
f6 = plot(1:36, piH.m6, '--r', 'LineWidth', 2);
f7 = plot(1:36, piH.m7, '-b*');
f8 = plot(1:36, piH.m8, '--mo');
f9 = plot(1:36, piH.m9, '-cd');
axis('tight')
ylim([-0.2 0.08])
box off
grid on
set(gca,'GridLineStyle','-','Layer','bottom')
set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36'},'XTick',[12 24 36 48])
legend([f6 f7 f8 f9], {'Orthogonalized uncertainty','Orthogonalized confidence', 'Gold price instrument','Bloom dummy identification'},'Location','SouthWest','FontSize',9)
legend boxoff
 
subplot(4,2,8)
plot1=plot(1:36,piL.inf','LineStyle','-','Color','k','LineWidth',1);hold on
plot1=plot(1:36,piL.sup'','LineStyle','-','Color','k','LineWidth',1);
plot(1:36, piL.m1, 'k', 'LineWidth', 2)
plot(1:36, piL.m6, '--r', 'LineWidth', 2)
plot(1:36, piL.m7, '-b*')
plot(1:36, piL.m8, '--mo')
plot(1:36, piL.m9, '-cd')
axis('tight')
ylim([-0.2 0.08])
box off
grid on
set(gca,'GridLineStyle','-','Layer','bottom')
set(gca,'XTickLabelMode', 'manual','XTickLabel',{'12','24','36'},'XTick',[12 24 36 48])

set(gcf, 'position', [0 0 600 1100]);
set(gcf,'paperpositionmode','auto')
print(gcf,'-depsc2','-loose','paper_results/fig6.eps');

close all
