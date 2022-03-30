% First run Develop_priors.m to develop prior samples
% Then run Master.m to develop posterior samples
% Then run PostProcess.m to create figures. 
close all
clear all

% 100 -yr GWP according to 2018 Ozone Assessment
GWP.H1301 = 6670;
GWP.H1211 = 1750;
GWP.cfc11 = 5160;
GWP.cfc12 = 10300;
GWP.cfc113 = 6080;
GWP.cfc114 = 8580;
GWP.cfc115 = 7310;
GWP.hcfc22 = 1780;
GWP.hcfc141b = 800;
GWP.hcfc142b = 2070;


ODP.H1301 = 0.5*(15.2+19);
ODP.H1211 = 0.5*(6.9+7.7);
ODP.cfc11 = 1;
ODP.cfc12 = 0.5*(0.73+0.81);
ODP.cfc113 = 0.815;
ODP.cfc114 = 0.5;
ODP.cfc115 = 0.26;
ODP.hcfc22 = 0.029;
ODP.hcfc141b = 0.5*(0.069+0.102);
ODP.hcfc142b = 0.5*(0.023+0.057);

LTconstant.H1301 = 72;
LTconstant.H1211 = 16;
LTconstant.cfc11 = 52;
LTconstant.cfc12 = 102;
LTconstant.cfc113 = 93;
LTconstant.cfc114 = 189;
LTconstant.cfc115 = 540;
LTconstant.hcfc22 = 11.9;
LTconstant.hcfc141b = 9.4;
LTconstant.hcfc142b = 18;

MFindx.H1301 = 14;
MFindx.H1211 = 12;
MFindx.cfc11 = 2;
MFindx.cfc12 = 3;
MFindx.cfc113 = 4;
MFindx.cfc114 = 5;
MFindx.cfc115 = 6;
MFindx.hcfc22 = 9;
MFindx.hcfc141b = 10;
MFindx.hcfc142b = 11;
MFindx.years = 1;


load('CFC11.mat'); 

ShortBank.cfc11 = median(CFC11post.shortBank); 
MediumBank.cfc11 = median(CFC11post.medBank); 
LongBank.cfc11 = median(CFC11post.longBank); 

MediumRF.cfc11 = CFC11post.rf_m; 
LongRF.cfc11 = CFC11post.rf_l;

ShortBank.cfc11(65) = 0;
MediumBank.cfc11(65) = median(CFC11post.medBank(:,64).*(1 - CFC11post.rf_m'));
LongBank.cfc11(65) = median(CFC11post.longBank(:,64).*(1 - CFC11post.rf_l'));

TotalBanks.cfc11 = CFC11post.shortBank + CFC11post.medBank + CFC11post.longBank;
TotalBanks.cfc11(:,65) = CFC11post.medBank(:,64).*(1 - CFC11post.rf_m') + CFC11post.longBank(:,64).*(1 - CFC11post.rf_l');


load('CFC12.mat'); 

ShortBank.cfc12a = median(CFC12post.shortBank1); 
ShortBank.cfc12b = median(CFC12post.shortBank2); 
MediumBank.cfc12 = median(CFC12post.medBank); 
LongBank.cfc12 = median(CFC12post.longBank); 

MediumRF.cfc12 = CFC12post.rf_m; 
LongRF.cfc12 = CFC12post.rf_l;

ShortBank.cfc12a(65) = 0; 
ShortBank.cfc12b(65) = 0; 
MediumBank.cfc12(65) = median(CFC12post.medBank(:,64).*(1 - CFC12post.rf_m'));
LongBank.cfc12(65) = median(CFC12post.longBank(:,64).*(1 - CFC12post.rf_l'));

TotalBanks.cfc12 = CFC12post.shortBank1 + CFC12post.shortBank2 + CFC12post.medBank + CFC12post.longBank;
TotalBanks.cfc12(:,65) = CFC12post.medBank(:,64).*(1 - CFC12post.rf_m') + CFC12post.longBank(:,64).*(1 - CFC12post.rf_l');

load('CFC113.mat'); 

ShortBank.cfc113 = median(CFC113post.shortBank); 
LongBank.cfc113 = median(CFC113post.longBank); 

LongRF.cfc113 = CFC113post.rf_l;

ShortBank.cfc113(65) = 0; 
LongBank.cfc113(65) = median(CFC113post.longBank(:,64).*(1 - CFC113post.rf_l'));

TotalBanks.cfc113 = CFC113post.shortBank + CFC113post.longBank;
TotalBanks.cfc113(:,65) = CFC113post.longBank(:,64).*(1 - CFC113post.rf_l');

load('CFC114.mat'); 

ShortBank.cfc114 = median(CFC114post.shortBank(:,end-64:end)); 
LongBank.cfc114 = median(CFC114post.longBank(:,end-64:end)); 

LongRF.cfc114 = CFC114post.rf_l;

TotalBanks.cfc114 = CFC114post.shortBank(:,end-64:end) + CFC114post.longBank(:,end-64:end);

load('CFC115.mat'); 

ShortBank.cfc115 = median(CFC115post.shortBank(:,end-64:end)); 
LongBank.cfc115 = median(CFC115post.longBank(:,end-64:end)); 

LongRF.cfc115 = CFC115post.rf_l;

TotalBanks.cfc115 = CFC115post.shortBank(:,end-64:end) + CFC115post.longBank(:,end-64:end);


load('HCFC22.mat'); 

ShortBank.hcfc22 = median(HCFC22post.shortBank); 
MediumBank.hcfc22 = median(HCFC22post.medBank); 
LongBank.hcfc22 = median(HCFC22post.longBank); 

MediumRF.hcfc22 = HCFC22post.rf_m; 
LongRF.hcfc22 = HCFC22post.rf_l;

TotalBanks.hcfc22 = HCFC22post.shortBank + HCFC22post.medBank + HCFC22post.longBank;

load('HCFC141b.mat') 

ShortBank.hcfc141b = median(HCFC141bpost.shortBank); 
MediumBank.hcfc141b = median(HCFC141bpost.medBank); 
LongBank.hcfc141b = median(HCFC141bpost.longBank); 

MediumRF.hcfc141b = HCFC141bpost.rf_m; 
LongRF.hcfc141b = HCFC141bpost.rf_l;

TotalBanks.hcfc141b = HCFC141bpost.shortBank + HCFC141bpost.medBank + HCFC141bpost.longBank;


load('HCFC142b.mat') 
ShortBank.hcfc142b = median(HCFC142bpost.shortBank); 
MediumBank.hcfc142b = median(HCFC142bpost.medBank); 
LongBank.hcfc142b = median(HCFC142bpost.longBank); 

MediumRF.hcfc142b = HCFC142bpost.rf_m; 
LongRF.hcfc142b = HCFC142bpost.rf_l;

TotalBanks.hcfc142b = HCFC142bpost.shortBank + HCFC142bpost.medBank + HCFC142bpost.longBank;

load('Halon1211.mat'); 

HalonRF.H1211 = median(halon1211post.rf);
HalonEmiss.H1211 = median(halon1211post.Emiss);
HalonBank.H1211 = median(halon1211post.Bank);

TotalBanks.H1211 = halon1211post.Bank;

load('Halon1301.mat'); 

HalonRF.H1301 = median(halon1301post.rf);
HalonEmiss.H1301 = median(halon1301post.Emiss);
HalonBank.H1301 = median(halon1301post.Bank);

TotalBanks.H1301 = halon1301post.Bank;

fighandle = figure; 
subplot(3,4,1)
Bank11(1,:) = ShortBank.cfc11;
Bank11(2,:) = MediumBank.cfc11;
Bank11(3,:) = LongBank.cfc11;
area([1955:2019],0.001*Bank11'); title('CFC-11 Bank'); ylabel('[Gg]')
legend('aerosols/oc foam','refrigeration','cc foam', 'Location','northwest');

subplot(3,4,2)
Bank12(1,:) = ShortBank.cfc12a + ShortBank.cfc12b;
Bank12(2,:) = MediumBank.cfc12;
Bank12(3,:) = LongBank.cfc12;
area([1955:2019],0.001*Bank12'); title('CFC-12 Bank'); ylabel('[Gg]')
legend('aerosols/oc foam','non-h ref.','ref.', 'Location','northwest');

subplot(3,4,3)
Bank113(1,:) = ShortBank.cfc113;
Bank113(2,:) = LongBank.cfc113;
area([1955:2019],0.001*Bank113'); title('CFC-113 Bank'); ylabel('[Gg]')
legend('solvents','heat pumps', 'Location','northwest');

subplot(3,4,4)
Bank114(1,:) = ShortBank.cfc114;
Bank114(2,:) = LongBank.cfc114;
area([1955:2019],0.001*Bank114'); title('CFC-114 Bank'); ylabel('[Gg]')
legend('aerosols','heat pumps', 'Location','northwest');

subplot(3,4,5)
Bank115(1,:) = ShortBank.cfc115;
Bank115(2,:) = LongBank.cfc115;
area([1955:2019],0.001*Bank115'); title('CFC-115 Bank'); ylabel('[Gg]')
legend('propellant','A/C', 'Location','northwest');


subplot(3,4,6)
Bank22(1,:) = ShortBank.hcfc22;
Bank22(2,:) = MediumBank.hcfc22;
Bank22(3,:) = LongBank.hcfc22;
area([1944:2019],Bank22'); title('HCFC-22 Bank'); ylabel('[Gg]')
legend('oc foam','non-h ref.','foam', 'Location','northwest');

subplot(3,4,7)
Bank141(1,:) = ShortBank.hcfc141b;
Bank141(2,:) = MediumBank.hcfc141b;
Bank141(3,:) = LongBank.hcfc141b;
area([1989:2019],Bank141'); title('HCFC-141b Bank'); ylabel('[Gg]')
legend('oc foam','non-h ref.','foam', 'Location','northwest');

subplot(3,4,8)
Bank142(1,:) = ShortBank.hcfc142b;
Bank142(2,:) = MediumBank.hcfc142b;
Bank142(3,:) = LongBank.hcfc142b;
area([1981:2019],Bank142'); title('HCFC-142b Bank'); ylabel('[Gg]')
legend('oc foam','non-h ref.','foam', 'Location','northwest');

subplot(3,4,9)
area([2019-57:2019],0.001*HalonBank.H1211); title('Halon 1211'); ylabel('[Gg]')
legend('fire extinguisher', 'Location','northwest');

subplot(3,4,10)
area([2019-57:2019],0.001*HalonBank.H1301); title('Halon 1301'); ylabel('[Gg]')
legend('fire extinguisher', 'Location','northwest');


figure_width = 18; % in inches
figure_height = 9; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(fighandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

FigName = 'Fig06';
Figstr = strcat('Figures/',FigName,'.pdf'); 
print(gcf, '-dpdf', Figstr);

fighandle = figure(2)
subplot(1,3,1); 
BankTotals(1:65,1) = 0.001*sum(Bank11,1)'; 
BankTotals(1:65,2) = 0.001*sum(Bank12,1)'; 
BankTotals(1:65,3) = 0.001*sum(Bank113,1)'; 
BankTotals(1:65,4) = 0.001*sum(Bank114,1)'; 
BankTotals(1:65,5) = 0.001*sum(Bank115,1)'; 
BankTotals(1:65,6) = sum(Bank22(:,12:end),1)'; 
BankTotals(35:65,7) = sum(Bank141,1)'; 
BankTotals(27:65,8) = sum(Bank142,1)'; 
BankTotals(65-57:65,9) = 0.001*HalonBank.H1211';
BankTotals(65-57:65,10) = 0.001*HalonBank.H1301';
h = area([1955:2019],BankTotals); ylabel('[Gg]'); 

HCFC_bankfrac = sum(BankTotals(65, 6:8))/sum(BankTotals(65, :))

h(1).FaceColor =      [    0    0.4470    0.7410];
h(2).FaceColor =      [0.8500    0.3250    0.0980];
h(3).FaceColor =      [0.3010    0.7450    0.9330];
h(4).FaceColor =      [0.4940    0.1840    0.5560];
h(5).FaceColor =      [0.5    0    1];
h(6).FaceColor =      [0.9290    0.6940    0.1250];
h(7).FaceColor =      [0.6350    0.0780    0.1840];
h(8).FaceColor =      [0    0    1];
h(9).FaceColor =      [1    0    1];
h(10).FaceColor =     [0.4660    0.6740    0.1880];
title('Banks by mass'); xlim([1955,2019])
legend('CFC-11', 'CFC-12','CFC-113', 'CFC-114','CFC-115','HCFC22','HCFC141b','HCFC142b','Halon1211','Halon1301','Location','northwest');

subplot(1,3,2); 
BankTotals(1:65,1) = GWP.cfc11*0.001*sum(Bank11,1)'; 
BankTotals(1:65,2) = GWP.cfc12*0.001*sum(Bank12,1)'; 
BankTotals(1:65,3) = GWP.cfc113*0.001*sum(Bank113,1)'; 
BankTotals(1:65,4) = GWP.cfc114*0.001*sum(Bank114,1)'; 
BankTotals(1:65,5) = GWP.cfc115*0.001*sum(Bank115,1)'; 
BankTotals(1:65,6) = GWP.hcfc22*sum(Bank22(:,12:end),1)'; 
BankTotals(35:65,7) = GWP.hcfc141b*sum(Bank141,1)'; 
BankTotals(27:65,8) = GWP.hcfc142b*sum(Bank142,1)'; 
BankTotals(65-57:65,9) = GWP.H1211*0.001*HalonBank.H1211';
BankTotals(65-57:65,10) = GWP.H1301*0.001*HalonBank.H1301';
h = area([1955:2019],BankTotals); ylabel('[Gg] of CO2eq'); 

% print the fraction of banks by GWP attributed to CFC11, 12 and HCFC22
CFC11_bankfrac = sum(BankTotals(65, 1))/sum(BankTotals(65, :))
CFC12_bankfrac = sum(BankTotals(65, 2))/sum(BankTotals(65, :))
HCFC22_bankfrac = sum(BankTotals(65, 6))/sum(BankTotals(65, :))

h(1).FaceColor =      [    0    0.4470    0.7410];
h(2).FaceColor =      [0.8500    0.3250    0.0980];
h(3).FaceColor =      [0.3010    0.7450    0.9330];
h(4).FaceColor =      [0.4940    0.1840    0.5560];
h(5).FaceColor =      [0.5    0    1];
h(6).FaceColor =      [0.9290    0.6940    0.1250];
h(7).FaceColor =      [0.6350    0.0780    0.1840];
h(8).FaceColor =      [0    0    1];
h(9).FaceColor =      [1    0    1];
h(10).FaceColor =     [0.4660    0.6740    0.1880];
title('Banks by CO2 eq GWP'); xlim([1955,2019])

subplot(1,3,3); 
BankTotals(1:65,1) = ODP.cfc11*0.001*sum(Bank11,1)'; 
BankTotals(1:65,2) = ODP.cfc12*0.001*sum(Bank12,1)'; 
BankTotals(1:65,3) = ODP.cfc113*0.001*sum(Bank113,1)'; 
BankTotals(1:65,4) = ODP.cfc114*0.001*sum(Bank114,1)'; 
BankTotals(1:65,5) = ODP.cfc115*0.001*sum(Bank115,1)'; 
BankTotals(1:65,6) = ODP.hcfc22*sum(Bank22(:,12:end),1)'; 
BankTotals(35:65,7) = ODP.hcfc141b*sum(Bank141,1)'; 
BankTotals(27:65,8) = ODP.hcfc142b*sum(Bank142,1)'; 
BankTotals(65-57:65,9) = ODP.H1211*0.001*HalonBank.H1211';
BankTotals(65-57:65,10) = ODP.H1301*0.001*HalonBank.H1301';
h = area([1955:2019],BankTotals); ylabel('[Gg] of CFC-11 ODP eq'); 

CFC1112_bankfrac = sum(BankTotals(65, 1:2))/sum(BankTotals(65, :))
Halon_bankfrac = sum(BankTotals(65, 9:10))/sum(BankTotals(65, :))

h(1).FaceColor =      [    0    0.4470    0.7410];
h(2).FaceColor =      [0.8500    0.3250    0.0980];
h(3).FaceColor =      [0.3010    0.7450    0.9330];
h(4).FaceColor =      [0.4940    0.1840    0.5560];
h(5).FaceColor =      [0.5    0    1];
h(6).FaceColor =      [0.9290    0.6940    0.1250];
h(7).FaceColor =      [0.6350    0.0780    0.1840];
h(8).FaceColor =      [0    0    1];
h(9).FaceColor =      [1    0    1];
h(10).FaceColor =     [0.4660    0.6740    0.1880];
title('Banks by CFC-11 ODP eq'); xlim([1955,2019])

figure_width = 13; % in inches
figure_height = 4.5; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(fighandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

FigName = 'Fig05';
Figstr = strcat('Figures/',FigName,'.pdf'); 
print(gcf, '-dpdf', Figstr);


BanksEnd(:,1) = TotalBanks.cfc11(:,end); 
BanksEnd(:,2) = TotalBanks.cfc12(:,end); 
BanksEnd(:,3) = TotalBanks.hcfc22(:,end); 
BanksEnd(:,4) = TotalBanks.hcfc141b(:,end); 
BanksEnd(:,5) = TotalBanks.hcfc142b(:,end); 
BanksEnd(:,6) = TotalBanks.H1211(:,end); 
BanksEnd(:,7) = TotalBanks.H1301(:,end); 

ReportVals(1,:) = 0.001*prctile(BanksEnd, 50); 
ReportVals(2,:) = 0.001*prctile(BanksEnd, 2.5); 
ReportVals(3,:) = 0.001*prctile(BanksEnd, 97.5); 

%% Mixing Ratios, Fig 1

load('AGAGE_published.mat'); 
clear MRobs_agage
tmp = MF;
Yrobs_agage = MF(:,1);
tmp_yr = floor(Yrobs_agage); 
Yrobs_agage = unique(tmp_yr);
for yy = 1:length(Yrobs_agage)
    ind = find(tmp_yr == Yrobs_agage(yy)); 
    MRobs_agage(yy,:) = mean(tmp(ind,2:11),1);
end

load('MixingRatios.mat'); 

fighandle = figure(3) 
subplot(3,4,1)
MRind = 2; 
MRobs = MixingRatios(:,MRind);
YRobs = MixingRatios(:,1);
ytmp3 = find(YRobs == 1980); 
ytmp4 = find(YRobs == 2018); 

MED = prctile(CFC11post.MF, 50);
LB = prctile(CFC11post.MF, 5);
UB = prctile(CFC11post.MF, 95);

p1 = boundedline(CFC11post.Years(2:end), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
hold on; 
p2 = plot(1980:2018 , MRobs(ytmp3:ytmp4), 'b', 'LineWidth',2); 
box on; title('CFC-11'); ylabel('Mole Fractions [ppt]');
xlim([1960,2020]); 
legend([p1, p2], 'posteriors', 'observed', 'Location','northwest'); 

subplot(3,4,2)
MRind = 3; 
MRobs = MixingRatios(:,MRind);
YRobs = MixingRatios(:,1);
ytmp3 = find(YRobs == 1980); 
ytmp4 = find(YRobs == 2018); 

MED = prctile(CFC12post.MF, 50);
LB = prctile(CFC12post.MF, 5);
UB = prctile(CFC12post.MF, 95);

p1 = boundedline(CFC12post.Years(2:end), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
hold on; 
p2 = plot(1980:2018 , MRobs(ytmp3:ytmp4), 'b', 'LineWidth',2); 
box on; title('CFC-12'); ylabel('Mole Fractions [ppt]');
xlim([1960,2020]); 

subplot(3,4,3)
MRind = 4; 
MRobs = MixingRatios(:,MRind);
YRobs = MixingRatios(:,1);
ytmp3 = find(YRobs == 1980); 
ytmp4 = find(YRobs == 2018); 

MED = prctile(CFC113post.MF, 50);
LB = prctile(CFC113post.MF, 5);
UB = prctile(CFC113post.MF, 95);

p1 = boundedline(CFC113post.Years+1, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
hold on; 
p2 = plot(Yrobs_agage , MRobs_agage(:,3), 'b', 'LineWidth',2); 
box on; title('CFC-113'); ylabel('Mole Fractions [ppt]');
xlim([1960,2020]); 

subplot(3,4,4)
MRind = 5; 
MRobs = MixingRatios(:,MRind);
YRobs = MixingRatios(:,1);
ytmp3 = find(YRobs == 2004); 
ytmp4 = find(YRobs == 2019); 

MED = prctile(CFC114post.MF, 50);
LB = prctile(CFC114post.MF, 5);
UB = prctile(CFC114post.MF, 95);

p1 = boundedline(CFC114post.Years(2:end-9), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
hold on; 
p2 = plot(Yrobs_agage , MRobs_agage(:,4), 'b', 'LineWidth',2); 
box on; title('CFC-114'); ylabel('Mole Fractions [ppt]');
xlim([1960,2020]); 

subplot(3,4,5)
MRind = 6; 
MRobs = MixingRatios(:,MRind);
YRobs = MixingRatios(:,1);
ytmp3 = find(YRobs == 2004); 
ytmp4 = find(YRobs == 2019); 

MED = prctile(CFC115post.MF, 50);
LB = prctile(CFC115post.MF, 5);
UB = prctile(CFC115post.MF, 95);

p1 = boundedline(CFC115post.Years(2:end-9), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
hold on; 
p2 = plot(Yrobs_agage , MRobs_agage(:,5), 'b', 'LineWidth',2); 
box on; title('CFC-115'); ylabel('Mole Fractions [ppt]');
xlim([1960,2020]); 

subplot(3,4,6)
MRind = 9; 
MRobs = MixingRatios(:,MRind);
YRobs = MixingRatios(:,1);
ytmp3 = find(YRobs == 1998); 
ytmp4 = find(YRobs == 2019); 

MED = prctile(HCFC22post.MF, 50);
LB = prctile(HCFC22post.MF, 5);
UB = prctile(HCFC22post.MF, 95);

p1 = boundedline(HCFC22post.Years, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
hold on; 
p2 = plot(Yrobs_agage , MRobs_agage(:,6), 'b', 'LineWidth',2); 
box on; title('HCFC-22'); ylabel('Mole Fractions [ppt]');
xlim([1960,2020]); 

subplot(3,4,7)
MRind = 10; 
MRobs = MixingRatios(:,MRind);
YRobs = MixingRatios(:,1);
ytmp3 = find(YRobs == 1998); 
ytmp4 = find(YRobs == 2019); 

MED = prctile(HCFC141bpost.MF, 50);
LB = prctile(HCFC141bpost.MF, 5);
UB = prctile(HCFC141bpost.MF, 95);

p1 = boundedline(HCFC141bpost.Years(2:32), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
hold on; 
p2 = plot(Yrobs_agage , MRobs_agage(:,7), 'b', 'LineWidth',2); 
box on; title('HCFC-141b'); ylabel('Mole Fractions [ppt]');
xlim([1960,2020]); 

subplot(3,4,8)
MRind = 11; 
MRobs = MixingRatios(:,MRind);
YRobs = MixingRatios(:,1);
ytmp3 = find(YRobs == 1998); 
ytmp4 = find(YRobs == 2019); 

MED = prctile(HCFC142bpost.MF, 50);
LB = prctile(HCFC142bpost.MF, 5);
UB = prctile(HCFC142bpost.MF, 95);

p1 = boundedline(HCFC142bpost.Years(1:39), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
hold on; 
p2 = plot(Yrobs_agage , MRobs_agage(:,8), 'b', 'LineWidth',2); 
box on; title('HCFC-142b'); ylabel('Mole Fractions [ppt]');
xlim([1960,2020]); 

subplot(3,4,9)
MRind = 12; 
MRobs = MixingRatios(:,MRind);
YRobs = MixingRatios(:,1);
ytmp3 = find(YRobs == 1994); 
ytmp4 = find(YRobs == 2019); 

MED = prctile(halon1211post.MF, 50);
LB = prctile(halon1211post.MF, 5);
UB = prctile(halon1211post.MF, 95);

p1 = boundedline(halon1211post.Years, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
hold on; 
p2 = plot(Yrobs_agage , MRobs_agage(:,9), 'b', 'LineWidth',2); 
box on; title('Halon-1211'); ylabel('Mole Fractions [ppt]');
xlim([1960,2020]); 

subplot(3,4,10)
MRind = 14; 
MRobs = MixingRatios(:,MRind);
YRobs = MixingRatios(:,1);
ytmp3 = find(YRobs == 1994); 
ytmp4 = find(YRobs == 2019); 

MED = prctile(halon1301post.MF, 50);
LB = prctile(halon1301post.MF, 5);
UB = prctile(halon1301post.MF, 95);

p1 = boundedline(halon1301post.Years, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
hold on; 
p2 = plot(Yrobs_agage , MRobs_agage(:,10), 'b', 'LineWidth',2); 
box on; title('Halon-1301'); ylabel('Mole Fractions [ppt]');
xlim([1960,2020]); 

figure_width = 18; % in inches
figure_height = 9; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(fighandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

FigName = 'Fig01';
Figstr = strcat('Figures/',FigName,'.pdf'); 
print(gcf, '-dpdf', Figstr);

%% Figure 2: Emissions
 

fighandle = figure(4) 
subplot(3,4,1)

Emiss = CFC11post.emiss_s + CFC11post.emiss_m + CFC11post.emiss_l;
MED = 0.001*prctile(Emiss, 50);
LB = 0.001*prctile(Emiss, 5);
UB = 0.001*prctile(Emiss, 95);

p1 = boundedline(CFC11post.Years(2:end), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
hold on; 
p2 = plot(1980:2018 , 0.001*CFC11post.ObsEmiss, 'b', 'LineWidth',2); 
box on; title('CFC-11'); ylabel('Emissions [Gg yr^{-1}]');
xlim([1960,2020]); 
legend([p1, p2], 'posteriors', 'obs-derived', 'Location','northeast'); 

subplot(3,4,2)
Emiss = CFC12post.emiss_s1 + CFC12post.emiss_s2 + CFC12post.emiss_m + CFC12post.emiss_l;
MED = 0.001*prctile(Emiss, 50);
LB = 0.001*prctile(Emiss, 5);
UB = 0.001*prctile(Emiss, 95);

p1 = boundedline(CFC12post.Years(2:end), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
hold on; 
p2 = plot(1980:2018 , 0.001*CFC12post.ObsEmiss, 'b', 'LineWidth',2); 
box on; title('CFC-12'); ylabel('Emissions [Gg yr^{-1}]');
xlim([1960,2020]); 

subplot(3,4,3)
Emiss = CFC113post.emiss_s + CFC113post.emiss_l + CFC113post.emiss_f;
MED = 0.001*prctile(Emiss, 50);
LB = 0.001*prctile(Emiss, 5);
UB = 0.001*prctile(Emiss, 95);

p1 = boundedline(CFC113post.Years+1, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
hold on; 
p2 = plot(2004:2018 , 0.001*CFC113post.ObsEmiss, 'b', 'LineWidth',2); 
box on; title('CFC-113'); ylabel('Emissions [Gg yr^{-1}]');
xlim([1960,2020]); 

subplot(3,4,4)
Emiss = CFC114post.emiss_s +  CFC114post.emiss_l;
MED = prctile(Emiss, 50);
LB = prctile(Emiss, 5);
UB = prctile(Emiss, 95);

p1 = boundedline(CFC114post.Years(2:end-9), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
hold on; 
p2 = plot(2004:2019, 0.001*CFC114post.ObsEmiss, 'b', 'LineWidth',2); 
box on; title('CFC-114'); ylabel('Emissions [Gg yr^{-1}]');
xlim([1960,2020]); 


subplot(3,4,5)
Emiss = CFC115post.emiss_s +  CFC115post.emiss_l;
MED = prctile(Emiss, 50);
LB = prctile(Emiss, 5);
UB = prctile(Emiss, 95);

p1 = boundedline(CFC115post.Years(2:end-9), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
hold on; 
p2 = plot(2004:2019, 0.001*CFC115post.ObsEmiss, 'b', 'LineWidth',2); 
box on; title('CFC-115'); ylabel('Emissions [Gg yr^{-1}]');
xlim([1960,2020]); 


subplot(3,4,6)
Emiss = HCFC22post.emiss_s + HCFC22post.emiss_m + HCFC22post.emiss_l + HCFC22post.emiss_f;
MED = prctile(Emiss, 50);
LB = prctile(Emiss, 5);
UB = prctile(Emiss, 95);

p1 = boundedline(HCFC22post.Years, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
hold on; 
p2 = plot(2004:2019, 0.001*HCFC22post.ObsEmiss, 'b', 'LineWidth',2); 
box on; title('HCFC-22'); ylabel('Emissions [Gg yr^{-1}]');
xlim([1960,2020]); 


subplot(3,4,7)
Emiss = HCFC141bpost.emiss_s + HCFC141bpost.emiss_m + HCFC141bpost.emiss_l;
MED = prctile(Emiss, 50);
LB = prctile(Emiss, 5);
UB = prctile(Emiss, 95);

p1 = boundedline(HCFC141bpost.Years(2:32), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
hold on; 
p2 = plot(2004:2019, 0.001*HCFC141bpost.ObsEmiss, 'b', 'LineWidth',2); 
box on; title('HCFC-141b'); ylabel('Emissions [Gg yr^{-1}]');
xlim([1960,2020]); 


subplot(3,4,8)
Emiss = HCFC142bpost.emiss_s + HCFC142bpost.emiss_m + HCFC142bpost.emiss_l + repmat(HCFC142bpost.de_f',1,39).*HCFC142bpost.prod_f(:,1:39);
MED = prctile(Emiss, 50);
LB = prctile(Emiss, 5);
UB = prctile(Emiss, 95);

p1 = boundedline(HCFC142bpost.Years(1:39), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
hold on; 
p2 = plot(2004:2019, 0.001*HCFC142bpost.ObsEmiss, 'b', 'LineWidth',2); 
box on; title('HCFC-142b'); ylabel('Emissions [Gg yr^{-1}]');
xlim([1960,2020]); 

subplot(3,4,9)

MED = 0.001*prctile(halon1211post.Emiss, 50);
LB = 0.001*prctile(halon1211post.Emiss, 5);
UB = 0.001*prctile(halon1211post.Emiss, 95);

p1 = boundedline(halon1211post.Years, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
hold on; 
p2 = plot(2004:2019 , 0.001*halon1211post.ObsEmiss, 'b', 'LineWidth',2); 
box on; title('Halon-1211'); ylabel('Emissions [Gg yr^{-1}]');
xlim([1960,2020]); 

subplot(3,4,10)
MED = 0.001*prctile(halon1301post.Emiss, 50);
LB = 0.001*prctile(halon1301post.Emiss, 5);
UB = 0.001*prctile(halon1301post.Emiss, 95);


p1 = boundedline(halon1301post.Years, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
hold on; 
p2 = plot(2004:2019 , 0.001*halon1301post.ObsEmiss, 'b', 'LineWidth',2); 
box on; title('Halon-1301'); ylabel('Emissions [Gg yr^{-1}]');
xlim([1960,2020]); 

figure_width = 18; % in inches
figure_height = 9; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(fighandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

FigName = 'Fig02';
Figstr = strcat('Figures/',FigName,'.pdf'); 
print(gcf, '-dpdf', Figstr);


%% Figure 3: Banks
 

fighandle = figure(5) 
subplot(3,4,1)

Banks = CFC11post.shortBank + CFC11post.medBank  + CFC11post.longBank ;
MED = 0.001*prctile(Banks, 50);
LB = 0.001*prctile(Banks, 5);
UB = 0.001*prctile(Banks, 95);

p1 = boundedline(CFC11post.Years(2:end), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
load('Draft1_JohnsBanks.mat', 'unnamed'); 
hold on; 
yr = unnamed(:,1); 
Johns_bank = (10^(-3))*unnamed(:,2); 
ind = yr>1995;
p2 = plot(yr(ind),Johns_bank(ind) ,'b', 'LineWidth',2); 
hold on;
p3 = plot([2008,2015],[1420,982],':o','LineWidth',2,'Color',[0 0 0],'MarkerSize',10,'MarkerEdgeColor',[0 0 0]);
%txt = 'WMO(2018)';
%text(2009,900,txt,'Color',[0 0 0],'FontSize',12)

box on; title('CFC-11'); ylabel('Banks [Gg]');
xlim([1960,2020]); 
legend([p1, p2, p3], 'posteriors', 'hybrid', 'WMO(2018)', 'Location','northwest'); 

subplot(3,4,2)
Banks = CFC12post.shortBank1 + CFC12post.shortBank2 + CFC12post.medBank + CFC12post.longBank;
MED = 0.001*prctile(Banks, 50);
LB = 0.001*prctile(Banks, 5);
UB = 0.001*prctile(Banks, 95);

p1 = boundedline(CFC12post.Years(2:end), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);

load('Draft1_JohnsBanks.mat', 'unnamed'); 
hold on; 
yr = unnamed(:,1); 
Johns_bank = (10^(-3))*unnamed(:,3); 
ind = yr>1995;
p4 = plot(yr(ind),Johns_bank(ind) ,'b', 'LineWidth',2); 

% plot(2008,1420,'+k','MarkerSize',8,'MarkerEdgeColor',[1 0.3 0.3]); %Table 5A-2 from WMO 2010 report
% txt = 'TEAP(2009)';
% text(2000.5,1358,txt,'Color',[1 0.3 0.3],'FontSize',12)

% Add WMO 2018 estimate
hold on;
p5 = plot([2008,2015],[394,47],':o','LineWidth',2,'Color',[0 0 0],'MarkerSize',10,'MarkerEdgeColor',[0 0 0]);
%txt = 'WMO(2018)';
%text(2009,900,txt,'Color',[0 0 0],'FontSize',12)
box on; title('CFC-12'); ylabel('Banks [Gg ]');
legend([p1, p4, p5], 'posteriors', 'hybrid', 'WMO(2018)', 'Location','northwest'); 
xlim([1960,2020]); 


subplot(3,4,3)
Bank = CFC113post.shortBank + CFC113post.longBank;
MED = 0.001*prctile(Bank, 50);
LB = 0.001*prctile(Bank, 5);
UB = 0.001*prctile(Bank, 95);

p1 = boundedline(CFC113post.Years+1, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
load('Draft1_JohnsBanks.mat', 'unnamed'); 
% hold on; 
% yr = unnamed(:,1); 
% Johns_bank = (10^(-3))*unnamed(:,4); 
% ind = yr>1995;
% p4 = plot(yr(ind),Johns_bank(ind) ,'b', 'LineWidth',2); 

box on; title('CFC-113'); ylabel('Banks [Gg]');
xlim([1960,2020]); 

subplot(3,4,4)
Bank = CFC114post.shortBank +  CFC114post.longBank;
MED = 0.001*prctile(Bank, 50);
LB = 0.001*prctile(Bank, 5);
UB = 0.001*prctile(Bank, 95);

p1 = boundedline(CFC114post.Years(2:end-9), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
box on; title('CFC-114'); ylabel('Bank [Gg]');
xlim([1960,2020]); 


subplot(3,4,5)
Bank = CFC115post.shortBank +  CFC115post.longBank;
MED = 0.001*prctile(Bank, 50);
LB = 0.001*prctile(Bank, 5);
UB = 0.001*prctile(Bank, 95);

p1 = boundedline(CFC115post.Years(2:end-9), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
box on; title('CFC-115'); ylabel('Bank [Gg]');
xlim([1960,2020]); 


subplot(3,4,6)
Bank = HCFC22post.shortBank + HCFC22post.medBank + HCFC22post.longBank;
MED = prctile(Bank, 50);
LB = prctile(Bank, 5);
UB = prctile(Bank, 95);

p1 = boundedline(HCFC22post.Years, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);

load('Draft1_JohnsBanks.mat', 'unnamed'); 
hold on; 
yr = unnamed(:,1); 
Johns_bank = (10^(-3))*unnamed(:,4); 
ind = yr>1995;
p4 = plot(yr(ind),Johns_bank(ind) ,'b', 'LineWidth',2); 

% ADD TEAP esttimate
hold on; 
p5 = plot([1993,2008],[302,1618],':o','LineWidth',2,'Color',[0 0 0],'MarkerSize',4,'MarkerEdgeColor',[0 0 0]);

% ADD TEAP 2009 estimate
hold on; 
p6 = plot(2008,2083,'+r','MarkerSize',10,'LineWidth',2); %Table 5A-2 from WMO 2010 report

legend([p1, p4, p5, p6], 'posterior', 'hybrid','TEAP (2009)','WMO (2007)', 'Location','northwest')

box on; title('HCFC-22'); ylabel('Bank [Gg]');
xlim([1960,2020]); 


subplot(3,4,7)
Bank = HCFC141bpost.shortBank + HCFC141bpost.medBank + HCFC141bpost.longBank;
MED = prctile(Bank, 50);
LB = prctile(Bank, 5);
UB = prctile(Bank, 95);

p1 = boundedline(HCFC141bpost.Years(2:32), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
load('Draft1_JohnsBanks.mat', 'unnamed'); 
hold on; 
yr = unnamed(:,1); 
Johns_bank = 0.001*unnamed(:,5); 
ind = yr>1995;
p4 = plot(yr(ind),Johns_bank(ind) ,'b', 'LineWidth',2); 

% ADD TEAP (2009) esttimate
hold on; 
p5 = plot([1994,2008],[115,941],':o','LineWidth',2,'Color',[0 0 0],'MarkerSize',4,'MarkerEdgeColor',[0 0 0]);

% ADD WMO (2007) estimate
hold on; 
p6 = plot(2008,961,'+r','MarkerSize',10,'LineWidth',2); %Table 5A-2 from WMO 2010 report

legend([p1, p4, p5, p6], 'posterior', 'hybrid', 'TEAP(2009)','WMO(2007)', 'Location','northwest')
box on; title('HCFC-141b'); ylabel('Emissions [Gg yr^{-1}]');
xlim([1960,2020]); 


subplot(3,4,8)
Bank = HCFC142bpost.shortBank + HCFC142bpost.medBank + HCFC142bpost.longBank;
MED = prctile(Bank, 50);
LB = prctile(Bank, 5);
UB = prctile(Bank, 95);

p1 = boundedline(HCFC142bpost.Years(1:39), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);

load('Draft1_JohnsBanks.mat', 'unnamed'); 
hold on; 
yr = unnamed(:,1); 
Johns_bank = 0.001*unnamed(:,6); 
ind = yr>1995;
p4 = plot(yr(ind),Johns_bank(ind) ,'b', 'LineWidth',2); 

% ADD TEAP 2009esttimate
hold on; 
p5 = plot([1993,2008],[95,273],':o','LineWidth',2,'Color',[0 0 0],'MarkerSize',4,'MarkerEdgeColor',[0 0 0]);

% ADD WMO(2007) estimate
hold on; 
p6 = plot(2008,211,'+r','MarkerSize',10,'LineWidth',2); %Table 5A-2 from WMO 2010 report

legend([p1, p4, p5, p6], 'posterior', 'hybrid', 'TEAP (2009)','WMO (2007)', 'Location','northwest')
box on; title('HCFC-142b'); ylabel('Bank [Gg]');
xlim([1960,2020]); 

subplot(3,4,9)

MED = 0.001*prctile(halon1211post.Bank, 50);
LB = 0.001*prctile(halon1211post.Bank, 5);
UB = 0.001*prctile(halon1211post.Bank, 95);

p1 = boundedline(halon1211post.Years, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
load('TEAP_HTOC2018Banks.mat')
hold on; 
p3 = plot(TEAP_HTOC2018banks.years, 0.001*TEAP_HTOC2018banks.halon1211,'k', 'LineWidth',2); 

load('Draft1_JohnsBanks.mat', 'unnamed'); 
hold on; 
yr = unnamed(:,1); 
Johns_bank = 0.001*unnamed(:,7); 
ind = yr>1995;
p4 = plot(yr(ind),Johns_bank(ind) ,'b', 'LineWidth',2); 

% ADD TEAP 2009esttimate
hold on; 
p5 = plot([1993,2008],[128,74],':o','LineWidth',2,'Color',[0 0 0],'MarkerSize',4,'MarkerEdgeColor',[0 0 0]);

% ADD WMO(2007) estimate
hold on; 
p6 = plot(2008,86,'+r','MarkerSize',10,'LineWidth',2); %Table 5A-2 from WMO 2010 report

legend([p1, p3, p4, p5, p6], 'posterior', 'HTOC 2018','hybrid','TEAP (2009)','WMO (2007)', 'Location','northwest'); 
box on; title('Halon-1211'); ylabel('Bank [Gg]');
xlim([1960,2020]); 

subplot(3,4,10)
MED = 0.001*prctile(halon1301post.Bank, 50);
LB = 0.001*prctile(halon1301post.Bank, 5);
UB = 0.001*prctile(halon1301post.Bank, 95);


p1 = boundedline(halon1301post.Years, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);
load('TEAP_HTOC2018Banks.mat')
hold on; 
p3 = plot(TEAP_HTOC2018banks.years, 0.001*TEAP_HTOC2018banks.halon1301,'k', 'LineWidth',2); 

load('Draft1_JohnsBanks.mat', 'unnamed'); 
hold on; 
yr = unnamed(:,1); 
Johns_bank = 0.001*unnamed(:,8); 
ind = yr>1995;
p4 = plot(yr(ind),Johns_bank(ind) ,'b', 'LineWidth',2); 

% ADD TEAP 2009esttimate
hold on; 
p5 = plot([1992,2008],[74,47],':o','LineWidth',2,'Color',[0 0 0],'MarkerSize',4,'MarkerEdgeColor',[0 0 0]);

% ADD WMO(2007) estimate
hold on; 
p6 = plot(2008,31,'+r','MarkerSize',10,'LineWidth',2); %Table 5A-2 from WMO 2010 report

legend([p1, p3, p4, p5, p6], 'posterior', 'HTOC 2018','hybrid','TEAP (2009)','WMO (2007)', 'Location','northwest'); 
box on; title('Halon-1301'); ylabel('Bank [Gg]');
xlim([1960,2020]); 

figure_width = 18; % in inches
figure_height = 9; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(fighandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

FigName = 'Fig03';
Figstr = strcat('Figures/',FigName,'.pdf'); 
print(gcf, '-dpdf', Figstr);

%% Figure 4: Emissions by source 

fighandle = figure(6) 
subplot(2,4,1)

Emiss = CFC11post.emiss_s;
MED = 0.001*prctile(Emiss, 50);
LB = 0.001*prctile(Emiss, 5);
UB = 0.001*prctile(Emiss, 95);

p1 = boundedline(CFC11post.Years(2:end), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);

Emiss = CFC11post.emiss_m;
MED = 0.001*prctile(Emiss, 50);
LB = 0.001*prctile(Emiss, 5);
UB = 0.001*prctile(Emiss, 95);
hold on;
p2 = boundedline(CFC11post.Years(2:end), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0, 0, 0.7 ]);

Emiss = CFC11post.emiss_l;
MED = 0.001*prctile(Emiss, 50);
LB = 0.001*prctile(Emiss, 5);
UB = 0.001*prctile(Emiss, 95);
hold on;
p3 = boundedline(CFC11post.Years(2:end), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0.7 ]);


box on; title('CFC-11'); ylabel('Emissions [Gg yr^{-1}]');
xlim([1950,2020]); ylim([0,400]);
legend([p1, p2, p3], 'Short Bank', 'Medium Bank', 'Long Bank', 'Location','northwest'); 

subplot(2,4,2)
Emiss = CFC12post.emiss_s1 + CFC12post.emiss_s2;
MED = 0.001*prctile(Emiss, 50);
LB = 0.001*prctile(Emiss, 5);
UB = 0.001*prctile(Emiss, 95);

p1 = boundedline(CFC12post.Years(2:end), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);

Emiss = CFC12post.emiss_m;
MED = 0.001*prctile(Emiss, 50);
LB = 0.001*prctile(Emiss, 5);
UB = 0.001*prctile(Emiss, 95);
hold on; 
p2 = boundedline(CFC12post.Years(2:end), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0, 0, 0.7]);

Emiss = CFC12post.emiss_l;
MED = 0.001*prctile(Emiss, 50);
LB = 0.001*prctile(Emiss, 5);
UB = 0.001*prctile(Emiss, 95);
hold on; 
p3 = boundedline(CFC12post.Years(2:end), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0.7]);
legend([p1, p2, p3], 'Short Bank', 'Medium Bank', 'Long Bank', 'Location','northwest'); 
ylim([0,500]);
box on; title('CFC-12'); ylabel('Emissions [Gg yr^{-1}]');
xlim([1950,2020]); 

subplot(2,4,3)
Emiss = CFC113post.emiss_s;
MED = 0.001*prctile(Emiss, 50);
LB = 0.001*prctile(Emiss, 5);
UB = 0.001*prctile(Emiss, 95);
p1 = boundedline(CFC113post.Years+1, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);

shortemiss113 = MED; 
Emiss = CFC113post.emiss_l;
MED = 0.001*prctile(Emiss, 50);
LB = 0.001*prctile(Emiss, 5);
UB = 0.001*prctile(Emiss, 95);
hold on;
p2 = boundedline(CFC113post.Years+1, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0.7]);
longemiss113 = MED; 

Emiss = CFC113post.emiss_f;
MED = 0.001*prctile(Emiss, 50);
LB = 0.001*prctile(Emiss, 5);
UB = 0.001*prctile(Emiss, 95);
hold on;
p3 = boundedline(CFC113post.Years+1, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0.7, 0.7]);
legend([p1, p2, p3], 'Short Bank', 'Long Bank', 'Feedstock', 'Location','northwest'); 
FeedStockEmiss(1) = mean(MED(end-8:end)); 
feedstockemiss113 = MED;

box on; title('CFC-113'); ylabel('Emissions [Gg yr^{-1}]');
xlim([1950,2020]); 

subplot(2,4,4)
Emiss = CFC114post.emiss_s;
MED = prctile(Emiss, 50);
LB = prctile(Emiss, 5);
UB = prctile(Emiss, 95);

p1 = boundedline(CFC114post.Years(2:end-9), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);

Emiss = CFC114post.emiss_l;
MED = prctile(Emiss, 50);
LB = prctile(Emiss, 5);
UB = prctile(Emiss, 95);
hold on;
p3 = boundedline(CFC114post.Years(2:end-9), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0.7 ]);
legend([p1, p3], 'Short Bank', 'Long Bank', 'Location','northwest'); 

box on; title('CFC-114'); ylabel('Emissions [Gg yr^{-1}]');
xlim([1950,2020]); 


subplot(2,4,5)
Emiss = CFC115post.emiss_s;
MED = prctile(Emiss, 50);
LB = prctile(Emiss, 5);
UB = prctile(Emiss, 95);

p1 = boundedline(CFC115post.Years(2:end-9), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);

Emiss =CFC115post.emiss_l;
MED = prctile(Emiss, 50);
LB = prctile(Emiss, 5);
UB = prctile(Emiss, 95);
hold on;
p2 = boundedline(CFC115post.Years(2:end-9), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0.7 ]);
legend([p1, p2], 'Short Bank', 'Long Bank', 'Location','northwest'); 
box on; title('CFC-115'); ylabel('Emissions [Gg yr^{-1}]');
xlim([1950,2020]); 


subplot(2,4,6)
Emiss = HCFC22post.emiss_s;
MED = prctile(Emiss, 50);
LB = prctile(Emiss, 5);
UB = prctile(Emiss, 95);

p1 = boundedline(HCFC22post.Years, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);

Emiss = HCFC22post.emiss_m;
MED = prctile(Emiss, 50);
LB = prctile(Emiss, 5);
UB = prctile(Emiss, 95);
hold on;
p2 = boundedline(HCFC22post.Years, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0, 0, 0.7 ]);

Emiss = HCFC22post.emiss_l;
MED = prctile(Emiss, 50);
LB = prctile(Emiss, 5);
UB = prctile(Emiss, 95);
hold on;
p3 = boundedline(HCFC22post.Years, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0.7]);

Emiss = HCFC22post.emiss_f;
MED = prctile(Emiss, 50);
LB = prctile(Emiss, 5);
UB = prctile(Emiss, 95);
hold on;
p4 = boundedline(HCFC22post.Years, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0.7, 0.7]);
legend([p1, p2, p3, p4], 'Short Bank', 'Medium Bank','Long Bank','Feedstock', 'Location','northwest'); 
box on; title('HCFC-22'); ylabel('Emissions [Gg yr^{-1}]');
xlim([1950,2020]); 
FeedStockEmiss(2) = mean(MED(end-8:end)); 

subplot(2,4,7)
Emiss = HCFC141bpost.emiss_s;
MED = prctile(Emiss, 50);
LB = prctile(Emiss, 5);
UB = prctile(Emiss, 95);
p1 = boundedline(HCFC141bpost.Years(2:32), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);

Emiss = HCFC141bpost.emiss_m;
MED = prctile(Emiss, 50);
LB = prctile(Emiss, 5);
UB = prctile(Emiss, 95);
hold on;
p2 = boundedline(HCFC141bpost.Years(2:32), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0, 0, 0.7 ]);

Emiss = HCFC141bpost.emiss_l;
MED = prctile(Emiss, 50);
LB = prctile(Emiss, 5);
UB = prctile(Emiss, 95);
hold on;
p3 = boundedline(HCFC141bpost.Years(2:32), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0.7 ]);
legend([p1, p2, p3], 'Short Bank', 'Medium Bank','Long Bank','Location','northwest'); 

box on; title('HCFC-141b'); ylabel('Emissions [Gg yr^{-1}]');
xlim([1950,2020]); 


subplot(2,4,8)
Emiss = HCFC142bpost.emiss_s;
MED = prctile(Emiss, 50);
LB = prctile(Emiss, 5);
UB = prctile(Emiss, 95);

p1 = boundedline(HCFC142bpost.Years(1:39), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0 ]);

Emiss = HCFC142bpost.emiss_s + HCFC142bpost.emiss_m + HCFC142bpost.emiss_l + repmat(HCFC142bpost.de_f',1,39).*HCFC142bpost.prod_f(:,1:39);
MED = prctile(Emiss, 50);
LB = prctile(Emiss, 5);
UB = prctile(Emiss, 95);
hold on;
p2 = boundedline(HCFC142bpost.Years(1:39), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0, 0, 0.7 ]);

Emiss = HCFC142bpost.emiss_l;
MED = prctile(Emiss, 50);
LB = prctile(Emiss, 5);
UB = prctile(Emiss, 95);
hold on;
p3 = boundedline(HCFC142bpost.Years(1:39), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0, 0.7 ]);

Emiss = repmat(HCFC142bpost.de_f',1,39).*HCFC142bpost.prod_f(:,1:39);
MED = prctile(Emiss, 50);
LB = prctile(Emiss, 5);
UB = prctile(Emiss, 95);
hold on;
p4 = boundedline(HCFC142bpost.Years(1:39), MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.7, 0.7, 0.7 ]);
legend([p1, p2, p3, p4], 'Short Bank', 'Medium Bank','Long Bank', 'Feedstock','Location','northwest'); 
FeedStockEmiss(3) = mean(MED(38-8:38)); 

box on; title('HCFC-142b'); ylabel('Emissions [Gg yr^{-1}]');
xlim([1950,2020]); 


figure_width = 18; % in inches
figure_height = 7; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(fighandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

FigName = 'Fig04';
Figstr = strcat('Figures/',FigName,'.pdf'); 
print(gcf, '-dpdf', Figstr);

% Table 2
FeedStockImportance = ["Mass";"GWP";"ODP"];
CFC113FeedStockEmiss = FeedStockEmiss(1)*[1; GWP.cfc113; ODP.cfc113];
HCFC22FeedStockEmiss = FeedStockEmiss(2)*[1; GWP.hcfc22; ODP.hcfc22];
HCFC142bFeedStockEmiss = FeedStockEmiss(3)*[1; GWP.hcfc142b; ODP.hcfc142b];

feedstockSummary = table(FeedStockImportance, CFC113FeedStockEmiss, HCFC22FeedStockEmiss, HCFC142bFeedStockEmiss) 

%% Table 3: Comparing posterior production estimates to reported production values
load('Prior_Production.mat')
% Reported Production for all types subject to banking (i.e. no feedstocks
% are included)

% CFC-11
Reported = sum(CFC11ReportedProd.BankedProd);
Post_prod = CFC11post.prod_s + CFC11post.prod_m + CFC11post.prod_l;
Post_tmp = median(sum(Post_prod'));
UnderReportedProd(1) = (Post_tmp - Reported)/Reported;
Post_tmp = prctile(sum(Post_prod'),16);
UnderReportedlow(1) = (Post_tmp - Reported)/Reported;
Post_tmp = prctile(sum(Post_prod'),84);
UnderReportedhigh(1) = (Post_tmp - Reported)/Reported;

% CFC-12
Reported = sum(CFC12ReportedProd.BankedProd);
Post_prod = CFC12post.prod_s1 + CFC12post.prod_s2 +CFC12post.prod_m + CFC12post.prod_l;
Post_tmp = median(sum(Post_prod'));
UnderReportedProd(2) = (Post_tmp - Reported)/Reported;
Post_tmp = prctile(sum(Post_prod'),16);
UnderReportedlow(2) = (Post_tmp - Reported)/Reported;
Post_tmp = prctile(sum(Post_prod'),84);
UnderReportedhigh(2) = (Post_tmp - Reported)/Reported;

% CFC-113
Reported = sum(CFC113ReportedProd.BankedProd);
Post_prod = CFC113post.prod_s + CFC113post.prod_l;
Post_tmp = median(sum(Post_prod'));
UnderReportedProd(3) = (Post_tmp - Reported)/Reported;
Post_tmp = prctile(sum(Post_prod'),16);
UnderReportedlow(3) = (Post_tmp - Reported)/Reported;
Post_tmp = prctile(sum(Post_prod'),84);
UnderReportedhigh(3) = (Post_tmp - Reported)/Reported;

% CFC-114
Reported = sum(CFC114ReportedProd.BankedProd(1:end-11));
Post_prod = CFC114post.prod_s(:, 1:end-11) + CFC114post.prod_l(:, 1:end-11);
Post_tmp = median(sum(Post_prod'));
UnderReportedProd(4) = (Post_tmp - Reported)/Reported;
Post_tmp = prctile(sum(Post_prod'),16);
UnderReportedlow(4) = (Post_tmp - Reported)/Reported;
Post_tmp = prctile(sum(Post_prod'),84);
UnderReportedhigh(4) = (Post_tmp - Reported)/Reported;

% CFC-115
Reported = nansum(CFC115ReportedProd.BankedProd(1:end-11));
Post_prod = CFC115post.prod_s(:,1:end-11) + CFC115post.prod_l(:,1:end-11);
Post_tmp = median(sum(Post_prod'));
UnderReportedProd(5) = (Post_tmp - Reported)/Reported;
Post_tmp = prctile(sum(Post_prod'),16);
UnderReportedlow(5) = (Post_tmp - Reported)/Reported;
Post_tmp = prctile(sum(Post_prod'),84);
UnderReportedhigh(5) = (Post_tmp - Reported)/Reported;


% HCFC 22
Reported = sum(HCFC22ReportedProd.BankedProd);
Post_prod = HCFC22post.prod_s + HCFC22post.prod_m + HCFC22post.prod_l;
Post_tmp = median(sum(Post_prod'));
UnderReportedProd(6) = (Post_tmp - Reported)/Reported;
Post_tmp = prctile(sum(Post_prod'),16);
UnderReportedlow(6) = (Post_tmp - Reported)/Reported;
Post_tmp = prctile(sum(Post_prod'),84);
UnderReportedhigh(6) = (Post_tmp - Reported)/Reported;

% HCFC 141b
Reported = sum(HCFC141bReportedProd.BankedProd(1:31));
Post_prod = HCFC141bpost.prod_s + HCFC141bpost.prod_m + HCFC141bpost.prod_l;
Post_tmp = median(sum(Post_prod'));
UnderReportedProd(7) = (Post_tmp - Reported)/Reported;
Post_tmp = prctile(sum(Post_prod'),16);
UnderReportedlow(7) = (Post_tmp - Reported)/Reported;
Post_tmp = prctile(sum(Post_prod'),84);
UnderReportedhigh(7) = (Post_tmp - Reported)/Reported;

% HCFC 142b
Reported = sum(HCFC142bReportedProd.BankedProd(1:end-11));
Post_prod = HCFC142bpost.prod_s(:,1:end-11) + HCFC142bpost.prod_m(:,1:end-11) + HCFC142bpost.prod_l(:,1:end-11);
Post_tmp = median(sum(Post_prod'));
UnderReportedProd(8) = (Post_tmp - Reported)/Reported;
Post_tmp = prctile(sum(Post_prod'),16);
UnderReportedlow(8) = (Post_tmp - Reported)/Reported;
Post_tmp = prctile(sum(Post_prod'),84);
UnderReportedhigh(8) = (Post_tmp - Reported)/Reported;

% Halon1211
Reported = sum(h1211ReportedProd.BankedProd);
Post_prod = halon1211post.Prod;
Post_tmp = median(sum(Post_prod'));
UnderReportedProd(9) = (Post_tmp - Reported)/Reported;
Post_tmp = prctile(sum(Post_prod'),16);
UnderReportedlow(9) = (Post_tmp - Reported)/Reported;
Post_tmp = prctile(sum(Post_prod'),84);
UnderReportedhigh(9) = (Post_tmp - Reported)/Reported;

% Halon1301
Reported = sum(h1301ReportedProd.BankedProd);
Post_prod = halon1301post.Prod;
Post_tmp = median(sum(Post_prod'));
UnderReportedProd(10) = (Post_tmp - Reported)/Reported;
Post_tmp = prctile(sum(Post_prod'),16);
UnderReportedlow(10) = (Post_tmp - Reported)/Reported;
Post_tmp = prctile(sum(Post_prod'),84);
UnderReportedhigh(10) = (Post_tmp - Reported)/Reported;

Chemical = ["CFC-11";"CFC-12";"CFC-113";"CFC-114";"CFC-115";"HCFC-22";"HCFC-141b";"HCFC-142b";"halon1211";"halon1301"];

UnderReported_minus1sigma = UnderReportedlow';
UnderReported_plus1sigma = UnderReportedhigh';

feedstockSummary = table(Chemical, UnderReportedProd', UnderReported_minus1sigma, UnderReported_plus1sigma) 
