%% HCFC-22 Priors

load('HCFC-22/Input/HCFC22Production.mat'); 
totalAFEAS = HCFC22Prod.shortBank + HCFC22Prod.mediumBank + HCFC22Prod.longBank;
y1 = HCFC22Prod.feedstockYear(1); 
nonA5Prod = totalAFEAS;
UNEP = HCFC22Prod.nonfeedstock; 
ind = find(HCFC22Prod.AFEASyears > (y1-1));
ScaleFrac(1:ind(1)-1) = 1; 
ScaleFrac(ind) = 0.001*HCFC22Prod.nonfeedstock(1:length(ind))./totalAFEAS(ind);
ScaleFrac(ScaleFrac<1) = 1;

y1 = HCFC22Prod.AFEASyears(1);
y2 = HCFC22Prod.AFEASyears(end);
y3 = HCFC22Prod.feedstockYear(end);

% Corrected short banks up until 2003
shortBankProd = zeros(length(y1:y3),1); 
mediumBankProd = zeros(length(y1:y3),1); 
longBankProd = zeros(length(y1:y3),1); 

feedstockProd = zeros(length(y1:y3),1); 

shortBankProd(1:length(y1:y2)) = HCFC22Prod.shortBank.*ScaleFrac'; 
mediumBankProd(1:length(y1:y2)) = HCFC22Prod.mediumBank.*ScaleFrac'; 
longBankProd(1:length(y1:y2)) = HCFC22Prod.longBank.*ScaleFrac'; 

% Extended production out until 2019

y2 = HCFC22Prod.AFEASyears(end); 
indy2 = find(HCFC22Prod.feedstockYear == y2); 

tmpind1 = length(y1:y2) + 1; 

Scaling = (1/HCFC22Prod.nonfeedstock(indy2))*HCFC22Prod.nonfeedstock(indy2:end);
shortBankProd(tmpind1-1:end) = shortBankProd(tmpind1-1)*Scaling;
mediumBankProd(tmpind1-1:end) = mediumBankProd(tmpind1-1)*Scaling;
longBankProd(tmpind1-1:end) = longBankProd(tmpind1-1)*Scaling;
feedstockProd(end-length(HCFC22Prod.feedstock)+1:end) = 0.001*HCFC22Prod.feedstock;

HCFC22ReportedProd.BankedProd =  shortBankProd + mediumBankProd + longBankProd;
HCFC22ReportedProd.yrs = y1:y3;
save('Prior_Production.mat','HCFC22ReportedProd')

% Short Production Prior Samps
N = 100000;
Nyrs = length(y1:y3); 
rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs,Nyrs,1);
exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

shortProdsamps = (repmat(0.2*shortBankProd',N,1)).*logninv(Um,zeros(N,Nyrs),0.5*ones(N,Nyrs))+0.90*repmat(shortBankProd',N,1);

rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs,Nyrs,1);
exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

mediumProdsamps = (repmat(0.2*mediumBankProd',N,1)).*logninv(Um,zeros(N,Nyrs),0.5*ones(N,Nyrs))+0.90*repmat(mediumBankProd',N,1);

rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs,Nyrs,1);
exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

longProdsamps = (repmat(0.2*longBankProd',N,1)).*logninv(Um,zeros(N,Nyrs),0.5*ones(N,Nyrs))+0.90*repmat(longBankProd',N,1);

rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs,Nyrs,1);
exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

feedstocksamps(:, 1:Nyrs) = (repmat(0.2*feedstockProd',N,1)).*logninv(Um,zeros(N,Nyrs),0.5*ones(N,Nyrs))+0.90*repmat(feedstockProd',N,1);


RFshortBank.y2 = 2*0.17*betarnd(4,4,1,N); 
RFshortBank.y1 = 1 - RFshortBank.y2; 

Years = y1:y3;
Regimey1_m = [0.2, 0.1, 0.07, 0.06]; 
RegimeYrs = [1944, 1977, 1984, 1993, 2019]; 

samps_tmp = betarnd(6,6,1,N);
for ii = 1:4
    Nyrs = RegimeYrs(ii+1) - RegimeYrs(ii)+1;
    y1ind = find(Years == RegimeYrs(ii)); 
    samp_vals = 2*Regimey1_m(ii)*samps_tmp;
    RFmediumBank.y1(y1ind:y1ind+Nyrs-1,:) = ones(Nyrs, 1)*samp_vals;%0.07*DEu(5,:);
end

m = 0.1; v = (0.50*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFmediumBank.LT = lognrnd(mu,sigma,1,N);

%RFmediumBank.LT =  exp(-1/4.5) + 2*(1-exp(-1/4.5))*betarnd(2,2,N,1) - (1-exp(-1/4.5));

m = 0.02; v = (0.50*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFlongBank.LT = lognrnd(mu,sigma,1,N);

% I think I misread this before - this has been changed to be 2% of
% feedstock is emitted.  But is it 4.5%?  
m = 0.02; v = (0.50*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFfeedstock.RF = lognrnd(mu,sigma,1,N);

save('HCFC22Priors.mat','longProdsamps','shortProdsamps','mediumProdsamps','feedstocksamps','RFfeedstock','RFlongBank','RFshortBank','RFmediumBank');

%% HCFC 141b
clear all
load('HCFC141b/Input/HCFC141bProduction.mat'); 

prod_s = HCFC141bProduction.short; 
prod_m = HCFC141bProduction.medium; 
prod_l = HCFC141bProduction.long; 

total1 = prod_s + prod_m + prod_l; 
total2 = 0.001*HCFC141bProduction.nonAF+0.001*HCFC141bProduction.AF;

years1 = HCFC141bProduction.AFEASyears; 
years2 = HCFC141bProduction.UNEPyears;

yr1 = min(years2(2), years1(1)); 
yr2 = max(years2(end), years1(end)); 

Years = yr1:yr2;
TotalProd = zeros(length(Years),3); 

yr1 = find(Years == years1(1)); 
yr2 = find(Years == years1(end)); 
TotalProd(yr1:yr2,1) = total1; 

yr1 = find(Years == years2(2)); 
yr2 = find(Years == years2(end)); 
TotalProd(yr1:yr2,2) = total2(2:end); 

TotalProd(:,3) = max(TotalProd(:,1), TotalProd(:,2)); 

ind = Years < years1(1); 
prod_s = [zeros(sum(ind),1); prod_s];
prod_m = [zeros(sum(ind),1); prod_m];
prod_l = [zeros(sum(ind),1); prod_l];

prod_s(ind) = prod_s(sum(ind)+1); 
prod_m(ind) = prod_m(sum(ind)+1); 
prod_l(ind) = prod_l(sum(ind)+1); 

ind = Years > years1(end);
prod_s = [prod_s; zeros(sum(ind),1)];
prod_m = [prod_m; zeros(sum(ind),1)];
prod_l = [prod_l; zeros(sum(ind),1)];

ind2 = find(Years == years1(end));
prod_s(ind) = prod_s(ind2); 
prod_m(ind) = prod_m(ind2); 
prod_l(ind) = prod_l(ind2); 

tmp_Total = prod_s + prod_m + prod_l; 
Scaling = TotalProd(:,3)./tmp_Total;

prod_s = prod_s.*Scaling; 
prod_m = prod_m.*Scaling;
prod_l = prod_l.*Scaling; 

HCFC141bReportedProd.BankedProd =  prod_s + prod_m + prod_l;
HCFC141bReportedProd.yrs = Years;
save('Prior_Production.mat','HCFC141bReportedProd','-append')

% Short Production Prior Samps
N = 100000;
Nyrs = length(Years); 
rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs,Nyrs,1);
exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

shortProdsamps(:,1:Nyrs) = (repmat(0.2*prod_s',N,1)).*logninv(Um,zeros(N,Nyrs),0.5*ones(N,Nyrs))+0.90*repmat(prod_s',N,1);

rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs,Nyrs,1);
exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

mediumProdsamps(:,1:Nyrs) = (repmat(0.2*prod_m',N,1)).*logninv(Um,zeros(N,Nyrs),0.5*ones(N,Nyrs))+0.90*repmat(prod_m',N,1);

rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs,Nyrs,1);
exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

longProdsamps(:,1:Nyrs) = (repmat(0.2*prod_l',N,1)).*logninv(Um,zeros(N,Nyrs),0.5*ones(N,Nyrs))+0.90*repmat(prod_l',N,1);

RFshortBank.y2 = 2*0.17*betarnd(4,4,1,N); 
RFshortBank.y1 = 1 - RFshortBank.y2; 

RFmediumBank.y1 = 0.15 + 0.3*betarnd(4,4,1,N); 

RFmediumBank.LT =  2*(1 - exp(-1/4.5))*betarnd(4,4,N,1);

RFlongBank.y1 = 2*0.1*betarnd(4,4,1,N); 

m = 0.045; v = (0.20*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFlongBank.LT = lognrnd(mu,sigma,1,N);

save('HCFC141bPriors.mat','longProdsamps','shortProdsamps','mediumProdsamps','RFlongBank','RFshortBank','RFmediumBank','Years');

%% HCFC141b fugitive scenario - doesn't work because total production is too small

clear all
excess_prod = 61*10^3; % Montzka 
AddtoUB = [linspace(0,excess_prod,13),repmat(excess_prod,1,6)]; % Added for CFC-11, removed here.

load('HCFC141b/Input/HCFC141bProduction.mat'); 

prod_s = HCFC141bProduction.short; 
prod_m = HCFC141bProduction.medium; 
prod_l = HCFC141bProduction.long; 

total1 = prod_s + prod_m + prod_l; 
total2 = 0.001*HCFC141bProduction.nonAF+0.001*HCFC141bProduction.AF;

years1 = HCFC141bProduction.AFEASyears; 
years2 = HCFC141bProduction.UNEPyears;

yr1 = min(years2(2), years1(1)); 
yr2 = max(years2(end), years1(end)); 

Years = yr1:yr2;
TotalProd = zeros(length(Years),3); 

yr1 = find(Years == years1(1)); 
yr2 = find(Years == years1(end)); 
TotalProd(yr1:yr2,1) = total1; 

yr1 = find(Years == years2(2)); 
yr2 = find(Years == years2(end)); 
TotalProd(yr1:yr2,2) = total2(2:end); 

TotalProd(:,3) = max(TotalProd(:,1), TotalProd(:,2)); 

ind = Years < years1(1); 
prod_s = [zeros(sum(ind),1); prod_s];
prod_m = [zeros(sum(ind),1); prod_m];
prod_l = [zeros(sum(ind),1); prod_l];

prod_s(ind) = prod_s(sum(ind)+1); 
prod_m(ind) = prod_m(sum(ind)+1); 
prod_l(ind) = prod_l(sum(ind)+1); 

ind = Years > years1(end);
prod_s = [prod_s; zeros(sum(ind),1)];
prod_m = [prod_m; zeros(sum(ind),1)];
prod_l = [prod_l; zeros(sum(ind),1)];

ind2 = find(Years == years1(end));
prod_s(ind) = prod_s(ind2); 
prod_m(ind) = prod_m(ind2); 
prod_l(ind) = prod_l(ind2); 

tmp_Total = prod_s + prod_m + prod_l; 
Scaling = TotalProd(:,3)./tmp_Total;

prod_s = prod_s.*Scaling; 
prod_m = prod_m.*Scaling;
prod_l = prod_l.*Scaling; 



% Short Production Prior Samps
N = 100000;
Nyrs = length(Years); 
rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs,Nyrs,1);
exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

shortProdsamps(:,1:Nyrs) = (repmat(0.3*prod_s',N,1)).*logninv(Um,zeros(N,Nyrs),0.5*ones(N,Nyrs))+0.80*repmat(prod_s',N,1);

rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs,Nyrs,1);
exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

mediumProdsamps(:,1:Nyrs) = (repmat(0.3*prod_m',N,1)).*logninv(Um,zeros(N,Nyrs),0.5*ones(N,Nyrs))+0.80*repmat(prod_m',N,1);

rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs,Nyrs,1);
exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

longProdsamps(:,1:Nyrs) = (repmat(0.3*prod_l',N,1)).*logninv(Um,zeros(N,Nyrs),0.5*ones(N,Nyrs))+0.80*repmat(prod_l',N,1);

RFshortBank.y2 = 2*0.17*betarnd(2,2,1,N); 
RFshortBank.y1 = 1 - RFshortBank.y2; 

RFmediumBank.y1 = 0.15 + 0.3*betarnd(2,2,1,N); 

RFmediumBank.LT =  2*(1 - exp(-1/4.5))*betarnd(2,2,N,1);

RFlongBank.y1 = 2*0.1*betarnd(2,2,1,N); 

m = 0.045; v = (0.50*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFlongBank.LT = lognrnd(mu,sigma,1,N);

save('HCFC141bPriors.mat','longProdsamps','shortProdsamps','mediumProdsamps','RFlongBank','RFshortBank','RFmediumBank','Years');


%% HCFC 142b 
clear all
load('HCFC142b/Input/HCFC142bProduction.mat'); 
load('HCFC142b/Input/HCFC142bFeedstock.mat'); 

prod_s = HCFC142bProduction.short; 
prod_m = HCFC142bProduction.medium; 
prod_l = HCFC142bProduction.long; 
prod_f = Feedstock.production; 

total1 = prod_s + prod_m + prod_l; 
total2 = 0.001*HCFC142bProduction.nonA5+0.001*HCFC142bProduction.A5;

years1 = HCFC142bProduction.AFEASyears; 
years2 = HCFC142bProduction.UNEPyears;

yr1 = min(years2(2), years1(1)); 
yr2 = max(years2(end), years1(end)); 

Years = yr1:yr2;
TotalProd = zeros(length(Years),3); 

yr1 = find(Years == years1(1)); 
yr2 = find(Years == years1(end)); 
TotalProd(yr1:yr2,1) = total1; 

yr1 = find(Years == years2(2)); 
yr2 = find(Years == years2(end)); 
TotalProd(yr1:yr2,2) = total2(2:end); 

TotalProd(:,3) = max(TotalProd(:,1), TotalProd(:,2)); 

ind = Years < years1(1); 
prod_s = [zeros(sum(ind),1); prod_s];
prod_m = [zeros(sum(ind),1); prod_m];
prod_l = [zeros(sum(ind),1); prod_l];

prod_s(ind) = prod_s(sum(ind)+1); 
prod_m(ind) = prod_m(sum(ind)+1); 
prod_l(ind) = prod_l(sum(ind)+1); 

ind = Years > years1(end);
prod_s = [prod_s; zeros(sum(ind),1)];
prod_m = [prod_m; zeros(sum(ind),1)];
prod_l = [prod_l; zeros(sum(ind),1)];

ind2 = find(Years == years1(end));
prod_s(ind) = prod_s(ind2); 
prod_m(ind) = prod_m(ind2); 
prod_l(ind) = prod_l(ind2); 

tmp_Total = prod_s + prod_m + prod_l; 
Scaling = TotalProd(:,3)./tmp_Total;

prod_s = prod_s.*Scaling; 
prod_m = prod_m.*Scaling;
prod_l = prod_l.*Scaling; 

HCFC142bReportedProd.BankedProd =  prod_s + prod_m + prod_l;
HCFC142bReportedProd.yrs = Years;
save('Prior_Production.mat','HCFC142bReportedProd','-append')

prod_f = zeros(size(prod_s)); 
ind1 = find(Years == Feedstock.years(1)); 
ind1b = find(Feedstock.years == Years(end)); 
prod_f(ind1:end) = 0.001*Feedstock.production(1:ind1b); 

% Short Production Prior Samps
N = 100000;
Nyrs = length(Years); 
rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs,Nyrs,1);
exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

shortProdsamps(:,1:Nyrs) = (repmat(0.3*prod_s',N,1)).*logninv(Um,zeros(N,Nyrs),0.5*ones(N,Nyrs))+0.90*repmat(prod_s',N,1);

rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs,Nyrs,1);
exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

mediumProdsamps(:,1:Nyrs) = (repmat(0.3*prod_m',N,1)).*logninv(Um,zeros(N,Nyrs),0.5*ones(N,Nyrs))+0.90*repmat(prod_m',N,1);

rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs,Nyrs,1);
exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

longProdsamps(:,1:Nyrs) = (repmat(0.3*prod_l',N,1)).*logninv(Um,zeros(N,Nyrs),0.5*ones(N,Nyrs))+0.90*repmat(prod_l',N,1);

rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs,Nyrs,1);
exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

Feedstocksamps(:,1:Nyrs) = (repmat(0.3*prod_f',N,1)).*logninv(Um,zeros(N,Nyrs),0.5*ones(N,Nyrs))+0.90*repmat(prod_f',N,1);


RFshortBank.y2 = 2*0.17*betarnd(4,4,1,N); 
RFshortBank.y1 = 1 - RFshortBank.y2; 

RFmediumBank.y1 = 0.15 + 0.3*betarnd(4,4,1,N); 

RFmediumBank.LT =  2*(1 - exp(-1/4.5))*betarnd(4,4,N,1);

RFlongBank.y1 = 0.15 + 2*0.175*betarnd(4,4,1,N); 

m = 0.03; v = (0.20*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFlongBank.LT = lognrnd(mu,sigma,1,N);

m = 0.02; v = (0.50*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFfeedstock.RF = lognrnd(mu,sigma,1,N);

save('HCFC142bPriors.mat','longProdsamps','shortProdsamps','mediumProdsamps','Feedstocksamps','RFfeedstock', 'RFlongBank','RFshortBank','RFmediumBank','Years');


%% CFC-114
clear all
load('CFC114/Input/CFC114Production.mat'); 

prod_s = CFC114Production.short; 
prod_l = CFC114Production.long; 

total1 = prod_s + prod_l; 
total2 = 0.001*CFC114Production.nonA5+0.001*CFC114Production.A5;

years1 = CFC114Production.AFEASyears; 
years2 = CFC114Production.UNEPyears;

yr1 = min(years2(2), years1(1)); 
yr2 = max(years2(end), years1(end)); 

Years = yr1:yr2;
TotalProd = zeros(length(Years),3); 

yr1 = find(Years == years1(1)); 
yr2 = find(Years == years1(end)); 
TotalProd(yr1:yr2,1) = total1; 

yr1 = find(Years == years2(2)); 
yr2 = find(Years == years2(end)); 
TotalProd(yr1:yr2,2) = total2(2:end); 

TotalProd(:,3) = max(TotalProd(:,1), TotalProd(:,2)); 

ind = Years < years1(1); 
prod_s = [zeros(sum(ind),1); prod_s];
prod_l = [zeros(sum(ind),1); prod_l];

prod_s(ind) = prod_s(sum(ind)+1); 
prod_l(ind) = prod_l(sum(ind)+1); 

ind = Years > years1(end);
prod_s = [prod_s; zeros(sum(ind),1)];
prod_l = [prod_l; zeros(sum(ind),1)];

ind2 = find(Years == years1(end));
prod_s(ind) = prod_s(ind2); 
prod_l(ind) = prod_l(ind2); 

tmp_Total = prod_s + prod_l; 
Scaling = TotalProd(:,3)./tmp_Total;

prod_s = prod_s.*Scaling; 
prod_l = prod_l.*Scaling; 

CFC114ReportedProd.BankedProd =  prod_s + prod_l;
CFC114ReportedProd.yrs = Years;
save('Prior_Production.mat','CFC114ReportedProd','-append')

% Short Production Prior Samps
N = 100000;
Nyrs = length(Years); 
rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs,Nyrs,1);
exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

shortProdsamps(:,1:Nyrs) = (repmat(0.2*prod_s',N,1)).*logninv(Um,zeros(N,Nyrs),0.5*ones(N,Nyrs))+0.90*repmat(prod_s',N,1);

rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs,Nyrs,1);
exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

longProdsamps(:,1:Nyrs) = (repmat(0.2*prod_l',N,1)).*logninv(Um,zeros(N,Nyrs),0.5*ones(N,Nyrs))+0.90*repmat(prod_l',N,1);


RFshortBank.y1 = betarnd(8,8,1,N); 
RFshortBank.y2 = 1 - RFshortBank.y1; 

m = 0.02; v = (m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFlongBank.y1 = lognrnd(mu,sigma,1,N);

m = 1 - exp(-1/20); v = (0.50*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFlongBank.LT = lognrnd(mu,sigma,1,N);

save('CFC114Priors.mat','longProdsamps','shortProdsamps','RFlongBank','RFshortBank','Years');

%% CFC-115
clear all
load('CFC115/Input/CFC115Production.mat'); 

prod_s = CFC115Production.short; 
prod_l = CFC115Production.long; 

total1 = prod_s + prod_l; 
total2 = 0.001*CFC115Production.nonA5+0.001*CFC115Production.A5;

years1 = CFC115Production.AFEASyears; 
years2 = CFC115Production.UNEPyears;

yr1 = min(years2(2), years1(1)); 
yr2 = max(years2(end), years1(end)); 

Years = yr1:yr2;
TotalProd = zeros(length(Years),3); 

yr1 = find(Years == years1(1)); 
yr2 = find(Years == years1(end)); 
TotalProd(yr1:yr2,1) = total1; 

yr1 = find(Years == years2(2)); 
yr2 = find(Years == years2(end)); 
TotalProd(yr1:yr2,2) = total2(2:end); 

TotalProd(:,3) = max(TotalProd(:,1), TotalProd(:,2)); 

ind = Years < years1(1); 
prod_s = [zeros(sum(ind),1); prod_s];
prod_l = [zeros(sum(ind),1); prod_l];

prod_s(ind) = prod_s(sum(ind)+1); 
prod_l(ind) = prod_l(sum(ind)+1); 

ind = Years > years1(end);
prod_s = [prod_s; zeros(sum(ind),1)];
prod_l = [prod_l; zeros(sum(ind),1)];

ind2 = find(Years == years1(end));
prod_s(ind) = prod_s(ind2); 
prod_l(ind) = prod_l(ind2); 

tmp_Total = prod_s + prod_l; 
Scaling = TotalProd(:,3)./tmp_Total;

prod_s = prod_s.*Scaling; 
prod_l = prod_l.*Scaling; 

CFC115ReportedProd.BankedProd =  prod_s + prod_l;
CFC115ReportedProd.yrs = Years;
save('Prior_Production.mat','CFC115ReportedProd','-append')

% Short Production Prior Samps
N = 100000;
Nyrs = length(Years); 
rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs,Nyrs,1);
exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

shortProdsamps(:,1:Nyrs) = (repmat(0.2*prod_s',N,1)).*logninv(Um,zeros(N,Nyrs),0.5*ones(N,Nyrs))+0.90*repmat(prod_s',N,1);

rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs,Nyrs,1);
exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

longProdsamps(:,1:Nyrs) = (repmat(0.2*prod_l',N,1)).*logninv(Um,zeros(N,Nyrs),0.5*ones(N,Nyrs))+0.90*repmat(prod_l',N,1);


RFshortBank.y1 = betarnd(8,8,1,N); 
RFshortBank.y2 = 1 - RFshortBank.y1; 


RFlongBank.y1 = 2*0.07*betarnd(4,4,1,N); 

m = 1 - exp(-1/10); v = (0.50*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFlongBank.LT = lognrnd(mu,sigma,1,N);

save('CFC115Priors.mat','longProdsamps','shortProdsamps','RFlongBank','RFshortBank','Years');

%% Halon

% According to McCulloch: Uncertainties
% Halon-1211:
%  prior to 1980: 13% for 1211 
% 1980-1983:  +1.5%, -0.5%
% 1983 - onwards: +2.5%, -1.5%

% Halon-1301
% Up until 1972: 8%
% After 1972: +1.5%, -0.5%
% Prior to 1972:  +/- 8% 
clear all
load('halons/Input/halonProd.mat'); 
Years = halon.productionYears(1):2020; 
N = 100000;
Nyrs = length(Years); 

Prod_1211 = zeros(size(Years)); 
Prod_1301 = zeros(size(Years)); 

Prod_1211(1:length(halon.production1211)) = halon.production1211;
Prod_1301(1:length(halon.production1301)) = halon.production1301;

h1211ReportedProd.BankedProd = Prod_1211;
h1211ReportedProd.yrs = Years;
save('Prior_Production.mat','h1211ReportedProd','-append')

h1301ReportedProd.BankedProd = Prod_1301;
h1301ReportedProd.yrs = Years;
save('Prior_Production.mat','h1301ReportedProd','-append')

% Distribution for 1211
Nyrs1 = length(halon.productionYears(1):1980);

rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs1,Nyrs1,1);
exp_val = repmat(abs(repmat([1:Nyrs1],Nyrs1,1)-repmat([1:Nyrs1]',1,Nyrs1)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs1), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

Prodsamps1211(:,1:Nyrs1) = (repmat(0.5*Prod_1211(1:Nyrs1),N,1)).*logninv(Um,zeros(N,Nyrs1),0.5*ones(N,Nyrs1))+0.80*repmat(Prod_1211(1:Nyrs1),N,1);

ind1 = Nyrs1+1; 

Nyrs1 = length(1981:1990);

rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs1,Nyrs1,1);
exp_val = repmat(abs(repmat([1:Nyrs1],Nyrs1,1)-repmat([1:Nyrs1]',1,Nyrs1)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs1), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

Prodsamps1211(:,ind1:ind1+Nyrs1-1) = (repmat(0.3*Prod_1211(ind1:ind1+Nyrs1-1),N,1)).*logninv(Um,zeros(N,Nyrs1),0.5*ones(N,Nyrs1))+0.95*repmat(Prod_1211(ind1:ind1+Nyrs1-1),N,1);

ind1 = ind1+Nyrs1; 

Nyrs1 = length(1991:2020);

rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs1,Nyrs1,1);
exp_val = repmat(abs(repmat([1:Nyrs1],Nyrs1,1)-repmat([1:Nyrs1]',1,Nyrs1)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs1), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

Prodsamps1211(:,ind1:ind1+Nyrs1-1) = (repmat(0.3*Prod_1211(ind1:ind1+Nyrs1-1),N,1)).*logninv(Um,zeros(N,Nyrs1),0.5*ones(N,Nyrs1))+0.95*repmat(Prod_1211(ind1:ind1+Nyrs1-1),N,1);

% Distribution for 1301
Nyrs1 = length(halon.productionYears(1):1972);

rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs1,Nyrs1,1);
exp_val = repmat(abs(repmat([1:Nyrs1],Nyrs1,1)-repmat([1:Nyrs1]',1,Nyrs1)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs1), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

Prodsamps1301(:,1:Nyrs1) = (repmat(0.5*Prod_1301(1:Nyrs1),N,1)).*logninv(Um,zeros(N,Nyrs1),0.5*ones(N,Nyrs1))+0.80*repmat(Prod_1301(1:Nyrs1),N,1);

ind1 = Nyrs1+1; 

Nyrs1 = length(1973:1990);

rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs1,Nyrs1,1);
exp_val = repmat(abs(repmat([1:Nyrs1],Nyrs1,1)-repmat([1:Nyrs1]',1,Nyrs1)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs1), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

Prodsamps1301(:,ind1:ind1+Nyrs1-1) = (repmat(0.3*Prod_1301(ind1:ind1+Nyrs1-1),N,1)).*logninv(Um,zeros(N,Nyrs1),0.5*ones(N,Nyrs1))+0.95*repmat(Prod_1301(ind1:ind1+Nyrs1-1),N,1);

ind1 = ind1+Nyrs1; 

Nyrs1 = length(1991:2020);

rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs1,Nyrs1,1);
exp_val = repmat(abs(repmat([1:Nyrs1],Nyrs1,1)-repmat([1:Nyrs1]',1,Nyrs1)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs1), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

Prodsamps1301(:,ind1:ind1+Nyrs1-1) = (repmat(0.3*Prod_1301(ind1:ind1+Nyrs1-1),N,1)).*logninv(Um,zeros(N,Nyrs1),0.5*ones(N,Nyrs1))+0.95*repmat(Prod_1301(ind1:ind1+Nyrs1-1),N,1);

RF1211.DEpre1990 = 0.4*betarnd(2,2,1,N); 
RF1211.RFpre1990 = 0.12*betarnd(2,2,1,N);

RF1211.DEpost1990 = 0.1*betarnd(2,2,1,N); 
RF1211.RFpost1990 = 0.08*betarnd(2,2,1,N); % Safeguarding the ozone layer chapter 9 says 4%

RF1301.DEpre1990 = 0.4*betarnd(2,2,1,N); 1
RF1301.RFpre1990 = 0.12*betarnd(2,2,1,N);

RF1301.DEpost1990 = 0.1*betarnd(2,2,1,N); 
RF1301.RFpost1990 = 0.08*betarnd(2,2,1,N);

save('halonPriors.mat','RF1211','RF1301','Prodsamps1301','Prodsamps1211','Years'); 

%% CFC-11
% Production priors

clear all
HomeDir = '/Users/meganlickley/Dropbox (MIT)/Research/Lifetimes/Code';
cd(HomeDir)

y1 = 1955; % Year Start Date
yswitch = 1988;
y2 = 2018; % Year End Date
Sim_yrs = y1:y2;
N = 5*10^5;

Prod_reported = NaN(size(Sim_yrs)); % Production is in Tonnes
load('CFC11/Input/WMO2002.mat') % Use until 1988 - WMO data is the adjusted AFEAS data
ytmp1 = find(year==y1); 
ytmp2 = find(year==yswitch);
Prod_reported(1:yswitch+1-y1) = Production(ytmp1:ytmp2);

load('CFC11/Input/a5_na5.mat') % Use Article 5 and non A5 prod total from 1989 onwards
ytmp1 = find(a5_na5(:,1) == yswitch+1);
ytmp2 = find(Sim_yrs == a5_na5(end,1));
Prod_reported(1989-y1+1:ytmp2) = sum(a5_na5(ytmp1:end,[2:3]),2);

% Total Production, Scale AFEAS production by sector
Prod_reported(ytmp2:end) = sum(a5_na5(end,[2:3]),2);

str = strcat(HomeDir,'/CFC11/Input/AFEAS_cfc11production.mat');
load(str,'closecell','nonhermetic','open_aero','yr')

ind = find(yr == y1); 
add_yrs = [yr(end)+1:Sim_yrs(end)]'; 

tmp_s = open_aero(ind:end) - open_aero(ind-1:end-1);
tmp_m = nonhermetic(ind:end) - nonhermetic(ind-1:end-1);
tmp_l = closecell(ind:end) - closecell(ind-1:end-1);

tmp_s = [tmp_s; tmp_s(end)*ones(size(add_yrs))];
tmp_m = [tmp_m; tmp_m(end)*ones(size(add_yrs))];
tmp_l = [tmp_l; tmp_l(end)*ones(size(add_yrs))];
tmp_total = tmp_s + tmp_m + tmp_l;

prod_s = (tmp_s./(tmp_total)).*Prod_reported'; 
prod_m = (tmp_m./(tmp_total)).*Prod_reported'; 
prod_l = (tmp_l./(tmp_total)).*Prod_reported'; 

CFC11ReportedProd.BankedProd = Prod_reported;
CFC11ReportedProd.yrs = y1:y2;
save('Prior_Production.mat','CFC11ReportedProd','-append')

% creating production distribution
for ii = 1:3
    
    switch ii
        case 1 
            Prod_tmp = prod_s';
        case 2
            Prod_tmp = prod_m'; 
        case 3
            Prod_tmp = prod_l';    
    end
    
    Nyrs = length(y1:yswitch);
    rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
    rho_tmp = repmat(rho,Nyrs,Nyrs,1);
    exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N);
    Rhom = rho_tmp.^exp_val;
    clear rho_tmp exp_val
    Zm = mvnrnd(zeros(N,Nyrs), Rhom, N);
    clear Rhom
    Um = normcdf(Zm,0,1);
    clear Zm
    Prod_samps(:,1:Nyrs) = (repmat(0.2*Prod_tmp(1:Nyrs),N,1)).*logninv(Um,zeros(N,Nyrs),0.5*ones(N,Nyrs))+0.95*repmat(Prod_tmp(1:Nyrs),N,1);
    clear Um

    Nyrs2 = length(Sim_yrs) - Nyrs; 
    rho_tmp = repmat(rho,Nyrs2,Nyrs2,1);

    exp_val = repmat(abs(repmat([1:Nyrs2],Nyrs2,1)-repmat([1:Nyrs2]',1,Nyrs2)),1,1,N);
    Rhom = rho_tmp.^exp_val;
    clear rho_tmp exp_val
    Zm = mvnrnd(zeros(N,Nyrs2), Rhom, N);
    clear Rhom
    Um = normcdf(Zm,0,1);
    clear Zm
    Prod_samps(:,Nyrs+1:length(Sim_yrs)) = (repmat(0.2*Prod_tmp(Nyrs+1:end),N,1)).*logninv(Um,zeros(N,Nyrs2),0.5*ones(N,Nyrs2))+0.95*repmat(Prod_tmp(Nyrs+1:end),N,1);
    clear Um

    switch ii
        case 1 
            shortProdsamps = Prod_samps;
        case 2
            mediumProdsamps = Prod_samps;
        case 3
            longProdsamps = Prod_samps;   
    end
end

DEshort = betarnd(8,6,1,N);
RFshortBank = 1-DEshort;

m = 0.0366; v = (0.50*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFlongBank = lognrnd(mu,sigma,1,N);

m = 0.07; v = (0.50*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
DEmedium= lognrnd(mu,sigma,1,N);

m = 10; v = (0.20*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFmediumLT = lognrnd(mu,sigma,1,N);

Years = 1955:2018;
save('CFC11Priors.mat','Years','longProdsamps','shortProdsamps','mediumProdsamps','RFlongBank','RFshortBank','RFmediumLT','DEshort','DEmedium','Prod_reported');


% CFC11 Fugitive Production Scenario

Nyrs3 = length(2000:y2);
excess_prod = 61*10^3; % Montzka 
AddtoUB = [linspace(0,excess_prod,13),repmat(excess_prod,1,6)];
Prod_fug = Prod_reported;
Prod_fug(end-Nyrs3+1:end) = Prod_reported(end-Nyrs3+1:end)+AddtoUB;

prod_s = (tmp_s./(tmp_total)).*Prod_fug'; 
prod_m = (tmp_m./(tmp_total)).*Prod_fug'; 
prod_l = (tmp_l./(tmp_total)).*Prod_fug'; 
prod_tot = prod_s + prod_m + prod_l;



% creating production distribution
for ii = 1:3
    
    rho_tmp = repmat(rho,Nyrs3,Nyrs3,1);
    exp_val = repmat(abs(repmat([1:Nyrs3],Nyrs3,1)-repmat([1:Nyrs3]',1,Nyrs3)),1,1,N);
    Rhom = rho_tmp.^exp_val;
    clear rho_tmp exp_val
    Zm = mvnrnd(zeros(N,Nyrs3), Rhom, N);
    clear Rhom
    Um = normcdf(Zm,0,1);
    switch ii
        case 1 
            Prod_tmp = prod_s';
            Prod_samps = shortProdsamps; 
        case 2
            Prod_tmp = prod_m'; 
            Prod_samps = mediumProdsamps; 
        case 3
            Prod_tmp = prod_l'; 
            Prod_samps = longProdsamps; 
    end
    
    Prod_samps(:,end-Nyrs3+1:end) = (repmat(Prod_tmp(end-Nyrs3+1:end),N,1)).*logninv(Um,zeros(N,Nyrs3),0.5*ones(N,Nyrs3));
    
    switch ii
        case 1 
            shortProdsamps = Prod_samps;
        case 2
            mediumProdsamps = Prod_samps;
        case 3
            longProdsamps = Prod_samps;   
    end
end

save('CFC11Priors.mat','Years','longProdsamps','shortProdsamps','mediumProdsamps','RFlongBank','RFshortBank','RFmediumLT','DEshort','DEmedium');


%% CFC-12 

clear all
HomeDir = '/Users/meganlickley/Dropbox (MIT)/Research/Lifetimes/Code';
cd(HomeDir)

y1 = 1955; % Year Start Date
yswitch = 1988;
y2 = 2018; % Year End Date
Sim_yrs = y1:y2;
N = 5*10^5;

Prod_reported = NaN(size(Sim_yrs)); % Production is in Tonnes
load('CFC12/Input/WMO2002.mat') % Use until 1988 - WMO data is the adjusted AFEAS data
ytmp1 = find(year==y1); 
ytmp2 = find(year==yswitch);
Prod_reported(1:yswitch+1-y1) = Production(ytmp1:ytmp2);

load('CFC12/Input/a5_na5.mat') % Use Article 5 and non A5 prod total from 1989 onwards
ytmp1 = find(a5_na5(:,1) == yswitch+1);
ytmp2 = find(Sim_yrs == a5_na5(end,1));
Prod_reported(1989-y1+1:ytmp2) = sum(a5_na5(ytmp1:end,[2:3]),2);
Prod_reported(ytmp2:end) = sum(a5_na5(end,[2:3]),2);

str = strcat(HomeDir,'/CFC12/Input/AFEAS_cfc12production.mat');
load(str,'closecell','nonhermetic','open_aero','refrigeration','yr')


ind = find(yr == y1); 
add_yrs = [yr(end)+1:Sim_yrs(end)]'; 

tmp_s1 = closecell(ind:end) - closecell(ind-1:end-1);
tmp_s2 = open_aero(ind:end) - open_aero(ind-1:end-1);
tmp_m = nonhermetic(ind:end) - nonhermetic(ind-1:end-1);
tmp_l = refrigeration(ind:end) - refrigeration(ind-1:end-1);

tmp_s1 = [tmp_s1; tmp_s1(end)*ones(size(add_yrs))];
tmp_s2 = [tmp_s2; tmp_s2(end)*ones(size(add_yrs))];
tmp_m = [tmp_m; tmp_m(end)*ones(size(add_yrs))];
tmp_l = [tmp_l; tmp_l(end)*ones(size(add_yrs))];
tmp_total = tmp_s1 + tmp_s2 + tmp_m + tmp_l;

prod_s1 = (tmp_s1./(tmp_total)).*Prod_reported'; 
prod_s2 = (tmp_s2./(tmp_total)).*Prod_reported'; 
prod_m = (tmp_m./(tmp_total)).*Prod_reported'; 
prod_l = (tmp_l./(tmp_total)).*Prod_reported'; 

CFC12ReportedProd.BankedProd = Prod_reported;
CFC12ReportedProd.yrs = y1:y2;
save('Prior_Production.mat','CFC12ReportedProd','-append')

% creating production distribution
for ii = 1:4
    
    switch ii
        case 1 
            Prod_tmp = prod_s1';
        case 2
            Prod_tmp = prod_s2';            
        case 3
            Prod_tmp = prod_m'; 
        case 4
            Prod_tmp = prod_l';    
    end
    
    Nyrs = length(y1:yswitch);
    rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
    rho_tmp = repmat(rho,Nyrs,Nyrs,1);
    exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N);
    Rhom = rho_tmp.^exp_val;
    clear rho_tmp exp_val
    Zm = mvnrnd(zeros(N,Nyrs), Rhom, N);
    clear Rhom
    Um = normcdf(Zm,0,1);
    clear Zm
    Prod_samps(:,1:Nyrs) = (repmat(0.2*Prod_tmp(1:Nyrs),N,1)).*logninv(Um,zeros(N,Nyrs),0.5*ones(N,Nyrs))+0.95*repmat(Prod_tmp(1:Nyrs),N,1);
    clear Um

    Nyrs2 = length(Sim_yrs) - Nyrs; 
    rho_tmp = repmat(rho,Nyrs2,Nyrs2,1);

    exp_val = repmat(abs(repmat([1:Nyrs2],Nyrs2,1)-repmat([1:Nyrs2]',1,Nyrs2)),1,1,N);
    Rhom = rho_tmp.^exp_val;
    clear rho_tmp exp_val
    Zm = mvnrnd(zeros(N,Nyrs2), Rhom, N);
    clear Rhom
    Um = normcdf(Zm,0,1);
    clear Zm
    Prod_samps(:,Nyrs+1:length(Sim_yrs)) = (repmat(0.1*Prod_tmp(Nyrs+1:end),N,1)).*logninv(Um,zeros(N,Nyrs2),0.5*ones(N,Nyrs2))+0.95*repmat(Prod_tmp(Nyrs+1:end),N,1);
    clear Um

    switch ii
        case 1 
            short1Prodsamps = Prod_samps;
        case 2
            short2Prodsamps = Prod_samps;            
        case 3
            mediumProdsamps = Prod_samps;
        case 4
            longProdsamps = Prod_samps;   
    end
end

% model aerosol and opencell foam together bc of AFEAS grouping
DEshort2 = betarnd(2,2,1,N);
RFshort2Bank = 1-DEshort2;

DEshort1 = betarnd(2,2,1,N);
RFshort1Bank = 1 - DEshort1;

m = 0.07; v = (0.5*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
DEmedium = lognrnd(mu,sigma,1,N);%0.07*DEu(5,:);

m = 10; v = (0.20*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFmediumLT = lognrnd(mu,sigma,1,N);

m = 0.02; v = (0.5*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
DElong = lognrnd(mu,sigma,1,N); %0.02*DEu(7,:);

m = 20; v = (0.20*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFlongLT = lognrnd(mu,sigma,1,N);

Years = 1955:2018;
save('CFC12Priors.mat','Years','longProdsamps','short1Prodsamps','short2Prodsamps','mediumProdsamps','RFlongLT','RFshort1Bank','RFshort2Bank','RFmediumLT','DEshort2','DEshort1','DEmedium','DElong','Prod_reported');

%% CFC-113 

clear all
HomeDir = '/Users/meganlickley/Dropbox (MIT)/Research/Lifetimes/Code';
cd(HomeDir)

y1 = 1955; % Year Start Date
yswitch = 1988;
y2 = 2018; % Year End Date
Sim_yrs = y1:y2;
N = 5*10^5;

Prod_reported = NaN(size(Sim_yrs)); % Production is in Tonnes
load('CFC113/Input/WMO2002.mat') % Use until 1988 - WMO data is the adjusted AFEAS data
ytmp1 = find(year==y1); 
ytmp2 = find(year==yswitch);
Prod_reported(1:yswitch+1-y1) = Production(ytmp1:ytmp2);

load('CFC113/Input/a5_na5.mat') % Use Article 5 and non A5 prod total from 1989 onwards
ytmp1 = find(a5_na5(:,1) == yswitch+1);
ytmp2 = find(Sim_yrs == a5_na5(end,1));
Prod_reported(1989-y1+1:ytmp2) = sum(a5_na5(ytmp1:end,[2:3]),2);
Prod_reported(ytmp2:end) = sum(a5_na5(end,[2:3]),2);

str = strcat(HomeDir,'/CFC113/Input/AFEAS_cfc113production.mat');
load(str,'longbank','shortbank','yr'); % in thousands of tonnes

ind = find(yr == y1); 
add_yrs = [yr(end)+1:Sim_yrs(end)]'; 

tmp_s = shortbank(ind:end) - shortbank(ind-1:end-1);
tmp_l = longbank(ind:end) - longbank(ind-1:end-1);

tmp_s = [tmp_s; tmp_s(end)*ones(size(add_yrs))];
tmp_l = [tmp_l; tmp_l(end)*ones(size(add_yrs))];
tmp_total = tmp_s + tmp_l;

prod_s = (tmp_s./(tmp_total)).*Prod_reported'; 
prod_l = (tmp_l./(tmp_total)).*Prod_reported'; 

CFC113ReportedProd.BankedProd = Prod_reported;
CFC113ReportedProd.yrs = y1:y2;
save('Prior_Production.mat','CFC113ReportedProd','-append')

str = strcat(HomeDir,'/CFC113/Input/Feedstock.mat'); % in tonnes

load(str)

prod_f = zeros(size(prod_s)); 
ind1 = find(Sim_yrs == Feedstock.years(1)); 
prod_f(ind1:end) = 0.001*Feedstock.tonnes(1:end-1); 


% creating production distribution
for ii = 1:2
    
    switch ii
        case 1 
            Prod_tmp = prod_s';
        case 2
            Prod_tmp = prod_l';    
    end
    
   % creating production distribution
    Nyrs = length(y1:yswitch);
    rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
    rho_tmp = repmat(rho,Nyrs,Nyrs,1);
    exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N);
    Rhom = rho_tmp.^exp_val;
    clear rho_tmp exp_val
    Zm = mvnrnd(zeros(N,Nyrs), Rhom, N);
    clear Rhom
    Um = normcdf(Zm,0,1);
    clear Zm
    Prod_samps(:,1:Nyrs) = (repmat(0.4*Prod_tmp(1:Nyrs),N,1)).*logninv(Um,zeros(N,Nyrs),0.5*ones(N,Nyrs))+0.70*repmat(Prod_tmp(1:Nyrs),N,1);
    clear Um

    Nyrs2 = length(Sim_yrs) - Nyrs; 
    rho_tmp = repmat(rho,Nyrs2,Nyrs2,1);

    exp_val = repmat(abs(repmat([1:Nyrs2],Nyrs2,1)-repmat([1:Nyrs2]',1,Nyrs2)),1,1,N);
    Rhom = rho_tmp.^exp_val;
    clear rho_tmp exp_val
    Zm = mvnrnd(zeros(N,Nyrs2), Rhom, N);
    clear Rhom
    Um = normcdf(Zm,0,1);
    clear Zm
    Prod_samps(:,Nyrs+1:length(Sim_yrs)) = (repmat(0.4*Prod_tmp(Nyrs+1:end),N,1)).*logninv(Um,zeros(N,Nyrs2),0.5*ones(N,Nyrs2))+0.70*repmat(Prod_tmp(Nyrs+1:end),N,1);
    clear Um

    switch ii
        case 1 
            shortProdsamps = Prod_samps;
        case 2         
            longProdsamps = Prod_samps;   
    end
end

Nyrs = length(Sim_yrs); 
rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs,Nyrs,1);
exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm

Feedstocksamps(:,1:Nyrs) = (repmat(0.3*prod_f',N,1)).*logninv(Um,zeros(N,Nyrs),0.5*ones(N,Nyrs))+0.95*repmat(prod_f',N,1);

DEshort = betarnd(12,12,1,N);
RFshort = 1-DEshort;

m = 0.02; v = (0.50*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
DElong= lognrnd(mu,sigma,1,N); %0.02*DEu(2,:);

m = 20; v = (0.2*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFlongBankLT = lognrnd(mu,sigma,1,N); 

m = 0.02; v = (0.50*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFfeedstock = lognrnd(mu,sigma,1,N);

Years = 1955:2018;
save('CFC113Priors.mat','Years','longProdsamps','shortProdsamps','RFlongBankLT','DElong','DEshort','RFshort','Prod_reported','RFfeedstock', 'Feedstocksamps');

