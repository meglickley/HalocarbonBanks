
% First run Develop_priors.m
% Then run Master.m
% Then run PostProcess.m

clear all
close all

load('HCFC22Priors.mat')

FileName = 'HCFC22';

Years = 1944:2019; 
HCFC22post.Years = Years;

Nyrs = length(Years);
% Lifetime estimates
LTind = 16; 
LT = sparc_lifetimes(LTind, 1944:2019); 

% Observations from agage.mit.edu
load('AGAGE_published.mat'); 
MRind = 7;
tmp = MF(:,MRind);
Yrobs = MF(:,1);
tmp_ind = ~isnan(tmp); 
tmp = tmp(tmp_ind); 
tmp_yr = Yrobs(tmp_ind); 
tmp_yr = floor(tmp_yr); 
YRobs = unique(tmp_yr);
for yy = 1:length(YRobs)
    ind = find(tmp_yr == YRobs(yy)); 
    MRobs(yy,1) = mean(tmp(ind));
end

clear MF


% Simulation model
% Molecular weight of HCFC-22: 86.47g/mol
Molecular_weight = 86.47;
ppt_to_tonnes =  164.5365*Molecular_weight;
N = 10^6;
Nresamps = 10^5;
Nyrs = length(Years(1):2019);

% Initiating simulation model
MF1 = 1; %ppt
t = 1;

% short banks
indx = datasample([1:size(shortProdsamps,1)],N);

prod_s = shortProdsamps(indx,:);
de_s   = RFshortBank.y1(indx); 

e_s(:,t) = prod_s(:,t).*de_s';
b_s(:,t) = prod_s(:,t).*(1 - de_s');

% medium banks
indx = datasample([1:size(shortProdsamps,1)],N);

prod_m = mediumProdsamps(indx,:);
de_m   = RFmediumBank.y1(:,indx); 
rf_m   = RFmediumBank.LT(indx); 

[b, e] = Bank_Emiss(prod_m(:,t), rf_m', de_m(t,:)',zeros(N,1));
e_m(:,t) = e; 
b_m(:,t) = b;

% long banks
indx = datasample([1:size(shortProdsamps,1)],N);

prod_l = longProdsamps(indx,:);
rf_l   = RFlongBank.LT(indx); 

[b, e] = Bank_Emiss(prod_l (:,t), rf_l', zeros(N,1),zeros(N,1));
e_l(:,t) = e; 
b_l(:,t) = b;

% feedstock emissions
indx = datasample([1:size(shortProdsamps,1)],N);
prod_f = feedstocksamps(indx,:);
rf_f   = RFfeedstock.RF(indx); 

e_f(:,t) = rf_f'.*prod_f(:,t); 

Emiss(:,t) = 1000*(e_s(:, t) + e_m(:, t) + e_l(:, t) + e_f(:, t)); 

MF(:,t+1) = exp(-1/LT(t))*MF1+(1/ppt_to_tonnes).*Emiss(:,t);


for t = 2:Nyrs

    e_s(:,t) = prod_s(:,t).*de_s' + b_s(:, t-1);
    b_s(:,t) = prod_s(:,t).*(1 - de_s');

    [b, e] = Bank_Emiss(prod_m(:,t),rf_m', de_m(t,:), b_m(:,t-1));
    
    e_m(:,t) = e; 
    b_m(:,t) = b;

    [b, e] = Bank_Emiss(prod_l(:,t), rf_l', zeros(N,1),  b_l(:,t-1));

    e_l(:,t) = e; 
    b_l(:,t) = b;

    e_f(:,t) = rf_f'.*prod_f(:,t); 

    Emiss(:,t) = 1000*(e_s(:, t) + e_m(:, t) + e_l(:, t) + e_f(:, t));
    
    LT_t(t) = LT(find(LTYears == Years(t))); 
    MF(:,t+1) = exp(-1/(LT_t(t)))*MF(:,t)+(1/ppt_to_tonnes)*Emiss(:,t);
end


Bank = b_s + b_m + b_l;

% Likelihood Section
ytmp1 = find(Years == 2004);
ytmp2 = find(Years == 2019);
N_tmp = length(ytmp1:ytmp2);

ytmp3 = find(YRobs == 2004); 
ytmp4 = find(YRobs == 2019); 

exp_val = abs(repmat([1:N_tmp],N_tmp,1)-repmat([1:N_tmp]',1,N_tmp));

% Uncertainty priors
MF_sigma_prior = (0.03*ones(N,1)+0.03*betarnd(2,2,N,1))*MRobs(ytmp3:ytmp4)';
parfor ii = 1:N
    
    rho_tmp = 0.98*ones(N_tmp,N_tmp);
    Rho = rho_tmp.^exp_val;
    Cov_matrix = Rho.*(MF_sigma_prior(ii,:)'*MF_sigma_prior(ii,:));
    
    likelihood1(ii) = mvnpdf(MF(ii, ytmp1+1:ytmp2+1)', MRobs(ytmp3:ytmp4), Cov_matrix);
end

Obs_Emiss = ppt_to_tonnes*(MRobs(ytmp3+1:ytmp4+1) - MRobs(ytmp3:ytmp4).*exp(-1./(LT_total(end-15:end)')));

likelihood1 = (1/sum(likelihood1))*likelihood1;

% check to make sure sufficient resampling from priors
figure; 
plot((1/sum(likelihood1))*cumsum(likelihood1))
title('Important Ratio HCFC-22');

% Resample from the priors based on the relative likelihood
ResampleIndex = nan(Nresamps,1);

% Updating priors based on SIR/likelihood 
IR_vec1 = 1/sum(likelihood1)*cumsum(likelihood1);

parfor ii = 1:Nresamps
    ResampleIndex(ii) = find(IR_vec1>rand,1);
end


figure; 
plot(Years, MF(ResampleIndex(1:100),2:end), 'k'); 
hold on; 
plot(2004:2019 , MRobs(ytmp3:ytmp4), '--r', 'LineWidth',2)


HCFC22post.rf_f = rf_f(ResampleIndex);
HCFC22post.rf_m = rf_m(ResampleIndex);
HCFC22post.rf_l = rf_l(ResampleIndex);
HCFC22post.de_s = de_s(ResampleIndex);
HCFC22post.de_m = de_m(ResampleIndex);

HCFC22post.prod_s = prod_s(ResampleIndex,:);
HCFC22post.prod_m = prod_m(ResampleIndex,:);
HCFC22post.prod_l = prod_l(ResampleIndex,:);
HCFC22post.prod_f = prod_f(ResampleIndex,:);
HCFC22post.emiss_s = e_s(ResampleIndex,:);
HCFC22post.emiss_l = e_l(ResampleIndex,:);
HCFC22post.emiss_m = e_m(ResampleIndex,:);
HCFC22post.emiss_f = e_f(ResampleIndex,:);

HCFC22post.shortBank = b_s(ResampleIndex,:);
HCFC22post.medBank = b_m(ResampleIndex,:);
HCFC22post.longBank = b_l(ResampleIndex,:);
HCFC22post.ObsEmiss = Obs_Emiss;
HCFC22post.Years = Years;
HCFC22post.MF = MF(ResampleIndex,2:end);

str = strcat(ScenarioName, '.mat'); 
save(str, 'HCFC22post'); 


%% HCFC-141b
clear all
load('HCFC141bPriors.mat')

% Lifetime estimates
LTind = 17; 
LTYears = 1944:2019;
LT = sparc_lifetimes(LTind, LTYears); 

FileName = 'HCFC141b';
load('HCFC141bPriors.mat');


load('AGAGE_published.mat'); 
MRind = 8;
tmp = MF(:,MRind);
Yrobs = MF(:,1);
tmp_ind = ~isnan(tmp); 
tmp = tmp(tmp_ind); 
tmp_yr = Yrobs(tmp_ind); 
tmp_yr = floor(tmp_yr); 
YRobs = unique(tmp_yr);
for yy = 1:length(YRobs)
    ind = find(tmp_yr == YRobs(yy)); 
    MRobs(yy,1) = mean(tmp(ind));
end

clear MF

% Simulation model
% Molecular weight of HCFC-141b: 116.94g/mol
Molecular_weight = 116.94;
ppt_to_tonnes =  164.5365*Molecular_weight;
N = 10^6;
Nresamps = 10^5;
Nyrs = length(Years(1):2019);


% Initiate Simulation model
MF1 = 0; %ppt
t = 1;

% short banks
indx = datasample([1:size(shortProdsamps,1)],N);

prod_s = shortProdsamps(indx,:);
de_s   = RFshortBank.y1(indx); 

e_s(:,t) = prod_s(:,t).*de_s';
b_s(:,t) = prod_s(:,t).*(1 - de_s');

% medium banks
indx = datasample([1:size(shortProdsamps,1)],N);

prod_m = mediumProdsamps(indx,:);
de_m   = RFmediumBank.y1(:,indx); 
rf_m   = RFmediumBank.LT(indx); 

[b, e] = Bank_Emiss(prod_m(:,t), rf_m', de_m',zeros(N,1));
e_m(:,t) = e; 
b_m(:,t) = b;

% long banks
indx = datasample([1:size(shortProdsamps,1)],N);

prod_l = longProdsamps(indx,:);
de_l   = RFlongBank.y1(:,indx); 
rf_l   = RFlongBank.LT(indx); 

[b, e] = Bank_Emiss(prod_l (:,t), rf_l', de_l',zeros(N,1));
e_l(:,t) = e; 
b_l(:,t) = b;


Emiss(:,t) = 1000*(e_s(:, t) + e_m(:, t) + e_l(:, t)); 

LT_t = LT(find(LTYears == Years(t))); 

MF(:,t+1) = exp(-1/LT_t)*MF1+(1/ppt_to_tonnes).*Emiss(:,t);
for t = 2:Nyrs
    
    e_s(:,t) = prod_s(:,t).*de_s' + b_s(:, t-1);
    b_s(:,t) = prod_s(:,t).*(1 - de_s');

    [b, e] = Bank_Emiss(prod_m(:,t),rf_m', de_m', b_m(:,t-1));
    
    e_m(:,t) = e; 
    b_m(:,t) = b;

    [b, e] = Bank_Emiss(prod_l(:,t), rf_l', de_l',  b_l(:,t-1));

    e_l(:,t) = e; 
    b_l(:,t) = b;

    Emiss(:,t) = 1000*(e_s(:, t) + e_m(:, t) + e_l(:, t)); 

    LT_t(t) = LT(find(LTYears == Years(t))); 
    MF(:,t+1) = exp(-1/LT_t(t))*MF(:,t)+(1/ppt_to_tonnes)*Emiss(:,t);
end

Bank = b_s + b_m + b_l;

% Likelihood Section
ytmp1 = find(Years == 2004);
ytmp2 = find(Years == 2019);
N_tmp = length(ytmp1:ytmp2);

ytmp3 = find(YRobs == 2004); 
ytmp4 = find(YRobs == 2019); 

exp_val = abs(repmat([1:N_tmp],N_tmp,1)-repmat([1:N_tmp]',1,N_tmp));


MF_sigma_prior = (0.035*ones(N,1)+0.025*betarnd(2,2,N,1))*MRobs(ytmp3:ytmp4)';

parfor ii = 1:N
    
    rho_tmp = 0.95*ones(N_tmp,N_tmp);
    Rho = rho_tmp.^exp_val;
    
    Cov_matrix = Rho.*(MF_sigma_prior(ii,:)'*MF_sigma_prior(ii,:));
    
    likelihood1(ii) = mvnpdf(MF(ii, ytmp1+1:ytmp2+1)', MRobs(ytmp3:ytmp4), Cov_matrix);
end

likelihood1 = (1/sum(likelihood1))*likelihood1;

Obs_Emiss = ppt_to_tonnes*(MRobs(ytmp3+1:ytmp4+1) - MRobs(ytmp3:ytmp4).*exp(-1./(LT_total(end-15:end)')));



figure; 
plot((1/sum(likelihood1))*cumsum(likelihood1))
title('Important Ratio HCFC-141b');

% Resample from the priors based on the relative likelihood
ResampleIndex = nan(Nresamps,1);

% Updating priors based on SIR/likelihood 
IR_vec1 = 1/sum(likelihood1)*cumsum(likelihood1);

parfor ii = 1:Nresamps
    ResampleIndex(ii) = find(IR_vec1>rand,1);
end


HCFC141bpost.rf_m = rf_m(ResampleIndex);
HCFC141bpost.rf_l = rf_l(ResampleIndex);
HCFC141bpost.de_s = de_s(ResampleIndex);
HCFC141bpost.de_l = de_l(ResampleIndex);
HCFC141bpost.de_m = de_m(ResampleIndex);

HCFC141bpost.prod_s = prod_s(ResampleIndex,1:31);
HCFC141bpost.prod_m = prod_m(ResampleIndex,1:31);
HCFC141bpost.prod_l = prod_l(ResampleIndex,1:31);
HCFC141bpost.emiss_s = e_s(ResampleIndex,:);
HCFC141bpost.emiss_l = e_l(ResampleIndex,:);
HCFC141bpost.emiss_m = e_m(ResampleIndex,:);

HCFC141bpost.shortBank = b_s(ResampleIndex,:);
HCFC141bpost.medBank = b_m(ResampleIndex,:);
HCFC141bpost.longBank = b_l(ResampleIndex,:);
HCFC141bpost.ObsEmiss = Obs_Emiss;
HCFC141bpost.Years = Years;
HCFC141bpost.MF = MF(ResampleIndex,2:end);

str = strcat(FileName, '.mat'); 
save(str, 'HCFC141bpost'); 

%% HCFC-142b
clear all
load('HCFC142bPriors.mat')
%Years = 1981:2019;

% Lifetime estimates
LTind = 18; 
LTYears = 1944:2020;
LT = sparc_lifetimes(LTind, LTYears); 

FileName = 'HCFC142b';


load('AGAGE_published.mat'); 
MRind = 9;
tmp = MF(:,MRind);
Yrobs = MF(:,1);
tmp_ind = ~isnan(tmp); 
tmp = tmp(tmp_ind); 
tmp_yr = Yrobs(tmp_ind); 
tmp_yr = floor(tmp_yr); 
YRobs = unique(tmp_yr);
for yy = 1:length(YRobs)
    ind = find(tmp_yr == YRobs(yy)); 
    MRobs(yy,1) = mean(tmp(ind));
end

clear MF

% Molecular weight of HCFC-142b: 100.495g/mol
Molecular_weight = 100.495;
ppt_to_tonnes =  164.5365*Molecular_weight;
N = 10^6;
Nresamps = 10^5;
Nyrs = length(Years(1):2020);


MF1 = 0; %ppt

t = 1;

indx = datasample([1:size(shortProdsamps,1)],N);

prod_s = shortProdsamps(indx,:);
de_s   = RFshortBank.y1(indx); 

e_s(:,t) = prod_s(:,t).*de_s';
b_s(:,t) = prod_s(:,t).*(1 - de_s');

indx = datasample([1:size(shortProdsamps,1)],N);

prod_m = mediumProdsamps(indx,:);
de_m   = RFmediumBank.y1(:,indx); 
rf_m   = RFmediumBank.LT(indx); 

[b, e] = Bank_Emiss(prod_m(:,t), rf_m', de_m',zeros(N,1));
e_m(:,t) = e; 
b_m(:,t) = b;

indx = datasample([1:size(shortProdsamps,1)],N);

prod_l = longProdsamps(indx,:);
de_l   = RFlongBank.y1(:,indx); 
rf_l   = RFlongBank.LT(indx); 

[b, e] = Bank_Emiss(prod_l (:,t), rf_l', de_l',zeros(N,1));
e_l(:,t) = e; 
b_l(:,t) = b;

% Feedstock
indx = datasample([1:size(shortProdsamps,1)],N);

prod_f = Feedstocksamps(indx,:);
de_f   = RFfeedstock.RF(indx); 

e_f(:,t) = de_f'.*prod_f(:,t); 

Emiss(:,t) = 1000*(e_s(:, t) + e_m(:, t) + e_l(:, t) + e_f(:,t)); 

LT_t = LT(find(LTYears == Years(t))); 

MF(:,t+1) = exp(-1/LT_t)*MF1+(1/ppt_to_tonnes).*Emiss(:,t);
for t = 2:Nyrs
    e_s(:,t) = prod_s(:,t).*de_s' + b_s(:, t-1);
    b_s(:,t) = prod_s(:,t).*(1 - de_s');

    [b, e] = Bank_Emiss(prod_m(:,t),rf_m', de_m', b_m(:,t-1));
    
    e_m(:,t) = e; 
    b_m(:,t) = b;

    [b, e] = Bank_Emiss(prod_l(:,t), rf_l', de_l',  b_l(:,t-1));

    e_l(:,t) = e; 
    b_l(:,t) = b;

    e_f(:,t) = de_f'.*prod_f(:,t); 
    
    Emiss(:,t) = 1000*(e_s(:, t) + e_m(:, t) + e_l(:, t) + e_f(:, t)); 

    LT_total(t) = LT(find(LTYears == Years(t))); 
    
    MF(:,t+1) = exp(-1/LT_total(t))*MF(:,t)+(1/ppt_to_tonnes)*Emiss(:,t);
end

Bank = b_s + b_m + b_l;

% Likelihood Section
ytmp1 = find(Years == 2004);
ytmp2 = find(Years == 2019);
N_tmp = length(ytmp1:ytmp2);

ytmp3 = find(YRobs == 2004); 
ytmp4 = find(YRobs == 2019); 

exp_val = abs(repmat([1:N_tmp],N_tmp,1)-repmat([1:N_tmp]',1,N_tmp));


MF_sigma_prior = (0.028*ones(N,1)+0.022*betarnd(2,2,N,1))*MRobs(ytmp3:ytmp4)';

corrparam = 0.96;

parfor ii = 1:N
    rho_tmp = corrparam*ones(N_tmp,N_tmp);
    Rho = rho_tmp.^exp_val;
    
    Cov_matrix = Rho.*(MF_sigma_prior(ii,:)'*MF_sigma_prior(ii,:));
    
    likelihood1(ii) = mvnpdf(MF(ii, ytmp1+1:ytmp2+1)', MRobs(ytmp3:ytmp4), Cov_matrix);
end

likelihood1 = (1/sum(likelihood1))*likelihood1;


Obs_Emiss = ppt_to_tonnes*(MRobs(ytmp3+1:ytmp4+1) - MRobs(ytmp3:ytmp4).*exp(-1./(LT_total(end-15:end)')));

figure; 
plot((1/sum(likelihood1))*cumsum(likelihood1))
title('Important Ratio HCFC-142b');

% Resample from the priors based on the relative likelihood
ResampleIndex = nan(Nresamps,1);

% Updating priors based on SIR/likelihood 
IR_vec1 = 1/sum(likelihood1)*cumsum(likelihood1);

parfor ii = 1:Nresamps
    ResampleIndex(ii) = find(IR_vec1>rand,1);
end

% Model each bank separately. 

figure; 
subplot(1,2,1); 
plot(Years(1:39), MF(1:100,2:end), 'k'); 
hold on; 
plot(2004:2019 , MRobs(ytmp3:ytmp4), '--r', 'LineWidth',2)

subplot(1,2,2); 
plot(Years(1:39), MF(ResampleIndex(1:1000),2:end), 'k'); 
hold on; 
plot(2004:2019 , MRobs(ytmp3:ytmp4), '--r', 'LineWidth',2)



HCFC142bpost.de_f = de_f(ResampleIndex);
HCFC142bpost.rf_m = rf_m(ResampleIndex);
HCFC142bpost.rf_l = rf_l(ResampleIndex);
HCFC142bpost.de_s = de_s(ResampleIndex);
HCFC142bpost.de_l = de_l(ResampleIndex);
HCFC142bpost.de_m = de_m(ResampleIndex);
HCFC142bpost.prod_s = prod_s(ResampleIndex,:);
HCFC142bpost.prod_m = prod_m(ResampleIndex,:);
HCFC142bpost.prod_l = prod_l(ResampleIndex,:);
HCFC142bpost.prod_f = prod_f(ResampleIndex,:);
HCFC142bpost.emiss_s = e_s(ResampleIndex,:);
HCFC142bpost.emiss_l = e_l(ResampleIndex,:);
HCFC142bpost.emiss_m = e_m(ResampleIndex,:);

HCFC142bpost.shortBank = b_s(ResampleIndex,:);
HCFC142bpost.medBank = b_m(ResampleIndex,:);
HCFC142bpost.longBank = b_l(ResampleIndex,:);
HCFC142bpost.ObsEmiss = Obs_Emiss;
HCFC142bpost.Years = Years;
HCFC142bpost.MF = MF(ResampleIndex,2:end);

str = strcat(FileName, '.mat'); 
save(str, 'HCFC142bpost'); 

%% CFC-113
clear all
load('CFC113Priors.mat')
%Years = 1944:2019;
FileName = 'CFC113';
% Lifetime estimates
LTind = 3; 
LTYears = 1935:2020;
LT = sparc_lifetimes(LTind, LTYears); 


load('AGAGE_published.mat'); 
MRind = 4;
tmp = MF(:,MRind);
Yrobs = MF(:,1);
tmp_ind = ~isnan(tmp); 
tmp = tmp(tmp_ind); 
tmp_yr = Yrobs(tmp_ind); 
tmp_yr = floor(tmp_yr); 
YRobs = unique(tmp_yr);
for yy = 1:length(YRobs)
    ind = find(tmp_yr == YRobs(yy)); 
    MRobs(yy,1) = mean(tmp(ind));
end

clear MF


% Molecular weight of CFC-113: 187.375g/mol
Molecular_weight = 187.375;
ppt_to_tonnes =  164.5365*Molecular_weight;
N = 10^6;
Nresamps = 10^5;
Nyrs = length(Years(1):2018);


% Initiate Simulation Model
MF1 = 0; %ppt

t = 1;

% Short banks
indx = datasample([1:size(shortProdsamps,1)],N);

prod_s = shortProdsamps(indx,:);
de_s   = DEshort(indx); 

e_s(:,t) = prod_s(:,t).*de_s';
b_s(:,t) = prod_s(:,t).*(1 - de_s');

% Long banks
indx = datasample([1:size(shortProdsamps,2)],N);

prod_l = longProdsamps(indx,:);
de_l   = DElong(:,indx); 
rf_l   = 1-exp(-1./RFlongBankLT(indx)); 

% Feedstocks
indx = datasample([1:size(shortProdsamps,2)],N);
de_f = RFfeedstock(indx); 
prod_f = 1000*Feedstocksamps(indx, :); 

[b, e] = Bank_Emiss(prod_l(:,t), rf_l', de_l',zeros(N,1));
e_l(:,t) = e; 
b_l(:,t) = b;


Emiss(:,t) = (e_s(:, t) +  e_l(:, t)); 

LT_t = LT(find(LTYears == Years(t))); 

MF(:,t+1) = exp(-1/LT_t)*MF1+(1/ppt_to_tonnes).*Emiss(:,t);
for t = 2:Nyrs
    e_s(:,t) = prod_s(:,t).*de_s' + b_s(:, t-1);
    b_s(:,t) = prod_s(:,t).*(1 - de_s');

    [b, e] = Bank_Emiss(prod_l(:,t), rf_l', de_l',  b_l(:,t-1));

    e_l(:,t) = e; 
    b_l(:,t) = b;
    
    e_f(:,t) = de_f'.*prod_f(:,t); 

    Emiss(:,t) = e_s(:, t) +  e_l(:, t) + e_f(:, t); 

    LT_t = LT(find(LTYears == Years(t))); 
    MF(:,t+1) = exp(-1/LT_t)*MF(:,t)+(1/ppt_to_tonnes)*Emiss(:,t);
end

Bank = b_s +  b_l;


% Likelihood Section
ytmp1 = find(Years == 2004);
ytmp2 = find(Years == 2018);
N_tmp = length(ytmp1:ytmp2);

ytmp3 = find(YRobs == 2004); 
ytmp4 = find(YRobs == 2018); 

Obs_Emiss = ppt_to_tonnes*(MRobs(ytmp3+1:ytmp4+1) - MRobs(ytmp3:ytmp4).*exp(-1./(LT(end-15:end-1)')));

exp_val = abs(repmat([1:N_tmp],N_tmp,1)-repmat([1:N_tmp]',1,N_tmp));


MF_sigma_prior = (0.01*ones(N,1)+0.03*betarnd(2,2,N,1))*MRobs(ytmp3:ytmp4)';

corrparam = 0.95;

parfor ii = 1:N
    rho_tmp = corrparam*ones(N_tmp,N_tmp);
    Rho = rho_tmp.^exp_val;
    
    Cov_matrix = Rho.*(MF_sigma_prior(ii,:)'*MF_sigma_prior(ii,:));
    
    likelihood1(ii) = mvnpdf(MF(ii, ytmp1:ytmp2)', MRobs(ytmp3:ytmp4), Cov_matrix);
end

likelihood1 = (1/sum(likelihood1))*likelihood1;

figure; 
plot((1/sum(likelihood1))*cumsum(likelihood1))
title('Important Ratio CFC-113');


% Resample from the priors based on the relative likelihood
ResampleIndex = nan(Nresamps,1);

% Updating priors based on SIR/likelihood 
IR_vec1 = 1/sum(likelihood1)*cumsum(likelihood1);

parfor ii = 1:Nresamps
    ResampleIndex(ii) = find(IR_vec1>rand,1);
end

% Model each bank separately. 

figure; 
subplot(1,2,1); 
plot(Years(1:end), MF(1:100,2:end), 'k'); 
hold on; 
plot(2004:2018 , MRobs(ytmp3:ytmp4), '--r', 'LineWidth',2)

subplot(1,2,2); 
plot(Years(1:end), MF(ResampleIndex(1:1000),2:end), 'k'); 
hold on; 
plot(2004:2018 , MRobs(ytmp3:ytmp4), '--r', 'LineWidth',2)



CFC113post.rf_l = rf_l(ResampleIndex);
CFC113post.de_s = de_s(ResampleIndex);
CFC113post.de_l = de_l(ResampleIndex);
CFC113post.de_f = de_f(ResampleIndex);

CFC113post.prod_s = prod_s(ResampleIndex,:);
CFC113post.prod_l = prod_l(ResampleIndex,:);
CFC113post.prod_f = prod_f(ResampleIndex,:);

CFC113post.emiss_s = e_s(ResampleIndex,:);
CFC113post.emiss_l = e_l(ResampleIndex,:);
CFC113post.emiss_f = e_f(ResampleIndex,:);

CFC113post.shortBank = b_s(ResampleIndex,:);

CFC113post.longBank = b_l(ResampleIndex,:);
CFC113post.ObsEmiss = Obs_Emiss;
CFC113post.Years = Years;
CFC113post.MF = MF(ResampleIndex,2:end);


str = strcat(FileName, '.mat'); 
save(str, 'CFC113post'); 

%% CFC-114

clear all
load('CFC114Priors.mat')
%Years = 1944:2019;
FileName = 'CFC114';
% Lifetime estimates
LTind = 4; 
LTYears = 1935:2020;
LT = sparc_lifetimes(LTind, LTYears); 


load('AGAGE_published.mat'); 
MRind = 5;
tmp = MF(:,MRind);
Yrobs = MF(:,1);
tmp_ind = ~isnan(tmp); 
tmp = tmp(tmp_ind); 
tmp_yr = Yrobs(tmp_ind); 
tmp_yr = floor(tmp_yr); 
YRobs = unique(tmp_yr);
for yy = 1:length(YRobs)
    ind = find(tmp_yr == YRobs(yy)); 
    MRobs(yy,1) = mean(tmp(ind));
end

clear MF


% Molecular weight of CFC-114: 170.92g/mol
Molecular_weight = 170.92;
ppt_to_tonnes =  164.5365*Molecular_weight;
N = 10^6;
Nresamps = 10^5;
Nyrs = length(Years(1):2020);

% Initiate Simulation Model
MF1 = 0; %ppt

t = 1;

% Short banks
indx = datasample([1:size(shortProdsamps,1)],N);

prod_s = shortProdsamps(indx,:);
de_s   = RFshortBank.y1(indx); 

e_s(:,t) = prod_s(:,t).*de_s';
b_s(:,t) = prod_s(:,t).*(1 - de_s');

% Long banks
indx = datasample([1:size(shortProdsamps,2)],N);

prod_l = longProdsamps(indx,:);
de_l   = RFlongBank.y1(:,indx); 
rf_l   = RFlongBank.LT(indx); 

[b, e] = Bank_Emiss(prod_l (:,t), rf_l', de_l',zeros(N,1));
e_l(:,t) = e; 
b_l(:,t) = b;


Emiss(:,t) = 1000*(e_s(:, t) +  e_l(:, t)); 

LT_t = LT(find(LTYears == Years(t))); 

MF(:,t+1) = exp(-1/LT_t)*MF1+(1/ppt_to_tonnes).*Emiss(:,t);
for t = 2:Nyrs
    e_s(:,t) = prod_s(:,t).*de_s' + b_s(:, t-1);
    b_s(:,t) = prod_s(:,t).*(1 - de_s');

    [b, e] = Bank_Emiss(prod_l(:,t), rf_l', de_l',  b_l(:,t-1));

    e_l(:,t) = e; 
    b_l(:,t) = b;

    Emiss(:,t) = 1000*(e_s(:, t) +  e_l(:, t)); 

    LT_t = LT(find(LTYears == Years(t))); 
    MF(:,t+1) = exp(-1/LT_t)*MF(:,t)+(1/ppt_to_tonnes)*Emiss(:,t);
end

Bank = b_s +  b_l;


% Likelihood Section
ytmp1 = find(Years == 2004);
ytmp2 = find(Years == 2019);
N_tmp = length(ytmp1:ytmp2);

ytmp3 = find(YRobs == 2004); 
ytmp4 = find(YRobs == 2019); 

Obs_Emiss = ppt_to_tonnes*(MRobs(ytmp3+1:ytmp4+1) - MRobs(ytmp3:ytmp4).*exp(-1./(LT(end-15:end)')));

exp_val = abs(repmat([1:N_tmp],N_tmp,1)-repmat([1:N_tmp]',1,N_tmp));

MF_sigma_prior = (0.01*ones(N,1)+0.02*betarnd(2,2,N,1))*MRobs(ytmp3:ytmp4)';

corrparam = 0.7;

parfor ii = 1:N
    
    rho_tmp = corrparam*ones(N_tmp,N_tmp);
    Rho = rho_tmp.^exp_val;
    Cov_matrix = Rho.*(MF_sigma_prior(ii,:)'*MF_sigma_prior(ii,:));
    
    likelihood1(ii) = mvnpdf(MF(ii, ytmp1+1:ytmp2+1)', MRobs(ytmp3:ytmp4), Cov_matrix);
end

likelihood1 = (1/sum(likelihood1))*likelihood1;

figure; 
plot((1/sum(likelihood1))*cumsum(likelihood1))
title('Important Ratio CFC-114');


% Resample from the priors based on the relative likelihood
ResampleIndex = nan(Nresamps,1);

% Updating priors based on SIR/likelihood 
IR_vec1 = 1/sum(likelihood1)*cumsum(likelihood1);

parfor ii = 1:Nresamps
    ResampleIndex(ii) = find(IR_vec1>rand,1);
end

% Model each bank separately. 

figure; 
subplot(1,2,1); 
plot(Years(11:end), MF(1:100,2:end), 'k'); 
hold on; 
plot(2004:2019 , MRobs(ytmp3:ytmp4), '--r', 'LineWidth',2)

subplot(1,2,2); 
plot(Years(11:end), MF(ResampleIndex(1:1000),2:end), 'k'); 
hold on; 
plot(2004:2019 , MRobs(ytmp3:ytmp4), '--r', 'LineWidth',2)


CFC114post.rf_l = rf_l(ResampleIndex);
CFC114post.de_s = de_s(ResampleIndex);
CFC114post.de_l = de_l(ResampleIndex);

CFC114post.prod_s = prod_s(ResampleIndex,:);

CFC114post.prod_l = prod_l(ResampleIndex,:);

CFC114post.emiss_s = e_s(ResampleIndex,:);
CFC114post.emiss_l = e_l(ResampleIndex,:);


CFC114post.shortBank = b_s(ResampleIndex,:);

CFC114post.longBank = b_l(ResampleIndex,:);
CFC114post.ObsEmiss = Obs_Emiss;
CFC114post.Years = Years;
CFC114post.MF = MF(ResampleIndex,2:end);


str = strcat(FileName, '.mat'); 
save(str, 'CFC114post'); 

%% CFC-115
% Lifetime of CFC-115 is so long that we assume it doesn't decay in this
% time period. 
clear all
load('CFC115Priors.mat')
%Years = 1944:2019;
ScenarioName = 'CFC115_original';
% Lifetime estimates
LTind = 5; 
LTYears = 1930:2019;
LT = sparc_lifetimes(LTind, LTYears); 


load('AGAGE_published.mat'); 
MRind = 6;
tmp = MF(:,MRind);
Yrobs = MF(:,1);
tmp_ind = ~isnan(tmp); 
tmp = tmp(tmp_ind); 
tmp_yr = Yrobs(tmp_ind); 
tmp_yr = floor(tmp_yr); 
YRobs = unique(tmp_yr);
for yy = 1:length(YRobs)
    ind = find(tmp_yr == YRobs(yy)); 
    MRobs(yy,1) = mean(tmp(ind));
end

clear MF


% Molecular weight of CFC-114: 153.961g/mol
Molecular_weight = 153.961;
ppt_to_tonnes =  164.5365*Molecular_weight;
N = 10^6;
Nresamps = 10^5;
Nyrs = length(Years(1):2020);

% Initiate Simulation Model
MF1 = 0; %ppt

t = 1;

% Short banks
indx = datasample([1:size(shortProdsamps,1)],N);

prod_s = shortProdsamps(indx,:);
prod_s(isnan(prod_s)) = 0;

de_s   = RFshortBank.y1(indx); 

e_s(:,t) = prod_s(:,t).*de_s';
b_s(:,t) = prod_s(:,t).*(1 - de_s');

% Long banks
indx = datasample([1:size(shortProdsamps,2)],N);

prod_l = longProdsamps(indx,:);
prod_l(isnan(prod_l)) = 0;

de_l   = RFlongBank.y1(:,indx); 
rf_l   = RFlongBank.LT(indx); 


[b, e] = Bank_Emiss(prod_l (:,t), rf_l', de_l',zeros(N,1));
e_l(:,t) = e; 
b_l(:,t) = b;

Emiss(:,t) = 1000*(e_s(:, t) +  e_l(:, t)); 


LT_t = LT(find(LTYears == Years(t))); 
if isempty(LT_t)
    MF(:,t+1) = MF1+(1/ppt_to_tonnes).*Emiss(:,t);
else
    MF(:,t+1) = exp(-1/LT_t)*MF1+(1/ppt_to_tonnes).*Emiss(:,t);
end

for t = 2:Nyrs
    e_s(:,t) = prod_s(:,t).*de_s' + b_s(:, t-1);
    b_s(:,t) = prod_s(:,t).*(1 - de_s');

    [b, e] = Bank_Emiss(prod_l(:,t), rf_l', de_l',  b_l(:,t-1));

    e_l(:,t) = e; 
    b_l(:,t) = b;

    Emiss(:,t) = 1000*(e_s(:, t) +  e_l(:, t)); 

    LT_t = LT(find(LTYears == Years(t))); 
    
    if isempty(LT_t)
        MF(:,t+1) = MF(:,t)+(1/ppt_to_tonnes).*Emiss(:,t);
    else
        MF(:,t+1) = exp(-1/LT_t)*MF(:,t)+(1/ppt_to_tonnes)*Emiss(:,t);
    end
end

Bank = b_s +  b_l;

% Likelihood Section
ytmp1 = find(Years == 2004);
ytmp2 = find(Years == 2019);
N_tmp = length(ytmp1:ytmp2);

ytmp3 = find(YRobs == 2004); 
ytmp4 = find(YRobs == 2019); 

exp_val = abs(repmat([1:N_tmp],N_tmp,1)-repmat([1:N_tmp]',1,N_tmp));
Obs_Emiss = ppt_to_tonnes*(MRobs(ytmp3+1:ytmp4+1) - MRobs(ytmp3:ytmp4).*exp(-1./(LT(end-15:end)')));


MF_sigma_prior = (0.05*ones(N,1)+0.03*betarnd(2,2,N,1))*MRobs(ytmp3:ytmp4)';

corrparam = 0.3; % Low correlation of errors

parfor ii = 1:N
    rho_tmp = corrparam*ones(N_tmp,N_tmp);
    Rho = rho_tmp.^exp_val;
    
    Cov_matrix = Rho.*(MF_sigma_prior(ii,:)'*MF_sigma_prior(ii,:));
    
    likelihood1(ii) = mvnpdf(MF(ii, ytmp1+1:ytmp2+1)', MRobs(ytmp3:ytmp4), Cov_matrix);
end

likelihood1 = (1/sum(likelihood1))*likelihood1;

figure; 
plot((1/sum(likelihood1))*cumsum(likelihood1))
title('Important Ratio CFC-115');


% Resample from the priors based on the relative likelihood
ResampleIndex = nan(Nresamps,1);

% Updating priors based on SIR/likelihood 
IR_vec1 = 1/sum(likelihood1)*cumsum(likelihood1);

parfor ii = 1:Nresamps
    ResampleIndex(ii) = find(IR_vec1>rand,1);
end

% Model each bank separately. 

figure; 
subplot(1,2,1); 
plot(Years(2:end-9), MF(1:100,2:end), 'k'); 
hold on; 
plot(2004:2019 , MRobs(ytmp3:ytmp4), '--r', 'LineWidth',2)

subplot(1,2,2); 
plot(Years(2:end-9), MF(ResampleIndex(1:1000),2:end), 'k'); 
hold on; 
plot(2004:2019 , MRobs(ytmp3:ytmp4), '--r', 'LineWidth',2)



CFC115post.rf_l = rf_l(ResampleIndex);
CFC115post.de_s = de_s(ResampleIndex);
CFC115post.de_l = de_l(ResampleIndex);

CFC115post.prod_s = prod_s(ResampleIndex,:);
CFC115post.prod_l = prod_l(ResampleIndex,:);

CFC115post.emiss_s = e_s(ResampleIndex,:);
CFC115post.emiss_l = e_l(ResampleIndex,:);


CFC115post.shortBank = b_s(ResampleIndex,:);
CFC115post.longBank = b_l(ResampleIndex,:);
CFC115post.ObsEmiss = Obs_Emiss;
CFC115post.Years = Years;
CFC115post.MF = MF(ResampleIndex,2:end);


str = strcat(FileName, '.mat'); 
save(str, 'CFC115post'); 
%% halon 1211


clear all
load('halonPriors.mat')

FileName = 'Halon1211';
% Lifetime estimates
LTind = 9; % halon-1301 is 10
LTYears = 1930:2020;
LT = sparc_lifetimes(LTind, LTYears); 


load('AGAGE_published.mat'); 
MRind = 10;
tmp = MF(:,MRind);
Yrobs = MF(:,1);
tmp_ind = ~isnan(tmp); 
tmp = tmp(tmp_ind); 
tmp_yr = Yrobs(tmp_ind); 
tmp_yr = floor(tmp_yr); 
YRobs = unique(tmp_yr);
for yy = 1:length(YRobs)
    ind = find(tmp_yr == YRobs(yy)); 
    MRobs(yy,1) = mean(tmp(ind));
end

clear MF


% Molecular weight of halon-1211: 165.36g/mol; 1301 is 148.9 g/mol
Molecular_weight = 165.36;
ppt_to_tonnes =  164.5365*Molecular_weight;
N = 10^6;
Nresamps = 10^5;
Nyrs = length(Years(1):2020);

% Initiate Simulation Model
MF1 = 0; %ppt
t = 1;

indx = datasample([1:size(Prodsamps1211,1)],N);

prod = Prodsamps1211(indx,:);
de   = RF1211.DEpre1990(indx); 
rf   = RF1211.RFpre1990(indx); 

e = prod(:,t).*de';
b(:,t) = prod(:,t).*(1 - de');

Bank(:,t) = b;
Emiss(:,t) = e; 


LT_t = LT(find(LTYears == Years(t))); 

if isempty(LT_t)
    MF(:,t+1) = MF1+(1/ppt_to_tonnes).*Emiss(:,t);
else
    MF(:,t+1) = exp(-1/LT_t)*MF1+(1/ppt_to_tonnes).*Emiss(:,t);
end

Nyrs1 = find(Years == 1990);
for t = 2:Nyrs1
    
    [b, e] = Bank_Emiss(prod(:,t), rf', de',  Bank(:,t-1));
    Bank(:,t) = b;
    Emiss(:,t) = e; 
    LT_t = LT(find(LTYears == Years(t))); 
    
    MF(:,t+1) = exp(-1/LT_t)*MF(:,t)+(1/ppt_to_tonnes)*Emiss(:,t);
end

de   = RF1211.DEpost1990(indx); 
rf   = RF1211.RFpost1990(indx); 
for t = Nyrs1+1:Nyrs
    
    [b, e] = Bank_Emiss(prod(:,t), rf', de',  Bank(:,t-1));
    
    Bank(:,t) = b;
    Emiss(:,t) = e; 
    
    LT_t = LT(find(LTYears == Years(t))); 
    
    MF(:,t+1) = exp(-1/LT_t)*MF(:,t)+(1/ppt_to_tonnes)*Emiss(:,t);
end

% Likelihood Section
ytmp1 = find(Years == 2004);
ytmp2 = find(Years == 2019);
N_tmp = length(ytmp1:ytmp2);

ytmp3 = find(YRobs == 2004); 
ytmp4 = find(YRobs == 2019); 

exp_val = abs(repmat([1:N_tmp],N_tmp,1)-repmat([1:N_tmp]',1,N_tmp));

Obs_Emiss = ppt_to_tonnes*(MRobs(ytmp3+1:ytmp4+1) - MRobs(ytmp3:ytmp4).*exp(-1./(LT(end-15:end)')));

MF_sigma_prior = (0.04*ones(N,1)+0.02*betarnd(2,2,N,1))*MRobs(ytmp3:ytmp4)';

corrparam = 0.7;

parfor ii = 1:N
    rho_tmp = corrparam*ones(N_tmp,N_tmp);
    Rho = rho_tmp.^exp_val;
    
    Cov_matrix = Rho.*(MF_sigma_prior(ii,:)'*MF_sigma_prior(ii,:));
    
    likelihood1(ii) = mvnpdf(MF(ii, ytmp1+1:ytmp2+1)', MRobs(ytmp3:ytmp4), Cov_matrix);
end

likelihood1 = (1/sum(likelihood1))*likelihood1;

figure; 
plot((1/sum(likelihood1))*cumsum(likelihood1))
title('Important Ratio halon-1211');

% Resample from the priors based on the relative likelihood
ResampleIndex = nan(Nresamps,1);

% Updating priors based on SIR/likelihood 
IR_vec1 = 1/sum(likelihood1)*cumsum(likelihood1);

parfor ii = 1:Nresamps
    ResampleIndex(ii) = find(IR_vec1>rand,1);
end



figure; 
subplot(1,2,1); 
plot(Years, MF(1:1000,1:end-1), 'k'); 
hold on; 
plot(2004:2019 , MRobs(ytmp3:ytmp4), '--r', 'LineWidth',2)

subplot(1,2,2); 
plot(Years, MF(ResampleIndex(1:1000),1:end-1), 'k'); 
hold on; 
plot(2004:2019 , MRobs(ytmp3:ytmp4), '--r', 'LineWidth',2)



halon1211post.rf = rf(ResampleIndex);
halon1211post.DE = de(ResampleIndex);
halon1211post.Prod = prod(ResampleIndex,:);

halon1211post.Emiss = Emiss(ResampleIndex,:);
halon1211post.Bank = Bank(ResampleIndex,:);
halon1211post.ObsEmiss = Obs_Emiss;
halon1211post.Years = Years;
halon1211post.MF = MF(ResampleIndex,2:end);


str = strcat(FileName, '.mat'); 
save(str, 'halon1211post'); 
%% Halon 1301

clear all
load('halonPriors.mat')
FileName = 'Halon1301';
%Years = 1944:2019;

% Lifetime estimates
LTind = 10; % halon-1301 is 10
LTYears = 1930:2020;
LT = sparc_lifetimes(LTind, LTYears); 


load('AGAGE_published.mat'); 
MRind = 11;
tmp = MF(:,MRind);
Yrobs = MF(:,1);
tmp_ind = ~isnan(tmp); 
tmp = tmp(tmp_ind); 
tmp_yr = Yrobs(tmp_ind); 
tmp_yr = floor(tmp_yr); 
YRobs = unique(tmp_yr);
for yy = 1:length(YRobs)
    ind = find(tmp_yr == YRobs(yy)); 
    MRobs(yy,1) = mean(tmp(ind));
end

clear MF

% Molecular weight of halon-1211: 165.36g/mol; 1301 is 148.9 g/mol
Molecular_weight = 148.9;
ppt_to_tonnes =  164.5365*Molecular_weight;
N = 10^6;
Nresamps = 10^5;
Nyrs = length(Years(1):2020);

% Initiate Simulation Model
MF1 = 0; %ppt
t = 1;

indx = datasample([1:size(Prodsamps1301,1)],N);

prod = Prodsamps1301(indx,:);
de   = RF1301.DEpre1990(indx); 
rf   = RF1301.RFpre1990(indx); 

e = prod(:,t).*de';
b(:,t) = prod(:,t).*(1 - de');

Bank(:,t) = b;
Emiss(:,t) = e; 

LT_t = LT(find(LTYears == Years(t))); 
if isempty(LT_t)
    MF(:,t+1) = MF1+(1/ppt_to_tonnes).*Emiss(:,t);
else
    MF(:,t+1) = exp(-1/LT_t)*MF1+(1/ppt_to_tonnes).*Emiss(:,t);
end

Nyrs1 = find(Years == 1990);

for t = 2:Nyrs1
    
    [b, e] = Bank_Emiss(prod(:,t), rf', de',  Bank(:,t-1));
    Bank(:,t) = b;
    Emiss(:,t) = e; 
    LT_t = LT(find(LTYears == Years(t))); 
    
    MF(:,t+1) = exp(-1/LT_t)*MF(:,t)+(1/ppt_to_tonnes)*Emiss(:,t);
end

de   = RF1301.DEpost1990(indx); 
rf   = RF1301.RFpost1990(indx); 
for t = Nyrs1+1:Nyrs
    
    [b, e] = Bank_Emiss(prod(:,t), rf', de',  Bank(:,t-1));
    
    Bank(:,t) = b;
    Emiss(:,t) = e; 
    
    LT_t = LT(find(LTYears == Years(t))); 
    
    MF(:,t+1) = exp(-1/LT_t)*MF(:,t)+(1/ppt_to_tonnes)*Emiss(:,t);
end

% Likelihood Section
ytmp1 = find(Years == 2004);
ytmp2 = find(Years == 2019);
N_tmp = length(ytmp1:ytmp2);

ytmp3 = find(YRobs == 2004); 
ytmp4 = find(YRobs == 2019); 

exp_val = abs(repmat([1:N_tmp],N_tmp,1)-repmat([1:N_tmp]',1,N_tmp));

Obs_Emiss = ppt_to_tonnes*(MRobs(ytmp3+1:ytmp4+1) - MRobs(ytmp3:ytmp4).*exp(-1./(LT(end-15:end)')));


MF_sigma_prior = (0.03*ones(N,1)+0.03*betarnd(2,2,N,1))*MRobs(ytmp3:ytmp4)';
corrparam = 0.7;

parfor ii = 1:N
    rho_tmp = corrparam*ones(N_tmp,N_tmp);
    Rho = rho_tmp.^exp_val;
    
    Cov_matrix = Rho.*(MF_sigma_prior(ii,:)'*MF_sigma_prior(ii,:));
    
    likelihood1(ii) = mvnpdf(MF(ii, ytmp1+1:ytmp2+1)', MRobs(ytmp3:ytmp4), Cov_matrix);
end

likelihood1 = (1/sum(likelihood1))*likelihood1;

figure; 
plot((1/sum(likelihood1))*cumsum(likelihood1))
title('Important Ratio halon-1301');


% Resample from the priors based on the relative likelihood
ResampleIndex = nan(Nresamps,1);

% Updating priors based on SIR/likelihood 
IR_vec1 = 1/sum(likelihood1)*cumsum(likelihood1);

parfor ii = 1:Nresamps
    ResampleIndex(ii) = find(IR_vec1>rand,1);
end


figure; 
subplot(1,2,1); 
plot(Years, MF(1:1000,1:end-1), 'k'); 
hold on; 
plot(2004:2019 , MRobs(ytmp3:ytmp4), '--r', 'LineWidth',2)

subplot(1,2,2); 
plot(Years, MF(ResampleIndex(1:1000),1:end-1), 'k'); 
hold on; 
plot(2004:2019 , MRobs(ytmp3:ytmp4), '--r', 'LineWidth',2)


halon1301post.rf = rf(ResampleIndex);
halon1301post.DE = de(ResampleIndex);
halon1301post.Prod = prod(ResampleIndex,:);

halon1301post.Emiss = Emiss(ResampleIndex,:);
halon1301post.Bank = Bank(ResampleIndex,:);
halon1301post.ObsEmiss = Obs_Emiss;
halon1301post.Years = Years;
halon1301post.MF = MF(ResampleIndex,2:end);


str = strcat(FileName, '.mat'); 
save(str, 'halon1301post'); 
%% CFC-11 by bank type

clear all

load('CFC11Priors.mat')
FileName= 'CFC11';

Years = 1955:2019; 
CFC11post.Years = Years;


Nyrs = length(Years);

% Lifetime estimates
LTind = 1; 
LT = sparc_lifetimes(LTind, 1955:2019); 

% Observations
load('MixingRatios.mat'); 
MRind = 2; 
MRobs = MixingRatios(:,MRind);
YRobs = MixingRatios(:,1);

ppt_to_tonnes =  22602.38457;
N = 10^6;
Nresamps = 10^5;
Nyrs = length(Years(1):2018);

% Initiate Simulation Model
MF1 = 3.29; %ppt
t = 1;

% short banks
indx = datasample([1:size(shortProdsamps,1)],N);

prod_s = shortProdsamps(indx,:);
de_s   = DEshort(indx); 

e_s(:,t) = prod_s(:,t).*de_s';
b_s(:,t) = prod_s(:,t).*(1 - de_s');

% medium banks
indx = datasample([1:size(shortProdsamps,1)],N);

prod_m = mediumProdsamps(indx,:);
de_m   = DEmedium(:,indx); 
rf_m   = 1-exp(-1./RFmediumLT(indx)); 

[b, e] = Bank_Emiss(prod_m(:,t), rf_m', de_m(t,:)',zeros(N,1));
e_m(:,t) = e; 
b_m(:,t) = b;

% long banks
indx = datasample([1:size(shortProdsamps,1)],N);

prod_l = longProdsamps(indx,:);
rf_l   = RFlongBank(indx); 

[b, e] = Bank_Emiss(prod_l (:,t), rf_l', zeros(N,1),5893.9*ones(N,1));
e_l(:,t) = e; 
b_l(:,t) = b;

indx = datasample([1:size(shortProdsamps,1)],N);

Emiss(:,t) = (e_s(:, t) + e_m(:, t) + e_l(:, t)); 

MF(:,t+1) = exp(-1/LT(t))*MF1+(1/ppt_to_tonnes).*Emiss(:,t);


for t = 2:Nyrs
    
    e_s(:,t) = prod_s(:,t).*de_s' + b_s(:, t-1);
    b_s(:,t) = prod_s(:,t).*(1 - de_s');

    [b, e] = Bank_Emiss(prod_m(:,t), rf_m', de_m, b_m(:,t-1));
    
    e_m(:,t) = e; 
    b_m(:,t) = b;

    [b, e] = Bank_Emiss(prod_l(:,t), rf_l', zeros(N,1),  b_l(:,t-1));

    e_l(:,t) = e; 
    b_l(:,t) = b;


    Emiss(:,t) = (e_s(:, t) + e_m(:, t) + e_l(:, t)); 

    MF(:,t+1) = exp(-1/(LT(t)))*MF(:,t)+(1/ppt_to_tonnes)*Emiss(:,t);
end


Bank = b_s + b_m + b_l;

% Likelihood Section
ytmp1 = find(Years == 1980);
ytmp2 = find(Years == 2018);
N_tmp = length(ytmp1:ytmp2);

ytmp3 = find(YRobs == 1980); 
ytmp4 = find(YRobs == 2018); 

exp_val = abs(repmat([1:N_tmp],N_tmp,1)-repmat([1:N_tmp]',1,N_tmp));


MF_sigma_prior = (0.02*ones(N,1)+0.02*betarnd(2,2,N,1))*MRobs(ytmp3:ytmp4)';


parfor ii = 1:N
    rho_tmp = 0.96*ones(N_tmp,N_tmp);
    Rho = rho_tmp.^exp_val;
    
    Cov_matrix = Rho.*(MF_sigma_prior(ii,:)'*MF_sigma_prior(ii,:));
    
    likelihood1(ii) = mvnpdf(MF(ii, ytmp1+1:ytmp2+1)', MRobs(ytmp3:ytmp4), Cov_matrix);
end

Obs_Emiss = ppt_to_tonnes*(MRobs(ytmp3:ytmp4) - MRobs(ytmp3-1:ytmp4-1).*exp(-1./(LT(end-38:end)')));

likelihood1 = (1/sum(likelihood1))*likelihood1;

figure; 
plot((1/sum(likelihood1))*cumsum(likelihood1))
title('Important Ratio CFC-11');


% Resample from the priors based on the relative likelihood
ResampleIndex = nan(Nresamps,1);

% Updating priors based on SIR/likelihood 
IR_vec1 = 1/sum(likelihood1)*cumsum(likelihood1);

parfor ii = 1:Nresamps
    ResampleIndex(ii) = find(IR_vec1>rand,1);
end


figure; 
plot(Years(2:end), MF(ResampleIndex(1:100),2:end), 'k'); 
hold on; 
plot(1980:2018 , MRobs(ytmp3:ytmp4), '--r', 'LineWidth',2)



CFC11post.rf_m = rf_m(ResampleIndex);
CFC11post.rf_l = rf_l(ResampleIndex);
CFC11post.de_s = de_s(ResampleIndex);
CFC11post.de_m = de_m(ResampleIndex);

CFC11post.prod_s = prod_s(ResampleIndex,:);
CFC11post.prod_m = prod_m(ResampleIndex,:);
CFC11post.prod_l = prod_l(ResampleIndex,:);

CFC11post.emiss_s = e_s(ResampleIndex,:);
CFC11post.emiss_l = e_l(ResampleIndex,:);
CFC11post.emiss_m = e_m(ResampleIndex,:);


CFC11post.shortBank = b_s(ResampleIndex,:);
CFC11post.medBank = b_m(ResampleIndex,:);
CFC11post.longBank = b_l(ResampleIndex,:);
CFC11post.ObsEmiss = Obs_Emiss;
CFC11post.Years = Years;
CFC11post.MF = MF(ResampleIndex,2:end);

str = strcat(ScenarioName, '.mat'); 
save(str, 'CFC11post'); 

%% CFC-12 by bank type

clear all

load('CFC12Priors.mat')

FileName= 'CFC12';

Years = 1955:2019; 
CFC12post.Years = Years;


Nyrs = length(Years);
% Lifetime estimates
LTind = 2; 
LT = sparc_lifetimes(LTind, 1955:2019); 

% Observations
load('MixingRatios.mat'); 
MRind = 3; 
MRobs = MixingRatios(:,MRind);
YRobs = MixingRatios(:,1);


ppt_to_tonnes =  19895.36;
N = 10^6;
Nresamps = 10^5;
Nyrs = length(Years(1):2018);

% Initiate Simulation Model
MF1 = 14.29; %ppt

t = 1;

% Short bank, type 1
indx = datasample([1:size(short1Prodsamps,1)],N);

prod_s1 = short1Prodsamps(indx,:);
de_s1   = DEshort1(indx); 

e_s1(:,t) = prod_s1(:,t).*de_s1';
b_s1(:,t) = prod_s1(:,t).*(1 - de_s1');

% Short bank, type 2
indx = datasample([1:size(short2Prodsamps,1)],N);

prod_s2 = short2Prodsamps(indx,:);
de_s2   = DEshort2(indx); 

e_s2(:,t) = prod_s2(:,t).*de_s2';
b_s2(:,t) = prod_s2(:,t).*(1 - de_s2');

% Medium bank
indx = datasample([1:size(short1Prodsamps,1)],N);

prod_m = mediumProdsamps(indx,:);
de_m   = DEmedium(:,indx); 
rf_m   = 1-exp(-1./RFmediumLT(indx)); 

[b, e] = Bank_Emiss(prod_m(:,t), rf_m', de_m',zeros(N,1));
e_m(:,t) = e; 
b_m(:,t) = b;

% Long bank
indx = datasample([1:size(short1Prodsamps,1)],N);

de_l   = DElong(:,indx); 
prod_l = longProdsamps(indx,:);
rf_l   = 1-exp(-1./RFlongLT(indx)); 

[b, e] = Bank_Emiss(prod_l (:,t), rf_l', de_l', 65198*ones(N,1));
e_l(:,t) = e; 
b_l(:,t) = b;

indx = datasample([1:size(short1Prodsamps,1)],N);

Emiss(:,t) = (e_s1(:, t) + e_s2(:, t) + e_m(:, t) + e_l(:, t)); 

MF(:,t+1) = exp(-1/LT(t))*MF1+(1/ppt_to_tonnes).*Emiss(:,t);


for t = 2:Nyrs
    
    e_s1(:,t) = prod_s1(:,t).*de_s1' + b_s1(:, t-1);
    b_s1(:,t) = prod_s1(:,t).*(1 - de_s1');
    
    e_s2(:,t) = prod_s2(:,t).*de_s2' + b_s2(:, t-1);
    b_s2(:,t) = prod_s2(:,t).*(1 - de_s2');


    [b, e] = Bank_Emiss(prod_m(:,t), rf_m', de_m', b_m(:,t-1));
    
    e_m(:,t) = e; 
    b_m(:,t) = b;

    [b, e] = Bank_Emiss(prod_l(:,t), rf_l', de_l',  b_l(:,t-1));

    e_l(:,t) = e; 
    b_l(:,t) = b;


    Emiss(:,t) = (e_s1(:, t) + e_s2(:, t) + e_m(:, t) + e_l(:, t)); 

    MF(:,t+1) = exp(-1/LT(t))*MF(:,t)+(1/ppt_to_tonnes)*Emiss(:,t);
end


Bank = b_s1 + b_s2 + b_m + b_l;

% Likelihood Section
ytmp1 = find(Years == 1980);
ytmp2 = find(Years == 2018);
N_tmp = length(ytmp1:ytmp2);

ytmp3 = find(YRobs == 1980); 
ytmp4 = find(YRobs == 2018); 

exp_val = abs(repmat([1:N_tmp],N_tmp,1)-repmat([1:N_tmp]',1,N_tmp));

MF_sigma_prior = (0.015*ones(N,1)+0.015*betarnd(2,2,N,1))*MRobs(ytmp3:ytmp4)';

parfor ii = 1:N
    rho_tmp = 0.98*ones(N_tmp,N_tmp);
    Rho = rho_tmp.^exp_val;
    
    Cov_matrix = Rho.*(MF_sigma_prior(ii,:)'*MF_sigma_prior(ii,:));
    
    likelihood1(ii) = mvnpdf(MF(ii, ytmp1+1:ytmp2+1)', MRobs(ytmp3:ytmp4), Cov_matrix);
end

Obs_Emiss = ppt_to_tonnes*(MRobs(ytmp3:ytmp4) - MRobs(ytmp3-1:ytmp4-1).*exp(-1./(LT(end-38:end)')));

likelihood1 = (1/sum(likelihood1))*likelihood1;

figure; 
plot((1/sum(likelihood1))*cumsum(likelihood1))
title('Important Ratio CFC-12');

% Resample from the priors based on the relative likelihood
ResampleIndex = nan(Nresamps,1);

% Updating priors based on SIR/likelihood 
IR_vec1 = 1/sum(likelihood1)*cumsum(likelihood1);

parfor ii = 1:Nresamps
    ResampleIndex(ii) = find(IR_vec1>rand,1);
end


CFC12post.rf_m = rf_m(ResampleIndex);
CFC12post.rf_l = rf_l(ResampleIndex);
CFC12post.de_s1 = de_s1(ResampleIndex);
CFC12post.de_s2 = de_s2(ResampleIndex);
CFC12post.de_m = de_m(ResampleIndex);

CFC12post.prod_s1 = prod_s1(ResampleIndex,:);
CFC12post.prod_s2 = prod_s2(ResampleIndex,:);
CFC12post.prod_m = prod_m(ResampleIndex,:);
CFC12post.prod_l = prod_l(ResampleIndex,:);

CFC12post.emiss_s1 = e_s1(ResampleIndex,:);
CFC12post.emiss_s2 = e_s2(ResampleIndex,:);
CFC12post.emiss_l = e_l(ResampleIndex,:);
CFC12post.emiss_m = e_m(ResampleIndex,:);

CFC12post.shortBank1 = b_s1(ResampleIndex,:);
CFC12post.shortBank2 = b_s2(ResampleIndex,:);
CFC12post.medBank = b_m(ResampleIndex,:);
CFC12post.longBank = b_l(ResampleIndex,:);
CFC12post.ObsEmiss = Obs_Emiss;
CFC12post.Years = Years;
CFC12post.MF = MF(ResampleIndex,2:end);


str = strcat(FileName, '.mat'); 
save(str, 'CFC12post'); 

