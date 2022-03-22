%% The data and codes for:
% Title:  Toward the general generative model of chromatographic retention 
% Authors: Agnieszka Kamedulska, £ukasz Kubik, Julia Jacyna, Wiktoria Struck-Lewicka, Micha³ Markuszewski, Pawe³ Wiczling
% Adress: Department of Biopharmaceutics and Pharmacodynamics, Medical University of Gdañsk, Gen. J. Hallera 107, 80-416 Gdañsk, Poland
% Data: 22.03.2022
%% Load data
clear all
data = readtable('Data\1-X_Bridge_Shield_C18_5cm.csv');
data.Mod  = categorical(data.Mod);
[~,~,j]=unique(data.Mod);  data.Mod2=2-j; % MeOH = 1, ACN = 2
dataNames = readtable('Data\4-compounds-names.csv');
dataACD = readtable('Data\2-ACD-pKas-logP.csv');
cov.pKaslit = dataACD{:,3:7};             % pKa values as predicted by ACD
cov.pKasliterror = dataACD{:,25:29};      % pKa error as predicted by ACD
cov.chargesA = abs(dataACD{:,13:18});     % number of ionized groups (anions)
cov.chargesB = abs(dataACD{:,19:24});     % number of ionized groups (cations)
cov.charges = cov.chargesA+cov.chargesB;  % absolute charge
cov.groupsA = diff(cov.chargesA,1,2);     % acidic group
cov.groupsB = -diff(cov.chargesB,1,2);    % basic group
cov.R = sum(cov.pKaslit<14,2);            % number of dissociation steps
cov.logP = dataACD.logP;                  % logP

%% Load checkmol data
functional_groups = readtable('Data\6-checkmol-functional-groups.csv');
functional_groups_names = readtable('Data\Legend-checkmol-functional-group-names.csv');

functional_groups=functional_groups(:,2:end);

% combine nr of caroboxylic acid and carboxyalic acid salt functional groups
functional_groups{:,76}=functional_groups{:,76}+functional_groups{:,77};       
functional_groups{functional_groups{:,202}>5.5,202} = 6; % heterocyclic compounds with more than 6 heterocycles are treated as if they have six

% exclude functional groups that repeat itself (some groups are nested)
idx_excluded = [1 2 3 6 27 28 37 47 48 51 55 61 62 67 73 74 75 77 80 91 99 109 116 117 121 125 129 142 153 154 160 161 168 173 178 181 182 186 187 191 196 199:204];
writetable(functional_groups_names(idx_excluded,:),'Tables/functional_groups_excluded.csv','Delimiter',',','QuoteStrings',false)
functional_groups_names(idx_excluded,:) = []; functional_groups(:,idx_excluded) = []; clear idx_excluded

% exclude functional groups not present on any analyte from the dataset
idx_not_present = find(sum(functional_groups{:,:})'==0);
writetable(functional_groups_names(idx_not_present,:),'Tables/functional_groups_not_present.csv','Delimiter',',','QuoteStrings',false)
functional_groups_names(idx_not_present,:) = []; functional_groups(:,idx_not_present) = []; clear idx_not_present

%% Filter Data
% Remove measurments with low score:
data(data.Score<95,:)=[];

% Remove analytes that have less than 42 measurments collected.
% There is (9+9+3)*4 = 84 measurements in total
[k,i,j]=unique(data.METID);
[cnt_uniquej, uniquej] = hist(j,unique(j));
idx = uniquej(cnt_uniquej<42);
data(ismember(j,idx),:)=[];
clear k i j cnt_uniquej uniquej idx

% Select measurment with the highest score (if several) 
data.METEXPID = data.METID.*100+data.EXPID;
data = sortrows(data,{'METID','METEXPID','Score'},{'ascend','ascend','ascend'});
[k,i,j]=unique(data.METEXPID,'last');
data = data(i,:);
clear i j k

% Remove dilevalol - it's repeated in the dataset
data(data.METID==72,:)=[];
%% Prepare data. Select analytes with max two dissociation steps
[~,i1,j]=unique(data.METID,'first');

data = data(cov.R(data.METID)<=2,:);

dataACD = dataACD(unique(data.METID),:);
dataNames = dataNames(unique(data.METID),:);

cov.chargesA = cov.chargesA(unique(data.METID),:);
cov.chargesB = cov.chargesB(unique(data.METID),:);
cov.charges = cov.charges(unique(data.METID),:);
cov.groupsA = cov.groupsA(unique(data.METID),:);
cov.groupsB = cov.groupsB(unique(data.METID),:);
cov.pKaslit = cov.pKaslit(unique(data.METID),:);
cov.pKasliterror = cov.pKasliterror(unique(data.METID),:);
cov.R = cov.R(unique(data.METID));
cov.logP = cov.logP(unique(data.METID),:);

functional_groups=functional_groups(unique(data.METID),:);
clear i1 j
%
cov.maxR = max(cov.R);
cov.charges  = cov.charges(:,1:cov.maxR+1);
cov.chargesA = cov.chargesA(:,1:cov.maxR+1);
cov.chargesB = cov.chargesB(:,1:cov.maxR+1);
cov.pKaslit  = cov.pKaslit(:,1:cov.maxR);
cov.groupsA  = cov.groupsA(:,1:cov.maxR);
cov.groupsB  = cov.groupsB(:,1:cov.maxR);
cov.pKasliterror = cov.pKasliterror(:,1:cov.maxR);

% exclude functional groups not present on any analyte from the dataset
idx_not_present = find(sum(functional_groups{:,:})'==0);
writetable(functional_groups_names(idx_not_present,:),'Tables/functional_groups_not_present_2.csv','Delimiter',',','QuoteStrings',false)
functional_groups_names(idx_not_present,:) = []; functional_groups(:,idx_not_present) = []; clear idx_not_present

%% Functional groups
[SortedSum,I] = sort(sum(functional_groups{:,:}>0.5));
figure('Color',[1 1 1])
subplot(1,2,1)
plot(1:1:30,SortedSum([1:1:30]),'-o')
xlabel('Functional group')
ylabel('                                                                      Number of analytes having at least one functional group of a given type')
view(90,90)
set(gca,'Xtick',[1:1:30],'XTickLabelRotation',0,'XTickLabel',functional_groups_names{I([1:1:30]),2})
set(gca,'Yscale','lin','FontSize',8)
subplot(1,2,2)
plot(31:1:60,SortedSum([31:1:60]),'-o')
view(90,90)
set(gca,'Xtick',[31:1:60],'XTickLabelRotation',0,'XTickLabel',functional_groups_names{I([31:1:60]),2})
set(gca,'Yscale','log','FontSize',8)
clear I SortedSum 
 savefig('Figures/FunctionalGroups.fig')
 set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
 print -dtiff -r300 Figures/FunctionalGroups.tif
 
%% Plot raw data (6 selected analytes)
uMETID=unique(data.METID);

uMETID_sample = [8 9 17 33 58 180]; % uMETID_sample = uMETID;

for i= 1:length(uMETID_sample);
Names = data.Name(data.METID==uMETID_sample(i));
plot_data(data,uMETID_sample(i))

annotation(gcf,'textbox',...
    [0.382142857142856 0.959328318066538 0.269642849639058 0.0369357038212866],...
    'String',Names(1),...
    'HorizontalAlignment','center',...
    'LineStyle','none');

h2 = findall(0,'type','axes'); set(h2,'ylim', [0 max(max(cell2mat(get(h2,'ylim'))))]); 

 savefig(['Figures/Individual/RawData' Names{1} '.fig'])
 set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
 print(gcf,['Figures/Individual/RawData' Names{1} '.tiff'],'-dtiff','-r300')

 close(gcf)
end

clear uMETID Names h2 i ktore uMETID_sample
%% Plot raw data (al analytes)
uMETID=unique(data.METID);
figure('Color',[1 1 1]);

for i=1:1:length(uMETID) 
plot_data(data,uMETID(i))
end

h1 = findobj(gcf,'Type', 'line'); set(h1,'LineStyle','-') 
h2 = findall(0,'type','axes'); set(h2,'ylim', [0 300])
%save
savefig('Figures/RawData.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/RawData.tif

clear uMETID h1 h2 i
%% Initialize variables and parameters
nObs = length(data.METID);
nAnalytes = length(unique(data.METID));
npH = length(unique(data.pH));

[~,i1,j]=unique(data.METID,'first');
[~,~,pHid]=unique(data.pH,'first');

% steps: MeOH (4-step aproximation), ACN  (10-step aproximation)
datastruct = struct(...
    'nAnalytes', nAnalytes, ...
    'nObs',nObs, ...
    'npH',npH, ...
    'analyte',j,...
    'pHid',pHid,...
    'steps',4.*(1-data.Mod2) + 10.*(data.Mod2),...
    'hplcparam',[data.tg data.td data.to data.te data.fio data.fik data.Mod2+1 data.pHo data.alpha1 data.alpha2 (data.Temp-25)/10],...   
    'mod', data.Mod2+1, ...
    'logPobs',cov.logP, ...
    'maxR',cov.maxR,...
    'R',cov.R,...
    'pKaslit',cov.pKaslit, ...
    'pKasliterror',cov.pKasliterror, ...
    'groupsA',cov.groupsA, ...
    'groupsB',cov.groupsB, ...
    'chargesA',cov.chargesA,...
    'chargesB',cov.chargesB,...
    'K', size(functional_groups,2),...
    'nrfungroups',functional_groups{:,:},...
	'trobs', data.RT);

clear i1 j expid npH pHid
%% Initialize
clear init0
% Initialize the values for each variable in each chain
for i=1:4
    S.logkwHat  =  normrnd(2.2,2,1);
    S.S1mHat     = normrnd(4,1,1) ;
	S.S1aHat     = normrnd(5,1,1) ;
    S.dlogkwHat = normrnd([-1 -1],0.125,1,2) ;
    S.dSmHat   = normrnd([0 0],0.5,1,2) ;
    S.dSaHat   = normrnd([0 0],0.5,1,2) ;
    S.S2mHat    = lognrnd(log(0.2),0.05,1,1) ; 
    S.S2aHat    = lognrnd(log(2),0.05,1,1) ; 
    S.beta  = normrnd([0.75 0.5 0.5],0.125,1,3) ;
    S.alphaAHat = normrnd([2 2],0.2,1,2) ;
    S.alphaBHat = normrnd(-[1 1],0.2,1,2) ;
    S.dlogkTHat  = normrnd(-0.087,0.022,1, 1);
    S.omegadlogkT  = lognrnd(log(0.022),0.2,1, 1);
    S.apH  = normrnd(0,0.1,1,2);
    S.sigma   = lognrnd(log(0.2),0.2,1, datastruct.nAnalytes);
    S.msigma  = lognrnd(log(0.2),0.2,1, 1);
    S.ssigma  = lognrnd(log(0.5),0.2,1, 1);
    S.omega = [1 1 1] .* exp(normrnd(0, 0.5, 1, 3));
    S.rho1 = [1 0.75 0.75
             0.75 1 0.75 
             0.75 0.75 1];
    S.L2 = [1 0
            0.75 0.6614];
    S.kappa = [0.25 0.25 0.25] .* exp(normrnd(0, 0.2, 1, 3));
    S.tau   = [0.5 0.5] .* exp(normrnd(0, 0.2, 1, 2));
    S.pilogkw = zeros(1,datastruct.K);
    S.piS1m = zeros(1,datastruct.K);
    S.piS1a = zeros(1,datastruct.K);
    S.sdpi = [0.1 0.1 0.1] .* exp(normrnd(0, 0.1, 1, 3));
    S.param =  [2+1.*datastruct.logPobs 4*ones(datastruct.nAnalytes,1)+0.5.*datastruct.logPobs 5*ones(datastruct.nAnalytes,1)+0.5.*datastruct.logPobs]; 
    S.dlogkwA = -1.*ones(datastruct.nAnalytes,datastruct.maxR+1);
    S.dlogkwB = -1.*ones(datastruct.nAnalytes,datastruct.maxR+1);
    S.dSmA = 0.*ones(datastruct.nAnalytes,datastruct.maxR+1);
    S.dSmB = 0.*ones(datastruct.nAnalytes,datastruct.maxR+1);
    S.dSaA = 0.*ones(datastruct.nAnalytes,datastruct.maxR+1);
    S.dSaB = 0.*ones(datastruct.nAnalytes,datastruct.maxR+1);
    S.dlogkT = normrnd(-0.0868,0.0217, 1, datastruct.nAnalytes);
    S.pKaw = datastruct.pKaslit;
    S.etaStd1 =zeros(2,datastruct.nAnalytes);
    S.etaStd2 =zeros(2,datastruct.nAnalytes);
    init0(i) = S;
end
clear S i i1 j kaHat kwHat nAnalytes nObs fi nExp
%% Use Stan for optimization
setenv('STAN_NUM_THREADS','6')
fprintf( 'Running Stan...\n' );
fito= stan('file','hplc-gra-redsum-qsrr-L.stan','data', datastruct,'method','optimize', ...
              'working_dir','Tmpstan','verbose', logical(1),'init',init0(1),'iter',1000,'warmup',1000, ...  
              'stan_home', 'C:\Users\biofarm\Documents\.cmdstanr\cmdstan-2.25.0');
fito.block()
save('fito.mat', 'fito','-v7.3')
%% Use Stan for prior predcitve check
% set environmentla variables
setenv('STAN_NUM_THREADS','6')
fprintf( 'Running Stan...\n' );
fitp= stan('file','hplc-gra-redsum-qsrr-L-priors.stan','data', datastruct, 'verbose', logical(1), ...
              'working_dir','Tmpstan','iter',1000,'warmup',100,'chains',4,'init',init0, ...
              'stan_home', 'C:\Users\biofarm\Documents\.cmdstanr\cmdstan-2.25.0');
fitp.block();
save('fitp.mat', 'fitp','-v7.3')
%% Use Stan for sampling 
% set environmentla variables
setenv('STAN_NUM_THREADS','6')
% set initial values for param based on opitmization
for i =1:4; init0(i).param =  fito.sim.samples.param; end
fprintf( 'Running Stan...\n' );
fit= stan('file','hplc-gra-redsum-qsrr-L.stan','data', datastruct, 'verbose', logical(1), ...
              'working_dir','Tmpstan','iter',1000,'warmup',1000,'chains',4,'init',init0, ...
              'stan_home', 'C:\Users\biofarm\Documents\.cmdstanr\cmdstan-2.25.0');
fit.block(); 
%% Summary of model parameters. Save to file
diary hplc-gra-redsum-qsrr-L.txt
fit.print(); 
diary off
save hplc-gra-redsum-qsrr-L.mat '-v7.3'
%%  Extract samples
samples = extract_stan_samples
save('samples-hplc-gra-redsum.mat', 'samples','-v7.3')
%% Simulate for better graphics
hplcparam_sim = readtable('Data\hplcparam_design.csv');
samples_sim = hplc_gra_sim(samples,datastruct,hplcparam_sim{:,:});
save('samples-hplc-gra-redsum-sim.mat', 'samples_sim','-v7.3')
%% Load saved data
hplcparam_sim = readtable('Data\hplcparam_design.csv');
load hplc-gra-redsum-qsrr-L.mat
load samples-hplc-gra-redsum.mat % 
load samples-hplc-gra-redsum-sim.mat
%% Goodness of Fit Plots, GOF
trPred_mean  = mean(samples.trPred);
trCond_mean  = mean(samples.trCond);

figure('Color', [1 1 1]);
subplot(3,1,1)
hold on
gscatter(trPred_mean,datastruct.trobs',datastruct.analyte)
xlabel('Population predicted t_{R,z}')
ylabel('Observed t_{R,z}')
plot(xlim,xlim,':')
xlim([0 350]);
ylim([0 350]);
legend off
subplot(3,1,2)
hold on
gscatter(trCond_mean,datastruct.trobs',datastruct.analyte)
plot(xlim,xlim,':')
xlabel('Individual Predicted t_{R,z}')
ylabel('Observed t_{R,z}')
legend off
xlim([0 350]);
ylim([0 350]);
subplot(3,1,3)
hold on
gscatter(data.EXPID ,datastruct.trobs'-trCond_mean,datastruct.analyte)
plot(xlim,[0 0],':')
xlabel('Experiment ID')
ylabel({'Residuals'})
legend off
xlim([1 84]);

savefig('Figures/GOF.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/GOF.tif

clear logkCond_mean logkPred_mean
clear trCond_mean trPred_mean
%% Trace plots:
samples_np = fit.extract('permuted',false);
Param = 'S2aHat';
[~,m]=size(samples_np(1).(Param));
for i=1:min(m,10)
figure('Color', [1 1 1]);
for z=1:4
hold on
plot(samples_np(z).(Param)(:,i),'-');
xlabel('Iteration')
end
ylabel([Param '(:,' num2str(i) ')'],'fontsize',12);
end
clear samples_np Param z i n m
%% Marginal posterior and prior distributions for the population-level parameters
load fitp.mat
samplesp =  fitp.extract;

h1 = figure;
Names = {'logkwHat' 1 '\theta_{logkwN}'
         'S1mHat' 1 '\theta_{S1mN}'
         'S1aHat' 1 '\theta_{S1aN}'
         'S2mHat' 1 '\theta_{S2m}'
         'S2aHat' 1 '\theta_{S2a}'
         'dlogkwHat' 1 '\theta_{dlogkwA}'
         'dlogkwHat' 2 '\theta_{dlogkwB}'
         'dSmHat' 1 '\theta_{dSmA}'
         'dSmHat' 2 '\theta_{dSmB}'
         'dSaHat' 1 '\theta_{dSaA}'
         'dSaHat' 2 '\theta_{dSaB}'
         'dlogkTHat' 1 '\theta_{dlogkT}'
         'beta' 1 '\beta_{logkwN}'
         'beta' 2 '\beta_{S1mN}'
         'beta' 3 '\beta_{S1aN}'
         'alphaAHat' 1 '\alpha_{mA}'
         'alphaAHat' 2 '\alpha_{aA}'
         'alphaBHat' 1 '\alpha_{mB}'
         'alphaBHat' 2 '\alpha_{aB}'
         'omega' 1 '\omega_{logkwN}'
         'omega' 2 '\omega_{S1mN}'
         'omega' 3 '\omega_{S1aN}'
         'rho1' 2 '\rho_1_{ [logkwN,S1mN]}'
         'rho1' 3 '\rho_1_{ [logkwN,S1aN]}'
         'rho1' 6 '\rho_1_{ [S1aN,S1mN]}'
         'omegadlogkT' 1 '\omega_{dlogkT}'
         'kappa' 1 '\kappa_{dlogkw}'
         'kappa' 2 '\kappa_{dSm}'
         'kappa' 3 '\kappa_{dSa}'
         'tau' 1 '\tau_{m}'
         'tau' 2 '\tau_{a}'
         'rho2' 2 '\rho_2_{ [\alpham,\alphaa]}'
         'apH' 1 'apH_{A}'
         'apH' 2 'apH_{B}'
         'msigma' 1 'm_{\sigma}'
         'ssigma' 1 's_{\sigma}'}
    
clear Param
for i=1:size(Names,1)
   temp = samplesp.(Names{i,1});
   Param(:,i) = squeeze(temp(:,Names{i,2}));
end

hold on
boxplot_pwhisker(Param,{'Labels',Names(:,3)},5,95);
set(gca, 'TickLabelInterpreter', 'tex');
view(90,90)
set(gca,'FontSize',10)
set(gca,'Position', [0.2343    0.1100    0.6707    0.8150])
ylabel('Marginal prior/posterior distributions','FontSize',10)

h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),[0.8 0.8 0.8],'FaceAlpha',.5,'EdgeColor',[0.8 0.8 0.8]);
end

h = findobj(gca,'type','line');
set(h,'Color',[0.8 0.8 0.8])

ax1=gca;
%  savefig('Figures/PriorDistribution.fig')
%  set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
%  print -dtiff -r300 Figures/PriorDistribution.tif

 clear temp Param
 
% Add posterior
h2 = figure
       
clear Param
for i=1:size(Names,1)
   temp = samples.(Names{i,1});
   Param(:,i) = temp(:,Names{i,2});
end

hold on
boxplot_pwhisker(Param,{'Labels',Names(:,3)},5,95);
set(gca, 'TickLabelInterpreter', 'tex');
view(90,90)
set(gca,'FontSize',10)
set(gca,'Position', [0.2343    0.1100    0.6707    0.8150])
ylabel('Marginal prior/posterior distributions','FontSize',10)

h = findobj(gca,'Tag','Box');
for j=1:1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),'b','FaceAlpha',.5,'EdgeColor','b');
end

h = findobj(gca,'type','line');
set(h,'Color','b')

ax2 = gca; ax2Chil = ax2.Children; % Get handles for all children from ax2
copyobj(ax2Chil, ax1);% Copy all ax2 objects to axis 1
close(h2)
set(ax1,'Ylim',[-1.5 6.5])

savefig('Figures/PosteriorDistribution.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18]) 
print -dtiff -r300 Figures/PosteriorDistribution.tif

%  
clear temp Param Names h1 h2 h ax1 ax2 i j ax2Chil
%% Calculate individual (analyte-specific) parameters

idata.logP =cov.logP;

  for j = 1:1:datastruct.nAnalytes
      samples.etap(:,j,1) = samples.param(:,j, 1) - samples.miu(:,j,1); % 
      samples.etap(:,j,2) = samples.param(:,j, 2) - samples.miu(:,j,2); % 
      samples.etap(:,j,3) = samples.param(:,j, 3) - samples.miu(:,j,3); % 
  end
  
idata.etap = squeeze(mean(samples.etap));    % logkw, S1, S2
idata.param = squeeze(mean(samples.param));  % logkw, S1, S2

idata.logkwx = squeeze(mean(samples.logkwx));
idata.logkmx = squeeze(mean(samples.logkwx-samples.S1mx));
idata.logkax = squeeze(mean(samples.logkwx-samples.S1ax));

idata.pKaw = squeeze(mean(samples.pKaw));
idata.alpham = squeeze(mean(samples.alpham));
idata.alphaa = squeeze(mean(samples.alphaa));

idata.pKam=squeeze(mean(samples.pKaw+samples.alpham));
idata.pKaa = squeeze(mean(samples.pKaw+samples.alphaa));

  for j = 1:1:datastruct.nAnalytes
 samples.dlogkw(:,j,:) =  (squeeze(samples.dlogkwA(:,j,:))) .* datastruct.chargesA(j,:) ...
                        + (squeeze(samples.dlogkwB(:,j,:))) .* datastruct.chargesB(j,:);
 samples.dS1m(:,j,:) =  squeeze(samples.dSmA(:,j,:)) .* datastruct.chargesA(j,:) ...
                       + squeeze(samples.dSmB(:,j,:)) .* datastruct.chargesB(j,:) ;
 samples.dS1a(:,j,:) = (  squeeze(samples.dSaA(:,j,:))) .* datastruct.chargesA(j,:) ...
                        + (squeeze(samples.dSaB(:,j,:))).* datastruct.chargesB(j,:) ; 
  end

 samples.dlogkm = samples.dlogkw - samples.dS1m;
 samples.dlogka = samples.dlogkw - samples.dS1a;

 idata.dlogkw = squeeze(mean(samples.dlogkw))./cov.charges;
 idata.dlogkm = squeeze(mean(samples.dlogkm))./cov.charges;
 idata.dlogka = squeeze(mean(samples.dlogka))./cov.charges;

 idata.dS1m = squeeze(mean(samples.dS1m))./cov.charges;
 idata.dS1a = squeeze(mean(samples.dS1a))./cov.charges;
 
 idata.chargesAB = cov.chargesB-cov.chargesA; % {-2,-1,0,1,2}
 idata.groupsAB = cov.groupsB-cov.groupsA;    % {-1 for Acids, 1 for Bases}

 idata.isdiss = 0.*cov.charges;
 for j = 1:1:datastruct.nAnalytes
    idata.isdiss(j,1:cov.R(j)+1)=1;
 end

clear j
%% Individual Parameters - Neutral Form
figure('Color', [1 1 1]);
xynames = {'logkwN_{i}','S1mN_{i}','S1aN_{i}','logP_i'};
gplotmatrix([idata.param(:,1) idata.param(:,2) idata.param(:,3) idata.logP],[],0*idata.logP,'kk',[],[],'on','stairs',xynames,xynames)

 h=get(gcf,'children');
 set(h(1),'Visible','off')
 savefig('Figures/IndividualParametersNeutralForm.fig')
 set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
 print -dtiff -r300 Figures/IndividualParametersNeutralForm.tif
clear h xynames
%% Effect of dissociation
figure('Color', [1 1 1]);
xynames = {'dlogkw_{r,i}','dS1m_{r,i}','dS1a_{r,i}'};
ktore = idata.isdiss(:)~=0;
gplotmatrix([idata.dlogkw(ktore) idata.dS1m(ktore) idata.dS1a(ktore)],[],idata.chargesAB(ktore),'rrybb',[],[],'on','stairs',xynames,xynames)

 h=get(gcf,'children');
 set(h(1),'Visible','off')
savefig('Figures/IndParamEffectDiss.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/IndParamEffectDiss.tif
clear h xynames ktore
%% pKas 
xynames = {'pKaw_{i,r}','pKam_{i,r}','pKaa_{i,r}','pKawlit_{i,r}'};
X=[idata.pKaw(:) idata.pKam(:) idata.pKaa(:) datastruct.pKaslit(:)]; X(X>12)=NaN; X(X<2)=NaN;
Y = idata.groupsAB(:);
[h,ax,bigax] = gplotmatrix(X(Y~=0,:),[],Y(Y~=0),'rbgkym',[],[],'on','stairs',xynames,xynames)
 set(ax,'XTick',2:2:12,'YTick',2:2:12)
 set(ax(1:4,1),'YTickLabel',2:2:12)
 set(ax(end,:),'XTickLabel',2:2:12)
 savefig('Figures/IndParampKas.fig')
  h=get(gcf,'children');
 set(h(1),'Visible','off')
 set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
 print -dtiff -r300 Figures/IndParampKas.tif
 clear h ax bigax xynames X Y

%% Sigmas
hist(mean(log(samples.sigma)))
xlabel('\sigma_i')
%% Effect of temperature
hist(mean(samples.dlogkT))
xlabel('dlogkT_i')
%% Influence of functional groups
figure('Color', [1 1 1]);
subplot(1,4,2)
hold on
boxplot_pwhisker(samples.pilogkw(:,:),{'Labels',functional_groups_names{:,2}},5,95);
plot(xlim,[0 0],':')
ylim([-1 1])
view(90,90)
set(gca,'FontSize',5)
set(gca,'Position', [0.2139    0.1100    0.2138    0.8150])
ylabel('\pi_{logkwN}','FontSize',8)
subplot(1,4,3)
hold on
boxplot_pwhisker(samples.piS1m(:,:),{'Labels',functional_groups_names{:,1}},5,95);
plot(xlim,[0 0],':')
ylim([-1 1])
view(90,90)
set(gca,'FontSize',5)
set(gca,'Position', [0.4854    0.1100    0.2178    0.8150])
ylabel('\pi_{S1mN}','FontSize',8)
subplot(1,4,4)
hold on
boxplot_pwhisker(samples.piS1a(:,:),{'Labels',functional_groups_names{:,1}},5,95);
plot(xlim,[0 0],':')
ylim([-1 1])
view(90,90)
set(gca,'FontSize',5)
set(gca,'Position', [0.7334    0.1100    0.1708    0.8150])
ylabel('\pi_{S1aN}','FontSize',8)

savefig('Figures/FunctionalGroupEffects.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/FunctionalGroupEffects.tif

%% Individual and population predictions (6 selected analytes)
prediction_type  = 'trObsCond'; % trObsPred || trObsCond

metidx = [8 9 17 33 58 180];
idx = find(ismember(unique(data.METID),metidx));
Names = dataNames{ismember(dataNames{:,1},metidx),2};

for i=1:length(metidx);
    
figure('Color', [1 1 1]);

plot_data(data,metidx(i))
plot_sim(samples_sim,hplcparam_sim,idx(i),prediction_type)

annotation(gcf,'textbox',...
    [0.382142857142856 0.959328318066538 0.269642849639058 0.0369357038212866],...
    'String',Names(i),...
    'HorizontalAlignment','center',...
    'LineStyle','none');

h2 = findall(0,'type','axes'); set(h2,'ylim', [0 max(max(cell2mat(get(h2,'ylim'))))]); 

savefig(['Figures/Individual/' prediction_type Names{i} '.fig'])
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print(gcf,['Figures/Individual/' prediction_type Names{i} '.tiff'],'-dtiff','-r300')

close(gcf)
end

clear h2 idx i metidx prediction_type Names

%% Individual and population predictions (all analytes)
prediction_type  = 'trObsCond'; % trObsPred || trObsCond

metidx = unique(data.METID);
idx = find(ismember(unique(data.METID),metidx));
Names = dataNames{ismember(dataNames{:,1},metidx),2};

for i=1:length(metidx);
    
figure('Color', [1 1 1]);

plot_data(data,metidx(i))
plot_sim(samples_sim,hplcparam_sim,idx(i),prediction_type)

annotation(gcf,'textbox',...
    [0.382142857142856 0.959328318066538 0.269642849639058 0.0369357038212866],...
    'String',Names(i),...
    'HorizontalAlignment','center',...
    'LineStyle','none');

h2 = findall(0,'type','axes'); set(h2,'ylim', [0 max(max(cell2mat(get(h2,'ylim'))))]); 

savefig(['Figures/All/' prediction_type Names{i} '.fig'])
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print(gcf,['Figures/All/' prediction_type Names{i} '.tiff'],'-dtiff','-r300')

close(gcf)
end

clear h2 idx i metidx prediction_type Names

%% Uncertainty chromatograms (population and individual predcitions) for expid = 47
map = colormap('lines');

figure('Color', [1 1 1]);
metidx = [8 9 17 33 58 180];
idx = find(ismember(unique(data.METID),metidx));
expid = [47];
expidx = find(hplcparam_sim.expid==expid);
subplot(2,1,1)
plot_uncertainity_chromatogram(samples_sim.trObsPred,expidx,idx)
tr = data.RT(ismember(data.METID,metidx)&data.EXPID==expid);
for i=1:length(tr); plot([tr(i) tr(i)], ylim,':','Color',map(i,:)); end
xlim([0 35])
title('Population Predictions')
subplot(2,1,2)
plot_uncertainity_chromatogram(samples_sim.trObsCond,expidx,idx)
tr = data.RT(ismember(data.METID,metidx)&data.EXPID==expid);
for i=1:length(tr); plot([tr(i) tr(i)], ylim,':','Color',map(i,:)); end
xlim([0 35])
legend1=legend(dataNames{ismember(dataNames{:,1},metidx),2})
set(legend1,...
    'Position',[0.723690473363513 0.155566707459585 0.20964285996982 0.17151941480726],...
    'EdgeColor',[1 1 1]);
title('Individual Predictions')
xlabel('t_R, min')
clear tr idx expid expidx metidx i map legend1
%% Limited Data predicitons. Initialize variables and parameters
metidx = [8 9 17 33 58 180];
expidx = [10 58 38];
hplcparam_sim((ismember(hplcparam_sim.expid,expidx)),:)
idx = find(ismember(unique(data.METID),metidx));
% Predictions based on selected measurments
data_red  = data(ismember(data.METID,metidx),:);
data_red(~ismember(data_red.EXPID,expidx),:)=[];

% Functional group effects
pilogkw = [0.0889968428595548 0.177471413490907;0.220089668251250 0.105968524228844;-0.0771815195591502 0.160325195284994;-0.0581799102143500 0.185308852260919;-0.0201489460052250 0.179988182045063;0.0307316725699825 0.188612904909725;0.0543621082853000 0.120183861760715;0.0324279323432999 0.182588420787220;-0.0559667868367500 0.171610973282035;-0.0415813552520000 0.127933551854748;-0.0180587034066250 0.180318487129278;-0.233834863103000 0.119557887197776;0.180337220964750 0.0846757399085927;0.222922921009750 0.108777856957033;-0.000501192305000001 0.131647570428321;0.0498735976307500 0.146940875915100;0.0513917052316499 0.132389150840822;0.287699338262500 0.111017971807696;0.124327135344325 0.0609063261016972;0.125937048681750 0.173498858874232;-0.0935018665000727 0.119702270923968;0.110744323446425 0.155684606307299;0.0988604463651751 0.154538855196298;0.0143156623489500 0.102276959128834;0.172177835703000 0.132273861408408;-0.0643482931813277 0.149964201815220;-0.0329482021954501 0.151168437643277;0.423955068749999 0.0885875918043772;0.0993908985895301 0.114247980097966;0.0298201324193250 0.182464612224784;0.262131928062750 0.162515143173171;0.125264699990745 0.139539749944581;0.108892420331975 0.106826055058630;0.0343301783096001 0.0934536247436581;0.143383295807107 0.151605033204463;-0.101369205041500 0.147919735532570;-0.0913377286217774 0.104111219920434;0.114594420303100 0.100462396653018;0.137289133807375 0.172043654639451;-0.159722669718600 0.165543368574018;-0.100180180674275 0.123600934733901;-0.172452433209500 0.143307896625322;-0.151108470096749 0.161005113067004;-0.0716024258115749 0.147898507889060;-0.410842843382499 0.140738135836306;0.000433213414249995 0.175875588073256;0.0476187438239748 0.0922153601330513;0.00506268302377501 0.182678051432591;0.0826298697423501 0.167920419044898;-0.00806687702650001 0.182219831488352;0.0879922134395501 0.181998003635673;0.160473326353290 0.166060658904532;-0.00545335401175000 0.180062747952312;0.0491115158007500 0.148776460100283;-0.0790377194335501 0.152635415304690;0.0804110515920003 0.186671844754081;0.0556822471320525 0.132266474042620;0.100240845222148 0.109456712990130;-0.0393841636146501 0.187559603506411;0.0723208501620699 0.182642665955130];
piS1m=[0.0217000559264250 0.154957706487726;0.106134236236000 0.103414295235257;0.0778085607394998 0.145230648283856;0.0655436159672500 0.159787278072099;0.00522921446779999 0.160644821033471;-0.00690202034879999 0.169823560721250;0.208429806273775 0.116246517233544;-0.0118745781861125 0.159477133102472;0.100122203701200 0.154500566094926;0.0644238794738526 0.120777670402343;0.0397175303883750 0.156209567867953;-0.00329354588375000 0.124583226054949;0.197713133829750 0.0897640046925205;0.0623477881566750 0.108635663842868;-0.00735698785372498 0.133018819334084;0.0655856143230000 0.135373389640694;-0.0676902528240000 0.127087050009747;0.147878240230000 0.110581200578122;0.244202023450001 0.0675280219418342;-0.0442913221941248 0.150922383973941;0.173501527884500 0.113705291932268;-0.0365348548789000 0.138043420187343;-0.105239627925450 0.145407961486631;0.0317578123057000 0.101285112390947;-0.00892374372052501 0.119905477353245;0.0149970674375000 0.135600616114253;-0.0448266668852501 0.130958829322927;0.301479356809751 0.0906550705894951;0.0594514271803000 0.112569001623926;0.0431497501975751 0.168617953959701;0.0485574303400775 0.140023445586041;0.0769855942661501 0.128243679718009;0.145813218572550 0.107365014749321;0.0279250855683500 0.0932651788037824;0.0303524569695750 0.129411974521634;0.178533757770950 0.139406109232078;-0.0842705905559500 0.104531826953638;0.154872312716500 0.0996907059774276;0.132141557586992 0.152110361897885;-0.0405228301690000 0.157263720169405;-0.140858090681658 0.118152826876805;0.156149595732875 0.133036941759429;0.134440491507875 0.138817968012840;0.139985368892750 0.133760811454953;-0.201352129304250 0.130151112570089;0.103264498479375 0.158431454882315;0.176517505312900 0.107967439587904;-0.142850163107775 0.176517622964707;0.0640129549986750 0.150302065640594;0.0308196197489499 0.163472973013075;0.0309345101660500 0.154618187499486;0.0439905966312500 0.141494845812313;0.0335519910824000 0.160873678513828;-0.0415121552322499 0.129682894598821;0.0865000891908250 0.135453887389572;-0.0418402961457725 0.162897804074594;0.107768039835750 0.128110626150208;0.144802405712934 0.113721157394385;0.193430843130750 0.169572416239922;0.0939918425327248 0.162166456557561];
piS1a =[0.0698588119030001 0.266997955870408;-0.0144656829931750 0.148733155306423;-0.0732117244065000 0.227686995189610;-0.161391291165675 0.275901409852163;0.0338917363063001 0.271877410824601;-0.0126832724365000 0.298346663477855;0.0203286872035000 0.177806814002461;0.111895183797900 0.271521660214816;-0.0993323177924252 0.241436538993382;0.0941978176669001 0.180301450303532;-0.0231392560610000 0.274114484568687;0.462180803630000 0.188957918025456;0.596185300750001 0.126331207619497;0.0747870440726349 0.160085514634547;0.215890568380500 0.201692251523540;0.0111592813609250 0.214000122466492;-0.174233959451578 0.189218817192267;0.301355481089251 0.159306379551736;0.370040966425000 0.0916786440894036;0.0844944851305248 0.246941048748585;0.216316038603325 0.168665961619620;0.0142375086965000 0.218564066263434;0.0210094238269000 0.223810574360652;-0.112519956732200 0.140485010893325;0.147743851832600 0.182825344199810;0.301887597099050 0.213382465113893;0.0615851100933250 0.203977937016960;0.350459336138750 0.127671681824427;0.196511031060970 0.162785467231247;-0.207562572348100 0.319149603000405;-0.240331138083249 0.215819497157118;-0.710979029074999 0.204438130939904;0.0522476095115002 0.153576361682664;0.0590822581768752 0.134308875041780;-0.196431292645900 0.195638042452994;-0.100291397251750 0.204207737155929;-0.197336832944275 0.145253622132466;0.144002680948375 0.145300767608426;-0.0980488794257999 0.264017658442989;0.474931054038000 0.248261603630060;-0.0508391006226000 0.172241816537243;0.564059040740000 0.210695763630808;0.0679908603261500 0.226609867181402;-0.351798112383883 0.206893336516015;0.142693050118950 0.187933820047448;-0.247756557571000 0.262197661810480;0.538594641175001 0.145058363009612;0.589749770952500 0.303034971268726;-0.0625016822134999 0.255047107307086;0.0519043844712500 0.283839486679962;-0.116257019861325 0.274911738051457;0.0535596520790248 0.222526855386509;0.0484423938490250 0.291964169639225;0.170557791367350 0.200459681659015;-0.230579356470800 0.217900375881488;0.0139409652806000 0.267260001573277;-0.106784199815388 0.193770940997105;-0.566092862437250 0.157219504776589;-0.400852275487326 0.285199083543308;-0.112792394288325 0.274935172540649];

clear metidx expidx idx
%% Exclude unnecesary data
nObs = length(data_red.METID);
nAnalytes = length(unique(data_red.METID));

npH = length(unique(data_red.pH));

[~,i1,j]=unique(data_red.METID,'first');
[~,~,pHid]=unique(data_red.pH,'first');

% steps: MeOH (4-step aproximation), ACN  (10-step aproximation)
datastruct_small = struct(...
    'nAnalytes', nAnalytes, ...
    'nObs',nObs, ...
    'npH',npH, ...
    'analyte',j,...
    'pHid',pHid,...
    'steps',4.*(1-data_red.Mod2) + 10.*(data_red.Mod2),...
    'hplcparam',[data_red.tg data_red.td data_red.to data_red.te data_red.fio data_red.fik data_red.Mod2+1 data_red.pHo data_red.alpha1 data_red.alpha2 (data_red.Temp-25)/10],...   
    'mod', data_red.Mod2+1, ...
    'logPobs',cov.logP(idx), ...
    'maxR',cov.maxR,...
    'R',cov.R(idx,:),...
    'pKaslit',cov.pKaslit(idx,:), ...
    'pKasliterror',cov.pKasliterror(idx,:), ...
    'groupsA',cov.groupsA(idx,:), ...
    'groupsB',cov.groupsB(idx,:), ...
    'chargesA',cov.chargesA(idx,:),...
    'chargesB',cov.chargesB(idx,:),...
    'mpilogkw', pilogkw(:,1)', ...
    'spilogkw', pilogkw(:,2)', ...
    'mpiS1m', piS1m(:,1)', ...
    'spiS1m', piS1m(:,2)', ...
    'mpiS1a', piS1a(:,1)', ...
    'spiS1a', piS1a(:,2)', ...
    'K', size(functional_groups,2),...
    'nrfungroups',functional_groups{idx,:},...
    'trobs', data_red.RT, ...
    'run_estimation', 1);

clear i1 j expid npH pHid nObs nAnalytes
%% Initialize
clear init0_simple
% Initialize the values for each variable in each chain
for i=1:4
    S.logkwHat  =  normrnd(3.1543,0.0871,1);
    S.S1mHat     = normrnd(4.5141,0.0971,1) ;
	S.S1aHat     = normrnd(5.5600,0.1348,1) ;
    S.dlogkwHat = normrnd([-0.7357,-0.9359],[0.0588,0.0461],1,2) ;
    S.dSmHat   = normrnd([0.3311,0.1098],[0.1030,0.0704],1,2) ;
    S.dSaHat   = normrnd([0.8910,-0.4577],[0.1057,0.0670],1,2) ;
    S.S2mHat    = normrnd(0.3741,0.0250,1,1) ; 
    S.S2aHat    = normrnd(0.8194,0.0342,1,1) ; 
    S.beta  = normrnd([0.7442,0.3553,0.3895],[0.0313,0.0393,0.0513],1,3) ;
    S.alphaAHat = normrnd([1.9736,2.1454],[0.1703,0.1844],1,2) ;
    S.alphaBHat = normrnd([-1.0005,-0.8940],[0.1405,0.1641],1,2) ;
    S.dlogkTHat  = normrnd(-0.0946,0.0026,1, 1);
    S.omegadlogkT  = normrnd(0.0344,0.0021,1, 1);
    S.apH  = normrnd([-0.0238,0.0851],[0.0010, 0.0009],1,2);
    S.sigma   = min(3.9,lognrnd(log(0.3671),1 ,1, datastruct_small.nAnalytes));
    S.msigma  = normrnd(0.3671, 0.0278,1, 1);
    S.ssigma  = normrnd(0.9989,0.0536,1, 1);
    S.omega = [0.6150,0.6762,0.9206] .* exp(normrnd([0.0393,0.0469,0.0631], 0.5, 1, 3));
    S.rho1 = [1 0.7811 0.7135 
             0.7811 1 0.9148 
             0.7135  0.9148 1];
    S.L2 = [1 0
            0.9408 0.3355];
    S.kappa = [0.5305,0.5526,0.5409] .* exp(normrnd([0.0276,0.0447,0.0422], 0.2, 1, 3));
    S.tau   = [2.2561,2.5580] .* exp(normrnd([0.1644,0.1841], 0.2, 1, 2));
    S.pilogkw = normrnd(pilogkw(:,1)',pilogkw(:,2)',1,datastruct_small.K);
    S.piS1m = normrnd(piS1m(:,1)',piS1m(:,2)',1,datastruct_small.K);
    S.piS1a = normrnd(piS1a(:,1)',piS1a(:,2)',1,datastruct_small.K);
    S.sdpi = [0.1964,0.1736,0.3162] .* exp(normrnd([0.0301,0.0288,0.0387], 0.1, 1, 3));
    S.param =  [2+0.75.*datastruct_small.logPobs 4*ones(datastruct_small.nAnalytes,1)+0.3.*datastruct_small.logPobs 5*ones(datastruct_small.nAnalytes,1)+0.3.*datastruct_small.logPobs]; 
    S.dlogkwA =  -0.7357.*ones(datastruct_small.nAnalytes,datastruct_small.maxR+1);
    S.dlogkwB = -0.9359.*ones(datastruct_small.nAnalytes,datastruct_small.maxR+1);
    S.dSmA = 0.3311.*ones(datastruct_small.nAnalytes,datastruct_small.maxR+1);
    S.dSmB = 0.1098.*ones(datastruct_small.nAnalytes,datastruct_small.maxR+1);
    S.dSaA = 0.8910.*ones(datastruct_small.nAnalytes,datastruct_small.maxR+1);
    S.dSaB = -0.4577.*ones(datastruct_small.nAnalytes,datastruct_small.maxR+1);
    S.dlogkT = normrnd(-0.0946,0.0217, 1, datastruct_small.nAnalytes);
    S.pKaw = datastruct_small.pKaslit;
    S.etaStd1 =zeros(2,datastruct_small.nAnalytes);
    S.etaStd2 =zeros(2,datastruct_small.nAnalytes);
    init0_simple(i) = S;
end
clear S i i1
%% Use Stan for sampling 
setenv('STAN_NUM_THREADS','1')
fprintf( 'Running Stan...\n' );
fit_small= stan('file','hplc-gra-redsum-qsrr-L-fixed.stan','data', datastruct_small, 'verbose', logical(1), ...
              'working_dir','Tmpstan','iter',1000,'warmup',1000,'chains',4,'init',init0_simple, ...
              'stan_home', 'C:\Users\biofarm\Documents\.cmdstanr\cmdstan-2.25.0');
fit_small.block(); 
%% Save
diary parameters_small.txt
fit_small.print(); 
diary off
%% Extract samples
samples_small_10_58_38 =  fit_small.extract;
samples_small_10_58_38 = hplc_gra_sim(samples_small_10_58_38,datastruct_small,hplcparam_sim{:,:});
save('samples_small_10_58_38.mat', 'samples_small_10_58_38', '-v7.3')
%%
load('samples_small_10_58_38.mat')
samples_small_sim=samples_small_10_58_38;
%% Individual predictions based on preliminary data
metidx = [8 9 17 33 58 180];

idx = find(ismember(unique(data_red.METID),metidx));
Names = dataNames{ismember(dataNames{:,1},metidx),2};

for i=1:length(metidx);
    
figure('Color', [1 1 1]);

plot_data(data,metidx(i))
plot_sim(samples_small_sim,hplcparam_sim,idx(i),'trObsCond')
annotation(gcf,'textbox',...
    [0.382142857142856 0.959328318066538 0.269642849639058 0.0369357038212866],...
    'String',Names(i),...
    'HorizontalAlignment','center',...
    'LineStyle','none');

h2 = findall(0,'type','axes'); set(h2,'ylim', [0 max(max(cell2mat(get(h2,'ylim'))))]); 

savefig(['Figures/Individual/LimDatatrObsCond' Names{i} '.fig'])
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print(gcf,['Figures/Individual/LimDatatrObsCond' Names{i} '.tiff'],'-dtiff','-r300')

close(gcf)
end

clear h2 idx Name legend1

%% Uncertainity chromatograms (selected)
map = colormap('lines');

expid_i = hplcparam_sim.expid(hplcparam_sim.tg==30&hplcparam_sim.Temp==0);

for i=1:length(expid_i);
    
expid = expid_i(i);    
figure('Color', [1 1 1]);
metidx = [8 9 17 33 58 180];
idx = find(ismember(unique(data.METID),metidx));
expidx = find(hplcparam_sim.expid==expid);
subplot(3,1,1)
plot_uncertainity_chromatogram(samples_sim.trObsPred,expidx,idx)
tr = data.RT(ismember(data.METID,metidx)&data.EXPID==expid);
for i=1:length(tr); plot([tr(i) tr(i)], ylim,':','Color',map(i,:)); end
xlim([0 35])
title('Population Predictions')
subplot(3,1,3)
plot_uncertainity_chromatogram(samples_sim.trObsCond,expidx,idx)
tr = data.RT(ismember(data.METID,metidx)&data.EXPID==expid);
for i=1:length(tr); plot([tr(i) tr(i)], ylim,':','Color',map(i,:)); end
xlim([0 35])
legend1=legend(dataNames{ismember(dataNames{:,1},metidx),2})
set(legend1,...
    'Position',[0.723690473363513 0.155566707459585 0.20964285996982 0.17151941480726],...
    'EdgeColor',[1 1 1]);

temp = hplcparam_sim(expidx,[1 7 8 11]);

if temp.Temp==0 & temp.Mod==1
xlabel(sprintf('t_R, min [tg=%d, MeOH, pH=%1.1f, 25^oC]',hplcparam_sim{expidx,[1 8]}))
end
if temp.Temp==0 & temp.Mod==2
xlabel(sprintf('t_R, min [tg=%d, ACN, pH=%1.1f, 25^oC]',hplcparam_sim{expidx,[1 8]}))
end
if temp.Temp==1 & temp.Mod==1
xlabel(sprintf('t_R, min [tg=%d, MeOH, pH=%1.1f, 35^oC]',hplcparam_sim{expidx,[1 8]}))
end
if temp.Temp==1 & temp.Mod==2
xlabel(sprintf('t_R, min [tg=%d, ACN, pH=%1.1f, 35^oC]',hplcparam_sim{expidx,[1 8]}))
end
title('Individual Predictions')
subplot(3,1,2)
idx = find(ismember(unique(data_red.METID),metidx));
expidx = find(hplcparam_sim.expid==expid);
plot_uncertainity_chromatogram(samples_small_sim.trObsCond,expidx,idx);
tr = data.RT(ismember(data.METID,metidx)&data.EXPID==expid);
for i=1:length(tr); plot([tr(i) tr(i)], ylim,':','Color',map(i,:)); end
xlim([0 35])
ylabel('Uncertainity chromatogram, pdf')
title('Limited Data Predictions')

savefig(['Figures/UncertainityChromatograms/UncertainityChromatogram-' num2str(expid) '.fig'])
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print(gcf,['Figures/UncertainityChromatograms/UncertainityChromatogram-' num2str(expid) '.tiff'],'-dtiff','-r300')

close(gcf)

end

clear expid_i idx legend1 Names
%%
ver
% MATLAB Version: 9.2.0.556344 (R2017a)
% MATLAB License Number: 261217
% Operating System: Microsoft Windows 10 Pro Version 10.0 (Build 19043)
% Java Version: Java 1.7.0_60-b19 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode
% ----------------------------------------------------------------------------------------------------
% MATLAB                                                Version 9.2         (R2017a)
% Bioinformatics Toolbox                                Version 4.8         (R2017a)
% Curve Fitting Toolbox                                 Version 3.5.5       (R2017a)
% Global Optimization Toolbox                           Version 3.4.2       (R2017a)
% MATLAB Compiler                                       Version 6.4         (R2017a)
% MATLAB Compiler SDK                                   Version 6.3.1       (R2017a)
% Optimization Toolbox                                  Version 7.6         (R2017a)
% Parallel Computing Toolbox                            Version 6.10        (R2017a)
% SimBiology                                            Version 5.6         (R2017a)
% Statistics and Machine Learning Toolbox               Version 11.1        (R2017a)
% Symbolic Math Toolbox                                 Version 7.2         (R2017a)