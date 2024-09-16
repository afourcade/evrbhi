%% Non parametric Cluster-based permutation testing with Fieldtrip
% Test differences in HEP between LA and HA conditions

% Fieldtrip initialization
addpath(genpath('C:\Users\Antonin\Documents\MATLAB\fieldtrip-master\'));
warning('off', 'MATLAB:rmpath:DirNotFound');
rmpath(genpath('C:\Users\Antonin\Documents\MATLAB\spm12\external\fieldtrip\')); 
ft_defaults

% Data path
data_path = 'D:\NeVRo\Data\';

% Choose condition, cluster sign, sig alpha & nb randomization
% If cond = avg, need to average HEPs over mov and nomov conditions, using
% script avg_nomov_mov.m beforehand
m_cond = 'avg'; %mov; nomov; avg
c_sign = 'neg'; %neg; pos
s_alpha = 0.025; %0.3 (max); 0.025 if two-tailed test
latency = [0.25 0.45]; % in line with previous HEP literature
numrandomization = 10000; %10000; 5000; 1000

% Load data
if strcmp(m_cond, 'nomov')
    load([data_path 'nomov\SBA\nomov_HEP_HA_LA.mat']); %nomov
elseif strcmp(m_cond, 'mov')
    load([data_path 'mov\SBA\mov_HEP_HA_LA.mat']); %mov
else
    load([data_path 'avg_mov_nomov\SBA\avg_HEP_HA_LA.mat']); %pooled
end

%% Statistical testing

% Load EEG layout and neighbours mapping
load('layout2.mat'); load('neighbours.mat');
%Initialize config
cfg = [];

% Type of testing: Choose (1) OR (2)
% (1) Cluster-based permutation testing parameters
cfg.channel     = {'all', '-HEOG', '-VEOG', '-ECG'}; % exclude EOG and ECG channels
cfg.latency = latency;
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT'; % within subnum design
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2; % at least two adjacent channels to enter consideration
cfg.neighbours = neighbours; % as defined above
cfg.tail = 0; % two-tailed test
% cfg.tail = 1; % one-tailed test HA > LA
cfg.clustertail = 0; % two-tailed test
% cfg.clustertail = 1; % one-tailed test HA > LA
cfg.alpha = s_alpha;
cfg.numrandomization = numrandomization;

% (2) (single electrode) paired t-tests, for ECG signal
% cfg.channel     = {'ECG'}; % if testing differences in ECG signal
% cfg.latency = [0.328 0.360]; % for ECG waveform testing
% cfg.avgovertime = 'no'; %'yes' if average
% cfg.parameter   = 'avg';
% cfg.method      = 'analytic';
% cfg.statistic   = 'ft_statfun_depsamplesT';
% cfg.tail = 0; % two-tailed test
% cfg.alpha       = s_alpha;
% cfg.correctm    = 'no';
% cfg.correctm = 'fdr';
% cfg.numrandomization = numrandomization;  

% prepare design matrix for comparison of two conditions
Nsubj = length(all_subjHA); % n=29
design = zeros(2,2*Nsubj);
design(1,:) = [1:Nsubj 1:Nsubj];
design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];

cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;

% Run statistics
stat= ft_timelockstatistics(cfg, all_subjHA{:}, all_subjLA{:});
%save stat stat

%% Grand average over HA and LA
% calculate the grand average for each condition
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
GA_LA        = ft_timelockgrandaverage(cfg, all_subjLA{:});
GA_HA         = ft_timelockgrandaverage(cfg, all_subjHA{:});

% calculate the grand average difference HA vs. LA
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
GA_HAvsLA    = ft_math(cfg, GA_HA, GA_LA);

%% Select electrodes which clustered 
poschans=stat.label(sum(stat.posclusterslabelmat==1,2)>0);
negchans=stat.label(sum(stat.negclusterslabelmat==1,2)>0);
% save chans poschans negchans

%% Cluster plot
cfg = [];
cfg.layout = 'layout2.mat';
cfg.alpha = s_alpha;
cfg.highlightcolorneg =[1 0 0];
% cfg.toi = [0.328 0.344 0.360];
cfg.subplotsize = ([1 1]);

ft_clusterplot(cfg,stat); colorbar;

%% Plot HA and LA time-courses of clustered  electrodes
cfg = [];
if strcmp(c_sign, 'pos')
    cfg.channel=poschans;
    pval = stat.posclusters(:).prob;
    cfg.ylim=[-0.4 0.6];
    a=[0.352 0.356]-0.001; %highlight times where sig. diff.
elseif strcmp(c_sign, 'neg')
    cfg.channel=negchans;
    pval = stat.negclusters(:).prob;
    cfg.ylim=[-1.1 0.2]; %neg
    a=[0.328 0.360]-0.001; %highlight times where sig. diff.
end

cfg.xlim =[-0.1 0.6];
% cfg.ylim=[-0.4 0.6]; %pos
% cfg.ylim=[-0.9 0.3]; %neg
cfg.linewidth= 2;
cfg.linestyle={'-' '-' ':'};

cfg.linecolor=[0.8500 0.3250 0.0980; 0 0.4470 0.7410; 0.40 0.1840 0.50];
figure; 
% figure('WindowState', 'maximized');
ft_singleplotER(cfg,GA_HA , GA_LA);

xlabel('Time (s)', 'FontSize', 25); ylabel('Potential (μV)', 'FontSize', 25);
title(['p-value: ' num2str(pval)]);
xline(0)
yline(0)
ax1=gca;
set(ax1,'XTick',-0.2:0.05:0.6,'FontSize',15)
set(ax1,'YTick',-1.1:0.1:0.2,'FontSize',15)
legend('HA', 'LA')
hold on
[id, x] = match_str(GA_HAvsLA.label, cfg.channel);
xall=GA_HA.time(GA_HA.time >=a(1) & GA_HA.time <=a(2));
y2=mean(GA_LA.avg(id, GA_LA.time >=a(1) & GA_LA.time <=a(2)));
y1=mean(GA_HA.avg(id, GA_LA.time >=a(1) & GA_LA.time <=a(2)));
yall=[y1, fliplr(y2)];
xall=[xall,sort(xall, 'descend')];
fill(xall,yall,[0.5 0.5 0.5], 'LineStyle','none')
alpha(0.4)

set(gcf,'Color','w');
% saveas(gcf,'HEPHALA_neg.svg')
% print(gcf, 'HEPHALA_neg.png', '-dpng', '-r900');
% saveas(gcf,'HEPHALA_pos.svg')
% print(gcf, 'HEPHALA_pos.png', '-dpng', '-r900');

%% Topographic plot HA vs. LA (uV)
load('layout2.mat')

cfg=[];
cfg.layout=layout;
cfg.style = 'straight'; cfg.comment= 'no'; cfg.marker='off';
cfg.highlight  = 'labels';
cfg.highlightcolor=[5 1 0];
cfg.highlightchannel=negchans;
% cfg.highlightchannel=[negchans; poschans];
cfg.zlim = [-0.41 0.41];
% cfg.zlim = 'maxmin';
%figure('units','normalized','outerposition',[0 0 1 1])
figure;
cfg.xlim = [0.328 0.360]; cfg.colorbar = 'southoutside';
cfg.comment    = 'xlim';
cfg.commentpos = 'title';
ax1=gca;
set(ax1,'FontSize',15,'Color','w')
ft_topoplotER(cfg,GA_HAvsLA);
%saveas(gcf,'topo_hepHALA.svg')
%print(gcf, 'topo_hepHALA.png', '-dpng', '-r900');

%% Topographic plot T values HA vs. LA
load('layout2.mat')

cfg=[];
cfg.parameter = 'stat';
cfg.layout=layout;
cfg.style = 'straight'; cfg.comment= 'no'; cfg.marker='off';
cfg.highlight  = 'on';
% cfg.highlightcolor=[5 1 0];
cfg.highlightchannel=negchans;
% cfg.highlightchannel=[negchans; poschans];
% cfg.zlim = [-3 3];
cfg.zlim = 'maxmin';
%figure('units','normalized','outerposition',[0 0 1 1])
figure;
cfg.xlim = [0.328 0.360]; cfg.colorbar = 'southoutside';
cfg.comment    = 'xlim';
cfg.commentpos = 'title';
ax1=gca;
set(ax1,'FontSize',15,'Color','w')
ft_topoplotER(cfg,stat);

%% OTHER PLOTS (JUST IN CASE)

%% Plot (single electrode) permutation testing results
% cfg = [];
% cfg.showlabels  = 'yes';
% cfg.style     = 'blank';
% cfg.layout    = 'layout2.mat';
% cfg.highlight = 'on';
% cfg.highlightchannel = find(stat.mask);
% cfg.comment   = 'no';
% figure; ft_topoplotER(cfg, GA_HAvsLA)
% title('Nonparametric: significant without multiple comparison correction')
% figure; ft_multiplotER(cfg, GA_HA, GA_LA, GA_HAvsLA)
% 
% %Plot electrode C4
% cfg = [];
% cfg.xlim =[-0.2 0.8];
% cfg.ylim=[-1.2 0.4];
% cfg.linewidth= 2;
% cfg.linestyle={'-' '-' ':'};
% cfg.linecolor=[0.8500 0.3250 0.0980; 0 0.4470 0.7410; 0.40 0.1840 0.50];
% cfg.channel = 'C4';
% figure; ft_singleplotER(cfg, GA_HA, GA_LA, GA_HAvsLA)
% xlabel('Time (s)', 'FontSize', 25); ylabel('Potential (uV)', 'FontSize', 25);
% xline(0)
% yline(0)
% ax1=gca;
% set(ax1,'XTick',-0.2:0.05:0.8,'FontSize',15)
% %set(ax1,'YTick',-0.4:0.4:0.8,'FontSize',15)
% set(ax1,'YTick',-0.8:0.1:0.6,'FontSize',15)
% legend('HA', 'LA', 'HAvsLA')
% 
% %Plot electrode CP2
% cfg = [];
% cfg.xlim =[-0.2 0.8];
% cfg.ylim=[-0.4 0.6];
% cfg.linewidth= 2;
% cfg.linestyle={'-' '-' ':'};
% cfg.linecolor=[0.8500 0.3250 0.0980; 0 0.4470 0.7410; 0.40 0.1840 0.50];
% cfg.channel = 'CP2';
% figure; ft_singleplotER(cfg, GA_HA, GA_LA, GA_HAvsLA)
% xlabel('Time (s)', 'FontSize', 25); ylabel('Potential (uV)', 'FontSize', 25);
% xline(0)
% yline(0)
% ax1=gca;
% set(ax1,'XTick',-0.2:0.05:0.8,'FontSize',15)
% %set(ax1,'YTick',-0.4:0.4:0.8,'FontSize',15)
% set(ax1,'YTick',-0.8:0.1:0.6,'FontSize',15)
% legend('HA', 'LA', 'HAvsLA')

%% FROM TUTORIAL

% figure;
% % define parameters for plotting
% timestep      = 0.05; %(in seconds)
% sampling_rate = all_subjLA{1,1}.fsample;
% sample_count  = length(stat.time);
% j = [latency(1):timestep:latency(2)];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
% m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in EEG samples
% 
% % get relevant values
% if strcmp(c_sign, 'pos')
%     pos_cluster_pvals = [stat.posclusters(:).prob];
%     pos_clust = find(pos_cluster_pvals < s_alpha);
%     pos       = ismember(stat.posclusterslabelmat, pos_clust);
% elseif strcmp(c_sign, 'neg')
%     neg_cluster_pvals = [stat.negclusters(:).prob];
%     neg_clust = find(neg_cluster_pvals < s_alpha);
%     neg       = ismember(stat.negclusterslabelmat, neg_clust);
% end
% 
% % First ensure the channels to have the same order in the average and in the statistical output.
% % This might not be the case, because ft_math might shuffle the order
% [i1,i2] = match_str(GA_HAvsLA.label, stat.label);
% 
% % plot
% for k = 1:(length(j)-1)
%    cfg.figure     = subplot(3,2,k);
%    cfg.xlim       = [j(k) j(k+1)];
%    cfg.zlim       = [-0.41 0.41];
%    if strcmp(c_sign, 'pos')
%         pos_int        = zeros(numel(GA_HAvsLA.label),1);
%         pos_int(i1)    = all(pos(i2, m(k):m(k+1)), 2);
%         cfg.highlightchannel = find(pos_int);
%    elseif strcmp(c_sign, 'neg')
%         neg_int        = zeros(numel(GA_HAvsLA.label),1);
%         neg_int(i1)    = all(neg(i2, m(k):m(k+1)), 2);
%         cfg.highlightchannel = find(neg_int);
%    end
%    cfg.highlight  = 'on';  
%    cfg.comment    = 'xlim';
%    cfg.commentpos = 'title';
%    cfg.layout     = 'layout2.mat';
%    ft_topoplotER(cfg, GA_HAvsLA);
% end

%% TOPOPLOT SIGNIFICANT TIMINGS

% figure;
% % define parameters for plotting
% sigwin_neg = [0.328 0.372]; %neg
% % sigwin_pos = [0.344 0.356]; %pos
% timestep      = sigwin_neg(2) - sigwin_neg(1); %(in seconds)
% sampling_rate = all_subjLA{1,1}.fsample;
% sample_count  = length(stat.time);
% j = [sigwin_neg(1):timestep:sigwin_neg(2)];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
% m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in EEG samples
% 
% % get relevant values
% if strcmp(c_sign, 'pos')
%     pos_cluster_pvals = [stat.posclusters(:).prob];
%     pos_clust = find(pos_cluster_pvals < s_alpha);
%     pos       = ismember(stat.posclusterslabelmat, pos_clust);
% elseif strcmp(c_sign, 'neg')
%     neg_cluster_pvals = [stat.negclusters(:).prob];
%     neg_clust = find(neg_cluster_pvals < s_alpha);
%     neg       = ismember(stat.negclusterslabelmat, neg_clust);
% end
% 
% % First ensure the channels to have the same order in the average and in the statistical output.
% % This might not be the case, because ft_math might shuffle the order
% [i1,i2] = match_str(GA_HAvsLA.label, stat.label);
% 
% % plot
% for k = 1:(length(j)-1)
%    cfg.figure     = subplot(3,2,k);
%    cfg.xlim       = [j(k) j(k+1)];
%    cfg.zlim       = [-0.41 0.41];
%    if strcmp(c_sign, 'pos')
%         pos_int        = zeros(numel(GA_HAvsLA.label),1);
%         pos_int(i1)    = all(pos(i2, m(k):m(k+1)), 2);
%         cfg.highlightchannel = find(pos_int);
%    elseif strcmp(c_sign, 'neg')
%         neg_int        = zeros(numel(GA_HAvsLA.label),1);
%         neg_int(i1)    = all(neg(i2, m(k):m(k+1)), 2);
%         cfg.highlightchannel = find(neg_int);
%    end
%    cfg.highlight  = 'labels';
%    cfg.highlightcolor=[1 1 1];
% %    cfg.highlight  = 'on';  
%    cfg.comment    = 'xlim';
%    cfg.commentpos = 'title';
%    cfg.layout     = 'layout2.mat';
%    ft_topoplotER(cfg, GA_HAvsLA); colorbar;
% end

%% Create electrodes layout & neighbours
% % Take biosemi template, add TP9 and TP10 manually; remove AF3, AF4, PO3,
% % PO4; saved in layout2.mat
% cfg.layout='D:\NeVRo\new_HEP_data_filtHP_0_3Hz\biosemi32.lay';
% cfg.channel     = {'all', '-AF3', '-AF4', '-PO3', '-PO4'}
% cfg.layout = 'layout2.mat';
% layout = ft_prepare_layout(cfg);
% cfg.layout=layout;
% ft_layoutplot(cfg)
%                                  
% cfg.method      = 'triangulation';
% cfg.channel     = {'all', '-HEOG', '-VEOG', '-ECG'}; % exclude EOG and ECG channels
% % cfg.feedback    = 'yes'; % show a neighbour plot
% neighbours= ft_prepare_neighbours(cfg, all_subjHA{1,1}); % define neighbouring channels
% cfg.neighbours = neighbours;
% ft_neighbourplot(cfg, all_subjHA{1,1});
% save('neighbours.mat', 'neighbours');
