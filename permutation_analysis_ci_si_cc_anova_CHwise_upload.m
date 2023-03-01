     
%% Selecting the interested frequency band using Permutation based on frequency cluster

% find bandSelect.m from 
% https://www.alivelearn.net/?p=3026
% Author : Xu Duan
% Date : 2023-02-23
% reference: The effect of speech-gesture asynchrony on the neural coupling of interlocutors in interpreter-mediated communication
% Social Cognitive and Affective Neuroscience

% refer to: Yuhang Long (Beijing Normal University)  longyuhangwork@163.com 
% ===Many THX===

%% load permutation results
clear; clc;
% ci rest
CI_geiding = load ('../CI_geiding/sychron/plot_A/permutation_result_synchronization_CI_geiding_rest_INS_only_nofilt_CHwise.mat'); % 106 * 26 * 26 * 20* 1000
% cc
SI_zhong = load ('../SI_zhong/sychron/plot_A/permutation_result_synchronization_SI_zhong_rest_INS_only_nofilt_CHwise.mat'); % 106 * 26 * 26 * 20 * 1000
% si restrain
SI_zhongying = load ('../SI_zhongying/sychron/plot_A/permutation_result_synchronization_SI_zhongying_rest_INS_only_nofilt_CHwise.mat'); % 106 * 26 * 26 * 20* 1000
CH_num = 26;
load period1.mat;   
for ii= 1 : 20
    species{ii,1} = 'INS';
end

%% reshape
aa = size(SI_zhongying.AC_INS_task);
SI_zhong_task_reshape = reshape(SI_zhong.AC_INS_task,aa(1),aa(2)*aa(3),aa(4),aa(5));
CI_geiding_task_reshape = reshape(CI_geiding.AC_INS_task,aa(1),aa(2)*aa(3),aa(4),aa(5));
CI_geiding_rest_reshape = reshape(CI_geiding.AC_INS_rest,aa(1),aa(2)*aa(3),aa(4),aa(5));
SI_zhongying_task_reshape = reshape(SI_zhongying.AC_INS_task,aa(1),aa(2)*aa(3),aa(4),aa(5));
SI_zhongying_restrain_reshape = reshape(SI_zhongying.AC_INS_restrain,aa(1),aa(2)*aa(3),aa(4),aa(5));

%%  3-way ANOVA for repeated measures ci--si--cc
%% ===real===
real.SI_zhong = SI_zhong_task_reshape(:,:,:,1); real.restrain = SI_zhongying_restrain_reshape(:,:,:,1);
real.CI_geiding = CI_geiding_task_reshape(:,:,:,1); real.rest = CI_geiding_rest_reshape(:,:,:,1);
real.SI_zhongying = SI_zhongying_task_reshape(:,:,:,1);
real.INS_zhong = real.SI_zhong - real.rest; real.INS_si = real.SI_zhongying - real.rest; real.INS_ci = real.CI_geiding - real.restrain;
p_value_real = zeros(size(real.INS_ci,1), size(real.INS_ci,2));
F_value_real = zeros(size(real.INS_ci,1), size(real.INS_ci,2));
tic;
for i = 1 : size(real.INS_ci,1)
    for j = 1 : size(real.INS_ci,2)
        tic;
        meas = [squeeze(real.INS_ci(i,j,:)), squeeze(real.INS_si(i,j,:)), squeeze(real.INS_zhong(i,j,:))];
        t = table(species,meas(:,1),meas(:,2),meas(:,3),...
            'VariableNames',{'species','meas1','meas2','meas3'});
        Meas = table([1 2 3]','VariableNames',{'Measurements'});
        rm = fitrm(t,'meas1-meas3~1','WithinDesign',Meas);
        ranovatbl = ranova(rm);
        p_value_real(i,j) = ranovatbl{1,5};
        F_value_real(i,j) = ranovatbl{1,4};
        toc;
    end
end
toc;
save('real_P_F_ci_si_zhong.mat','p_value_real','F_value_real');

%% ===re-sample===
loop = 1000;
p_value_perm = zeros(size(real.INS_ci,1), size(real.INS_ci,2),loop-1);
F_value_perm = zeros(size(real.INS_ci,1), size(real.INS_ci,2),loop-1);
Meas = table([1 2 3]','VariableNames',{'Measurements'});
perm_INS_CI = CI_geiding_task_reshape - SI_zhongying_restrain_reshape;
perm_INS_SI = SI_zhongying_task_reshape - CI_geiding_rest_reshape;
perm_INS_zhong = SI_zhong_task_reshape - CI_geiding_rest_reshape;
parfor i_loop = 2 : 1000
    fprintf('i_loop = %d', i_loop);
    tic;
    % Inter-brain neural synchronization    
    p_value_temp = zeros(size(SI_zhong_task_reshape,1), size(SI_zhong_task_reshape,2));
    F_value_temp = zeros(size(SI_zhong_task_reshape,1), size(SI_zhong_task_reshape,2));   
    perm_INS_CI_lop = squeeze(perm_INS_CI(:,:,:,i_loop)); perm_INS_SI_lop = squeeze(perm_INS_SI(:,:,:,i_loop)); perm_INS_zhong_lop = squeeze(perm_INS_zhong(:,:,:,i_loop)); 
    for i = 1 : size(SI_zhong_task_reshape,1)
        for j = 1 : size(SI_zhong_task_reshape,2)
            % tic;
            meas = [squeeze(perm_INS_CI_lop(i,j,:)), squeeze(perm_INS_SI_lop(i,j,:)), squeeze(perm_INS_zhong_lop(i,j,:))];
            t = table(species,meas(:,1),meas(:,2),meas(:,3),...
                'VariableNames',{'species','meas1','meas2','meas3'});
            rm = fitrm(t,'meas1-meas3~1','WithinDesign',Meas);
            ranovatbl = ranova(rm);
            p_value_temp(i,j) = ranovatbl{1,5};
            F_value_temp(i,j) = ranovatbl{1,4};
            % toc;
        end
    end
    p_value_perm(:,:,i_loop-1) = p_value_temp;
    F_value_perm(:,:,i_loop-1) = F_value_temp;
    toc;
end
save('perm_P_F_ci_si_zhong.mat','p_value_perm','F_value_perm');

%% !!!!!!!!!!!!pause!!!!!!!!!!!
%% Analysis
load('real_P_F_ci_si_zhong.mat');
load('perm_P_F_ci_si_zhong.mat');
% Hypo2
F_value_real = abs(F_value_real);
F_value_perm = abs(F_value_perm);

%% OS
CH_num = 26;
period = period1;
period_peri = period(period>3);
p_value_real_peri = p_value_real(period>3,:);  F_value_real_peri = F_value_real(period>3,:);
F_selected = bandSelect(p_value_real_peri, 0.05, 0.05); 
F_real_mean = zeros(1,length(F_selected.fs));
for i_sel = 1 : size(F_selected,1)
    F_real_mean(i_sel) = mean(F_value_real_peri(F_selected.fs{i_sel,1},F_selected.ch(i_sel,1)));
end

%% PS
p_value_perm_all = p_value_perm;
F_value_perm_all = F_value_perm;
loop = 1000;
loop = loop-1; F_value_perm_mean = zeros(1,loop);
p_value_perm_all_peri = p_value_perm_all(period>3,:,:);  F_value_perm_all_peri = F_value_perm_all(period>3,:,:);
for i_lop = 1 : loop
    tic;
    T_selected = bandSelect(p_value_perm_all_peri(:,:,i_lop), 0.05, 0.05);
    numm = [];
    for i = 1 :length(T_selected.fs)
        numm(i) = length(T_selected.fs{i,1});
    end
    [aa,bb] = sort(numm,'descend');
    F_value_perm_mean(i_lop) = mean(F_value_perm_all_peri(T_selected.fs{bb(1),1},T_selected.ch(bb(1),1),i_lop));
     toc;
end
save('F_value_perm_mean_ci_si_zhong.mat','F_value_perm_mean');
load 'F_value_perm_mean_ci_si_zhong.mat';

%% ====================plot======================
figure;
y2 = prctile(sort(F_value_perm_mean), 99);
hist(sort(F_value_perm_mean),30); hold on;
h = findobj(gca,'Type','patch');
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'w';
figure;
hist(sort(F_value_perm_mean),30); h = findobj(gca,'Type','patch');
h.FaceColor = [0 0.5 0.5]; h.EdgeColor = 'w';
aa =get(gca,'XTickLabel');  loweredge = str2num(aa{1,1});
bb =get(gca,'XTickLabel'); higheredge = str2num(bb{end,1});

set(gca,'XLim',[loweredge higheredge]); % set(gca,'XTick',[-7:0.5:1]); set(gca,'XTickLabel',[-7:0.5:1]);
hold on;
limy = 140;
% grey box 2
v2 = [y2 0; higheredge 0; higheredge limy; y2 limy];
f2 = [1 2 3 4];
patch('Faces',f2,'Vertices',v2,'FaceColor',[96 96 96]/255,'FaceAlpha',.3,'EdgeColor','none');
hold on;
[~, sel] = find(F_real_mean > y2);
for i = 1 : length(sel)
    plot( [F_real_mean(sel(i)) F_real_mean(sel(i))], get(gca, 'YLim'), '-r', 'LineWidth', 0.5);
end
xlabel('F-value'); ylabel('Number of permutation samples');
% title('The distribution of 1000 PS');

F_all = F_selected(sel,2);
F_comp = [];
for i = 1:size(F_all,1)
    eval(['F.f',num2str(i),' = cell2mat(F_all.fs(i));']);
    eval(['freq_select = union(F.f',num2str(i),',F_comp)']);
    F_comp = freq_select;
end


