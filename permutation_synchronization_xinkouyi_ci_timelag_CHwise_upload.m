


%% =====================================================================================
%% ===                                  CI                                           ===
%% =====================================================================================
% Author : Xu Duan
% Date : 2023-02-23
% reference: The effect of speech-gesture asynchrony on the neural coupling of interlocutors in interpreter-mediated communication
% Social Cognitive and Affective Neuroscience

clear; clc;

rootpath_CI = '..\CI_geiding\';
plotrootpath = strcat('..\CI_geiding\','sychron\plot_A\');% load CI_geding NIRS_KIT, badchannel_process,filter 后的数据
sub_num = 20;
% load A C
for groupA = 1 : sub_num
    eval(['A_nirsdata_group',num2str(groupA),' = load (strcat(rootpath_CI,''A\3processed_nofilt\'', ''A',num2str(groupA),'_CI_geiding.mat'')); ']);
    eval(['C_nirsdata_group',num2str(groupA),' = load (strcat(rootpath_CI,''C\3processed_nofilt\'', ''C',num2str(groupA),'_CI_geiding.mat'')); ']);
end
Fs = 1/A_nirsdata_group1.nirsdata.T;
CH_num = A_nirsdata_group1.nirsdata.nch;
% mark
for i_group = 1 : sub_num
    eval(['[mark_group1(:,2,i_group),~] = find(A_nirsdata_group',num2str(i_group),'.nirsdata.vector_onset==1);' ]);
    eval(['mark_group1(:,1,i_group) = A_nirsdata_group',num2str(i_group),'.nirsdata.vector_onset(mark_group1(:,2,i_group)); ']);%mark第一列是mark的名称，第二列是对应的时间点
    eval(['[mark_group6(:,2,i_group),aaa] = find(A_nirsdata_group',num2str(i_group),'.nirsdata.vector_onset==6);' ]);
    eval(['mark_group6(:,1,i_group) = A_nirsdata_group',num2str(i_group),'.nirsdata.vector_onset(mark_group6(:,2,i_group)); ']);%mark第一列是mark的名称，第二列是对应的时间点
end
tag_temp = squeeze(mark_group1(:,2,:));
tags(1:5,:) = tag_temp(1:5,:); tags(6,:) = squeeze(mark_group6(1,2,:)); tags(7:11,:) = tag_temp(6:10,:); tags(12,:) = squeeze(mark_group6(2,2,:));
half_1 = min(tags(6,:)-tags(1,:)); half_2 = min(tags(12,:)-tags(7,:));

% start_mean_task = 30*Fs; stop_mean_task = 30*Fs+half_1+half_2;
start_mean_task_resam = ceil(30*Fs/10); stop_mean_task_resam = ceil((30*Fs+half_1+half_2)/10);
A_nirsdata_oxy_allgroup = zeros(CH_num,ceil((half_1+30*Fs+half_2+5*Fs)/10),sub_num); C_nirsdata_oxy_allgroup = zeros(CH_num,ceil((half_1+30*Fs+half_2+5*Fs)/10),sub_num);

for i_group = 1 : sub_num
    eval(['A_nirsdata_oxy_allgroup_10hz = (A_nirsdata_group',num2str(i_group),'.nirsdata.oxyData([tags(1,i_group)-30*Fs:tags(1,i_group)+half_1-1 tags(7,i_group):tags(7,i_group)+half_2-1+5*Fs],:));']);
    eval(['C_nirsdata_oxy_allgroup_10hz = (C_nirsdata_group',num2str(i_group),'.nirsdata.oxyData([tags(1,i_group)-30*Fs:tags(1,i_group)+half_1-1 tags(7,i_group):tags(7,i_group)+half_2-1+5*Fs],:));']);
    A_nirsdata_oxy_allgroup(:,:,i_group) = (resample(A_nirsdata_oxy_allgroup_10hz,1,10))';
    C_nirsdata_oxy_allgroup(:,:,i_group) = (resample((C_nirsdata_oxy_allgroup_10hz),1,10))';
end

% load restrain
rootpath_restrain = '..\rest_rain\';
for groupA = 1 : sub_num
    eval(['A_rest_nirsdata_group',num2str(groupA),' = load (strcat(rootpath_restrain,''A\3processed_nofilt\'', ''A',num2str(groupA),'_restrain.mat'')); ']);
    eval(['C_rest_nirsdata_group',num2str(groupA),' = load (strcat(rootpath_restrain,''C\3processed_nofilt\'', ''C',num2str(groupA),'_restrain.mat'')); ']);
end
for i_group = 1 : sub_num
    eval(['point_reference(:,i_group) = find(A_rest_nirsdata_group',num2str(i_group),'.nirsdata.vector_onset==9);']);
end

start_mean_rest = 10*Fs; stop_mean_rest = 290*Fs-1;
start_mean_rest_resam = 10*Fs/10; stop_mean_rest_resam = ceil((290*Fs-1)/10);
A_nirsdata_rest_allgroup = zeros(CH_num,300*Fs/10,sub_num); C_nirsdata_rest_allgroup = zeros(CH_num,300*Fs/10,sub_num);
for i_group = 1 : sub_num
    eval(['A_nirsdata_rest_allgroup_10hz = (A_rest_nirsdata_group',num2str(i_group),'.nirsdata.oxyData(point_reference(1,i_group) : point_reference(1,i_group)+300*Fs-1,:));']);
    eval(['C_nirsdata_rest_allgroup_10hz = (C_rest_nirsdata_group',num2str(i_group),'.nirsdata.oxyData(point_reference(1,i_group) : point_reference(1,i_group)+300*Fs-1,:));']);
    A_nirsdata_rest_allgroup(:,:,i_group) = (resample((A_nirsdata_rest_allgroup_10hz),1,10))';
    C_nirsdata_rest_allgroup(:,:,i_group) = (resample((C_nirsdata_rest_allgroup_10hz),1,10))';
end

for i_group = 1 : sub_num
    eval(['clear A_nirsdata_group',num2str(i_group),' C_nirsdata_group',num2str(i_group),' A_rest_nirsdata_group',num2str(i_group),' C_rest_nirsdata_group',num2str(i_group),';']);
end

%% permutation test
perm_num = 1;
MaxScale = 128; period_length = 106; % MaxScale = 256; period_length = 106;
task_timelength = 12000; rest_timelength = 4000;   
period_selcet = (41:49)+40; CH1 = [16]; CH2 = [8]; laggs = [17];
% period_selcet = [25:33]+40; CH1 = 17; CH2 = 8; laggs = 0;
lagg = laggs; erwei = abs(lagg)+1;
template = reshape(1:676, 26,26);  
load persample.mat;   
AC_INS_task = zeros(sub_num,perm_num,length(CH1));  AC_INS_restrain = zeros(sub_num,perm_num,length(CH1));
longg_task = size(C_nirsdata_oxy_allgroup,2); longg_rest = size(C_nirsdata_rest_allgroup,2);   
%% ===============period_select===================
%% =======================CH========================
for i_sel = 1 : length(CH1) 
    for i_rand = 1 : 1000
        fprintf('i_rand = %d', i_rand);
        tic;
        for i_group = 1 :  sub_num
            %% ==================================== lagg ====================================
            signaltask_2 = C_nirsdata_oxy_allgroup(:,erwei:end,persample(i_group,2,i_rand));
            signaltask_1 = A_nirsdata_oxy_allgroup(:,1:longg_task-abs(lagg),persample(i_group,1,i_rand));
            signalrest_2 = C_nirsdata_rest_allgroup(:,erwei:end,persample(i_group,2,i_rand));
            signalrest_1 = A_nirsdata_rest_allgroup(:,1:longg_rest-abs(lagg),persample(i_group,1,i_rand));          
            %  sychronication between A and C
            % task
            [Rsq1, period1, scale1, ~, ~] = wtc(signaltask_1(CH1(i_sel),:)', signaltask_2(CH2(i_sel),:)','S0',0.3, 'MaxScale', MaxScale, 'mcc', 0);
            temp1 = zscore(mean(Rsq1(:, : ),2));
            AC_INS_task(i_group,i_rand,i_sel) = mean(temp1(period_selcet),1);
            % rest
            [Rsq2, ~, ~, ~, ~] = wtc(signalrest_1(CH1(i_sel),:)', signalrest_2(CH2(i_sel),:)','S0',0.3, 'MaxScale', MaxScale, 'mcc', 0);
            temp2 = zscore(mean(Rsq2(:, : ),2));
            AC_INS_restrain(i_group,i_rand,i_sel) = mean(temp2(period_selcet),1);
              
        end
        toc;
    end
end

%% save
save(strcat(plotrootpath,'permutation_result_synchronization_CI_geiding_rest_INS_only_nofilt_CHwise_timelag_ifg.mat'), 'AC_INS_task', 'AC_INS_restrain'); %'period1', 'scale1' 

