

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
% tags(1:5,:) = tag_temp(1:5,:); tags(6,:) = squeeze(mark_group6(1,2,:)); tags(7:11,:) = tag_temp(6:10,:); tags(12,:) = squeeze(mark_group6(2,2,:));
tags(1:5,:) = tag_temp(1:5,:); tags(6,:) = tag_temp(5,:) + 910; tags(7:11,:) = tag_temp(6:10,:); tags(12,:) = tag_temp(10,:) + 910;

half_1 = min(tags(6,:)-tags(1,:)); half_2 = min(tags(12,:)-tags(7,:));
tags_seg1(1:5,1) = (tags([1:5])-tags(1))'+1; tags_seg1(6:10,1) = (tags([7:11])-tags(7))' + half_1; tags_seg1(:,2) = tags_seg1(:,1) + 45*Fs;
tags_seg2(1:5,1) = (tags([1:5])-tags(1))' +45*Fs; tags_seg2(6:10,1) = (tags([7:11])-tags(7))' + half_1 + 45*Fs; tags_seg2(:,2) = tags_seg2(:,1) + 45*Fs; tags_seg2_resh = ceil(reshape(tags_seg2,1,20)/10);
tags_seg1_resh = []; tags_seg2_resh = [];
for i = 1 : 10
    tags_temp1 = ceil(tags_seg1/10);
    tags_seg1_resh = [tags_seg1_resh, tags_temp1(i,1):tags_temp1(i,2)];
    tags_temp2 = ceil(tags_seg2/10);
    tags_seg2_resh = [tags_seg2_resh, tags_temp2(i,1):tags_temp2(i,2)];
end

% start_mean_task = 30*Fs; stop_mean_task = 30*Fs+half_1+half_2;
start_mean_task_resam = ceil(30*Fs/10); stop_mean_task_resam = ceil((30*Fs+half_1+half_2)/10);
A_nirsdata_oxy_allgroup = zeros(CH_num,ceil((half_1+half_2)/10),sub_num); C_nirsdata_oxy_allgroup = zeros(CH_num,ceil((half_1+half_2)/10),sub_num);

for i_group = 1 : sub_num
    eval(['A_nirsdata_oxy_allgroup_10hz = (A_nirsdata_group',num2str(i_group),'.nirsdata.oxyData([tags(1,i_group):tags(1,i_group)+half_1-1 tags(7,i_group):tags(7,i_group)+half_2-1],:));']);
    eval(['C_nirsdata_oxy_allgroup_10hz = (C_nirsdata_group',num2str(i_group),'.nirsdata.oxyData([tags(1,i_group):tags(1,i_group)+half_1-1 tags(7,i_group):tags(7,i_group)+half_2-1],:));']);
    A_nirsdata_oxy_allgroup(:,:,i_group) = (resample(A_nirsdata_oxy_allgroup_10hz,1,10))';
    C_nirsdata_oxy_allgroup(:,:,i_group) = (resample((C_nirsdata_oxy_allgroup_10hz),1,10))';
end

% load rest
rootpath_restrain = '..\rest\';
for groupA = 1 : sub_num
    eval(['A_rest_nirsdata_group',num2str(groupA),' = load (strcat(rootpath_restrain,''A\3processed_nofilt\'', ''A',num2str(groupA),'_rest.mat'')); ']);
    eval(['C_rest_nirsdata_group',num2str(groupA),' = load (strcat(rootpath_restrain,''C\3processed_nofilt\'', ''C',num2str(groupA),'_rest.mat'')); ']);
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
MaxScale = 128; period_select = [25:33] + 40;% [41:49] + 40; % MaxScale = 256; period_length = 106;
time_lag = -35:35;  duan = length(time_lag);
% i_ch = 15; j_ch = 24;
i_ch = 17; j_ch = 8;
AC_INS_task_seg1 = zeros(sub_num,length(time_lag));           AC_INS_rest = zeros(sub_num,length(time_lag));
AC_INS_task_seg2 = zeros(sub_num,length(time_lag));           
% AC_sig95_task = zeros(period_length,task_timelength,CH_num,sub_num);       AC_sig95_rest = zeros(period_length,rest_timelength,CH_num,sub_num);
load persample.mat;
for i_group = 1 : sub_num
    fprintf('i_group = %n', i_group);
    tic;
    for i_lag = 1 : 71
        lagg = time_lag(i_lag);
        erwei = abs(lagg)+1;       
        % seg1
        signaltask_1 = [];  signaltask_2= []; Rsq1= [];
        signaltask_A_seg1 = A_nirsdata_oxy_allgroup(:,tags_seg1_resh,:); signaltask_C_seg1 = C_nirsdata_oxy_allgroup(:,tags_seg1_resh,:);
        longg_task = size(signaltask_C_seg1,2);  longg_rest = size(C_nirsdata_rest_allgroup,2); 
        if lagg > 0
            signaltask_2 = signaltask_C_seg1(:,erwei:end,i_group);
            signaltask_1 = signaltask_A_seg1(:,1:longg_task-abs(lagg),i_group);
            signalrest_2 = C_nirsdata_rest_allgroup(:,erwei:end,i_group);
            signalrest_1 = A_nirsdata_rest_allgroup(:,1:longg_rest-abs(lagg),i_group);
        else
            signaltask_1 = signaltask_A_seg1(:,erwei:end,i_group);
            signaltask_2 = signaltask_C_seg1(:,1:longg_task-abs(lagg),i_group);
            signalrest_1 = A_nirsdata_rest_allgroup(:,erwei:end,i_group);
            signalrest_2 = C_nirsdata_rest_allgroup(:,1:longg_rest-abs(lagg),i_group);
        end
        [Rsq1, period1, scale1, ~, ~] = wtc(signaltask_1(i_ch,:)', signaltask_2(j_ch,:)','S0',0.3, 'MaxScale', MaxScale, 'mcc', 0);
        AC_INS_task_seg1(i_group, i_lag) = mean(mean(Rsq1(period_select,:),2),1);
        [Rsqrest, ~, ~, ~, ~] = wtc(signalrest_1(i_ch,:)', signalrest_2(j_ch,:)','S0',0.3, 'MaxScale', MaxScale, 'mcc', 0);
        AC_INS_rest(i_group, i_lag) = mean(mean(Rsqrest(period_select,:),2),1);
        % seg2
        signaltask_1 = [];  signaltask_2= []; Rsq1= [];
        signaltask_A_seg2 = A_nirsdata_oxy_allgroup(:,tags_seg2_resh,:); signaltask_C_seg2 = C_nirsdata_oxy_allgroup(:,tags_seg2_resh,:);
        longg_task = size(signaltask_C_seg2,2);
        if lagg > 0
            signaltask_2 = signaltask_C_seg2(:,erwei:end,i_group);
            signaltask_1 = signaltask_A_seg2(:,1:longg_task-abs(lagg),i_group);
        else
            signaltask_1 = signaltask_A_seg2(:,erwei:end,i_group);
            signaltask_2 = signaltask_C_seg2(:,1:longg_task-abs(lagg),i_group);
        end
        [Rsq1, period1, scale1, ~, ~] = wtc(signaltask_1(i_ch,:)', signaltask_2(j_ch,:)','S0',0.3, 'MaxScale', MaxScale, 'mcc', 0);
        AC_INS_task_seg2(i_group, i_lag) = mean(mean(Rsq1(period_select,:),2),1);
    end
end

%% save
save(strcat(plotrootpath,'result_synchronization_CI_geiding_rest_INS_only_nofilt_CHwise_segment_timelag_IFG.mat'), 'AC_INS_task_seg1', 'AC_INS_task_seg2','AC_INS_rest',...
    'tags'); %'period1', 'scale1'






