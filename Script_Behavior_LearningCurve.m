
    IDi=2;
    Fs = 1000;
    Go_trial_type = 1;
    Bin = 20;
    Step = 20;
%% Section 1: transfer CSV to mat
    RawData = table2array(readtable('Data_Behavior_LearningCurve.csv'));
    trial_time = RawData(:,1);
    mice_trial = RawData(:,2*IDi+1);
    mice_lick = RawData(:,2*IDi);

%% Section 2: calculate performance & plot
    index_lick = find(mice_lick>4);
    index_lick2 = find(diff(index_lick)>3);
    Time_lick_onset = [index_lick(1);index_lick(index_lick2+1)];
    Time_lick_offset = [index_lick(index_lick2);index_lick(end)];
    if mice_lick(1)>1
        Time_lick_offset(1) = [];
    elseif mice_lick(end)>1
        Time_lick_onset(end) = [];
    end
    clear index_lick index_lick2;
    
    index_trial = find(mice_trial>1);
    index_trial2 = find(diff(index_trial)>3);
    Time_trial_onset = [index_trial(1);index_trial(index_trial2+1)];
    Time_trial_onset = Time_trial_onset(1:239);
    
    % 获得一些定性的参数，例如trial类型和hit cr miss fa
    Trial_type = zeros(length(Time_trial_onset),1);
    Trial_output = zeros(length(Time_trial_onset),1);
    for iTrial = 1:length(Time_trial_onset)
        if mice_trial(Time_trial_onset(iTrial)+Fs*0.5)>4 & mice_trial(Time_trial_onset(iTrial)+Fs*1.2)>4
            Trial_type(iTrial,1) = 2;
        elseif mice_trial(Time_trial_onset(iTrial)+Fs*0.5)<3 & mice_trial(Time_trial_onset(iTrial)+Fs*1.2)<3
            Trial_type(iTrial,1) = 1;
        end
        index_trial_onset = Time_trial_onset(iTrial,1);
        index_resp_lick = find(Time_lick_onset>=(1.5*Fs+index_trial_onset)&Time_lick_onset<(index_trial_onset+3*Fs));
        if ~isempty(index_resp_lick)
            if Trial_type(iTrial,1)==Go_trial_type
                Trial_output(iTrial,1) = 1;
            else
                Trial_output(iTrial,1) = 4;
            end
        else
            if Trial_type(iTrial,1)==Go_trial_type
                Trial_output(iTrial,1) = 5;
            else
                Trial_output(iTrial,1) = 3;
            end
        end
    end
    output_nMouse(:,1) = Trial_type;
    % 可以计算学习曲线
    nBin = ceil((length(find(Trial_type>0))-Bin)/Step)+1;
    index_output =  Trial_output(find(Trial_type>0));
   for iBin = 1:nBin
       if iBin~=nBin
           Bin_output = [];
           Bin_output = index_output((iBin-1)*Step+1:(iBin-1)*Step+Bin);
           nHit = length(find(Bin_output==1));nMiss = length(find(Bin_output==5));
           nCR = length(find(Bin_output==3));nFA = length(find(Bin_output==4));
           dprime_nMouse(iBin,1) = Fun_calculate_dprime_ForMDtoOFCfiber(nHit,nHit+nMiss,nCR,nCR+nFA);
           Hit_nMouse(iBin,1) = nHit/(nHit+nMiss);
           CR_nMouse(iBin,1) = nCR/(nCR+nFA);
           FA_nMouse(iBin,1) = 1-CR_nMouse(iBin,1);
       else
           Bin_output = [];
           Bin_output = index_output((iBin-1)*Step+1:end);
           nHit = length(find(Bin_output==1));nMiss = length(find(Bin_output==5));
           nCR = length(find(Bin_output==3));nFA = length(find(Bin_output==4));
           dprime_nMouse(iBin,1) = Fun_calculate_dprime_ForMDtoOFCfiber(nHit,nHit+nMiss,nCR,nCR+nFA);
           Hit_nMouse(iBin,1) = nHit/(nHit+nMiss);
           CR_nMouse(iBin,1) = nCR/(nCR+nFA);
           FA_nMouse(iBin,1) = 1-CR_nMouse(iBin,1);
       end
   end
            

%% 
figure
plot(1:1:length(dprime_nMouse),dprime_nMouse,'-o','LineWidth',1);
xlabel('Block #');ylabel('d''');

figure
yyaxis left
plot(1:1:length(dprime_nMouse),Hit_nMouse,'-o','LineWidth',1);
xlabel('Block #');ylabel('Hit rate');
set(gca,'YLim',[-0.1 1.1]);
hold on
yyaxis right;
plot(1:1:length(dprime_nMouse),1-CR_nMouse,'-o','LineWidth',1);
set(gca,'YLim',[-0.1 1.1]);
ylabel('FA rate');

