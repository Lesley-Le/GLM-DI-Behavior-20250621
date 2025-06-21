
load('Data_GLM_MDneu41.mat');
%% Basic parameter settings
warning off;
Block_length= 50;
Step_length = 25;
bin_ms = 250;% ms
psth_bin_ms = 100;% ms
Fs_TXT = 1000;
Fs_recording = 30000;


%% Section1: calculate learning curve
% Input: Trigger,Basic,TrialInformation;
Trial_num=length(TrialInformation);
Block_num=floor((Trial_num-Block_length)/Step_length)+1;

Trial_type = vertcat(TrialInformation.TrialType);
Trial_output = vertcat(TrialInformation.TrialOutput);
Bin = Block_length;Step = Step_length;
nBin = ceil((length(find(Trial_type>0))-Bin)/Step)+1;
index_output = Trial_output;
dprime = [];
outcome_rate = [];
Hit_rate = [];
FA_rate = [];    
for iBin = 1:nBin
    if iBin~=nBin
        Bin_output = [];
        Bin_output = index_output((iBin-1)*Step+1:(iBin-1)*Step+Bin);
        nHit = length(find(Bin_output==1));nMiss = length(find(Bin_output==5));
        nCR = length(find(Bin_output==3));nFA = length(find(Bin_output==4));
        dprime(iBin,1) = calculate_dprime(nHit,nHit+nMiss,nCR,nCR+nFA);
        outcome_rate(iBin,1) =(nHit/(nHit + nMiss) + nCR/(nCR +nFA)) /2;
        Hit_rate(iBin,1) =nHit/(nHit +nMiss);
        FA_rate(iBin,1) =nFA/(nCR +nFA);            
    else
        Bin_output = [];
        Bin_output = index_output((iBin-1)*Step+1:end);
        nHit = length(find(Bin_output==1));nMiss = length(find(Bin_output==5));
        nCR = length(find(Bin_output==3));nFA = length(find(Bin_output==4));
        dprime(iBin,1) = calculate_dprime(nHit,nHit+nMiss,nCR,nCR+nFA);
        outcome_rate(iBin,1) =(nHit/(nHit + nMiss) + nCR/(nCR +nFA)) /2;
        Hit_rate(iBin,1) =nHit/(nHit +nMiss);
        FA_rate(iBin,1) =nFA/(nCR +nFA);              
    end
end
clear Bin_output iBin index_output nCR nFA nHit nMiss Trial_num Trial_output Trial_type;

[~,nBin_to_max_dprime] = max(dprime);
if nBin == nBin_to_max_dprime
    nTrial_to_max_dprime = length(TrialInformation);
else
    nTrial_to_max_dprime = (nBin_to_max_dprime-1)*Step+Bin;
end

Block_num = nBin_to_max_dprime;
Trial_num = nTrial_to_max_dprime;

if ~isempty(find(FA_rate==0,1))
    if find(FA_rate==0,1)-1<=Block_num
        Block_num = find(FA_rate==0,1)-1;
        Trial_num = (Block_num-1)*Step+Bin;
    end
end

figure('Color','w','Position',[100 100 400 300])
plot(1:Block_num,dprime(1:Block_num),'-o');
xlabel('Block #');ylabel('d''');
Behav_perf.dprime = dprime;
Behav_perf.FA_rate = FA_rate;
Behav_perf.Hit_rate = Hit_rate;
Behav_perf.Block_num = Block_num;
clearvars -except Basic Bin bin_ms Behav_perf Block_length Fs_TXT psth_bin_ms Step Step_length TrialInformation Trigger Fs_recording...
    SpikeTrain SpikeAligned BasicER;
%% Section2: plot Spike PSTH
Go_psth_ave_ea_block = [];
Go_psth_err_ea_block =[];
NoGo_psth_ave_ea_block =[];
NoGo_psth_err_ea_block = [];
for bln=1:Behav_perf.Block_num
    Block_start_trial = 1+(bln-1)*Step_length;
    Block_end_trial = (bln-1)*Step_length+Block_length;
    Trial_seq = Block_start_trial:1:Block_end_trial;
    Trial_seq_type = vertcat(TrialInformation(Trial_seq).TrialType);

    Go_seq = Trial_seq(Trial_seq_type==2);
    NoGo_seq = Trial_seq(Trial_seq_type==1);

    %
    index_trial = Go_seq;
    index_psth = [];
    for iTrial = 1:length(index_trial)
        spk = SpikeAligned(index_trial(iTrial)).spikeTimeRelative/Fs_recording;
        spk_psth = histcounts(spk,-1:(psth_bin_ms/1000):3)/(psth_bin_ms/1000);
        index_psth = [index_psth;spk_psth];
    end
    Go_psth_nTrial{1,bln} = index_psth;
    Go_psth_ave_ea_block(bln,:) = mean(index_psth,1);
    Go_psth_err_ea_block(bln,:) = std(index_psth,1)/sqrt(length(Go_seq));
    % 
    index_trial = NoGo_seq;
    index_psth = [];
    for iTrial = 1:length(index_trial)
        spk = SpikeAligned(index_trial(iTrial)).spikeTimeRelative/Fs_recording;
        spk_psth = histcounts(spk,-1:(psth_bin_ms/1000):3)/(psth_bin_ms/1000);
        index_psth = [index_psth;spk_psth];
    end
    NoGo_psth_nTrial{1,bln} = index_psth; 
    NoGo_psth_ave_ea_block(bln,:) = mean(index_psth,1);
    NoGo_psth_err_ea_block(bln,:) = std(index_psth,1)/sqrt(length(NoGo_seq));
            
end
max3 = max([Go_psth_ave_ea_block;NoGo_psth_ave_ea_block],[],'all');
max3 = 80;
% Plot PSTH
figure('color','w','Position', [100 100 1400 800]);
for j = 1:Behav_perf.Block_num

    subplot(3,4,j);
    hold on;
    he1=errorbar_2([1:40], Go_psth_ave_ea_block(j,:), Go_psth_err_ea_block(j,:),'r');
    he2=errorbar_2([1:40], NoGo_psth_ave_ea_block(j,:), NoGo_psth_err_ea_block(j,:),'b');
    set(he1(2),'LineWidth',0.5);
    set(he2(2),'LineWidth',0.5);
    set(he1(1),'LineWidth',0.5);
    set(he2(1),'LineWidth',0.5);
    plot([10 10],[0 max3],'g-.');
    plot([15 15],[0 max3],'k-.');
    plot([25 25],[0 max3],'k-.');
    if j==1
        hls = legend([he1(1) he2(1)], 'Go','NoGo');
        set(hls,'fontsize',8,'Position',[0.0511 0.8434 0.0475 0.0565]);
        ylabel('Spikes/s');
    end
    xlabel('Time (s)');
    set(gca,'ylim',[0 max3],'Box','off','fontsize',11);
    set(gca,'xtick',[0:10:40],'XTickLabel',{'-1','0','1','2','3'});
    tit_txt = sprintf('d''= %0.2g', Behav_perf.dprime(j));
    ht = title(tit_txt);
    set(ht,'fontsize',12);
end

clear Block_end_trial bln Block_start_trial Go_seq NoGo_seq he1 he2 hls ht index_psth index_trial index_output;
clear iTrial j max3 spk spk_psth tit_txt Trial_seq Trial_seq_type;

PSTH.Go_psth_nBlock = Go_psth_nTrial;
PSTH.Go_psth_mean_nBlock = Go_psth_ave_ea_block;
PSTH.Go_psth_sem_nBlock = Go_psth_err_ea_block;
PSTH.NoGo_psth_nBlock = NoGo_psth_nTrial;
PSTH.NoGo_psth_mean_nBlock = NoGo_psth_ave_ea_block;
PSTH.NoGo_psth_sem_nBlock = NoGo_psth_err_ea_block;
PSTH.bin_size_in_ms = psth_bin_ms;
clear Go_psth_nTrial Go_psth_ave_ea_block Go_psth_err_ea_block NoGo_psth_nTrial NoGo_psth_ave_ea_block NoGo_psth_err_ea_block;

%% Section3: GLM
Stim_Delay_binned=(Basic.Stimulus_duration+Basic.Delay_duration)/bin_ms;
Resp_win_binned=Basic.Response_window/bin_ms;
GoCue_win_binned=Basic.GoCue_duration/bin_ms;
    
    for neu=1
        for bln=1:Behav_perf.Block_num
            Block_start_trial = 1+(bln-1)*Step_length;
            Block_end_trial = (bln-1)*Step_length+Block_length;
            Trial_seq = Block_start_trial:1:Block_end_trial;
            Trial_seq_type = vertcat(TrialInformation(Trial_seq).TrialType);
           

            nTrial_X_Y = [];
            clear Y Spk Y_binned;
            for iTrial = 1:Block_length
                %Spike
                Trial_start_ER = floor(BasicER.trigger(Trial_seq(iTrial)));%trial start是该trial开始
                Trial_end_ER = floor(BasicER.trigger(Trial_seq(iTrial)+1));%trial end是到下一个trial开始前
                Trial_dur_ER = Trial_end_ER-Trial_start_ER;
                Y=zeros(Trial_dur_ER,1);
                Spk=SpikeTrain-Trial_start_ER+1; % FR_depth_modulated_OK(neu,1),第1列是细胞id
                Spk=Spk(Spk>0&Spk<Trial_dur_ER);
                Y(Spk)=1;%
                Y_binned=fcoarse_bin(Y',bin_ms*BasicER.Fs/1000);
                nBin_Y = length(Y_binned);
                %Lick Rate
                Trial_start_txt = floor(Trigger.Trial_start(Trial_seq(iTrial)));
                Trial_end_txt = floor(Trigger.Trial_end(Trial_seq(iTrial)));
                Trial_dur_txt = Trial_end_txt-Trial_start_txt;
                LickRate=zeros(Trial_dur_txt,1);
                Lick=Trigger.Lick_absolute-Trial_start_txt+1;
                Lick=Lick(Lick>0&Lick<Trial_dur_txt);
                LickRate(Lick)=1;%
                Lick_binned=fcoarse_bin(LickRate',bin_ms*Fs_TXT/1000);
                if length(Lick_binned)~=nBin_Y
                    nBin_min = min(length(Lick_binned),nBin_Y);
                else 
                    nBin_min = nBin_Y;
                end
                Lick_binned = Lick_binned(1:nBin_min);
                Y_binned = Y_binned(1:nBin_min);
                nTrial_X_Y(iTrial).Y_binned = Y_binned;
                %go & nogo
                go_binned = zeros(nBin_min,1);nogo_binned = zeros(nBin_min,1);
                if TrialInformation(Trial_seq(iTrial)).TrialType==2 % go
                    go_binned(1:Stim_Delay_binned,1)=1;
                elseif TrialInformation(Trial_seq(iTrial)).TrialType==1 % no go
                    nogo_binned(1:Stim_Delay_binned,1)=1;
                end
                nTrial_X_Y(iTrial).go_binned = go_binned;
                nTrial_X_Y(iTrial).nogo_binned = nogo_binned;
                nTrial_X_Y(iTrial).LickRate_binned = Lick_binned;
                % outcome
                hit_binned = zeros(nBin_min,1);CR_binned = zeros(nBin_min,1);
                miss_binned = zeros(nBin_min,1);fa_binned = zeros(nBin_min,1);
                switch TrialInformation(Trial_seq(iTrial)).TrialOutput
                    case 1 %Hit
                        hit_binned(Stim_Delay_binned+1:Stim_Delay_binned+Resp_win_binned,1)=1;
                    case 3 %CR
                        CR_binned(Stim_Delay_binned+1:Stim_Delay_binned+Resp_win_binned,1)=1;
                    case 4 %FA
                        fa_binned(Stim_Delay_binned+1:Stim_Delay_binned+Resp_win_binned,1)=1;
                    case 5 %Miss
                        miss_binned(Stim_Delay_binned+1:Stim_Delay_binned+Resp_win_binned,1)=1;
                end
                nTrial_X_Y(iTrial).hit_binned = hit_binned;
                nTrial_X_Y(iTrial).CR_binned = CR_binned;
                nTrial_X_Y(iTrial).fa_binned = fa_binned;
                nTrial_X_Y(iTrial).miss_binned = miss_binned;

                % gocue
                gocue_binned = zeros(nBin_min,1);
                gocue_binned(Stim_Delay_binned+1:Stim_Delay_binned+GoCue_win_binned,1)=1;
                nTrial_X_Y(iTrial).gocue_binned = gocue_binned;
                
            end

            clear GLM_Y GLM_X;
            GLM_Y = vertcat(nTrial_X_Y.Y_binned);
            GLM_X(:,1) = vertcat(nTrial_X_Y.go_binned);
            GLM_X(:,2) = vertcat(nTrial_X_Y.nogo_binned);
            GLM_X(:,3) = vertcat(nTrial_X_Y.LickRate_binned);
            GLM_X(:,4) = vertcat(nTrial_X_Y.hit_binned);
            GLM_X(:,5) = vertcat(nTrial_X_Y.CR_binned);
            GLM_X(:,6) = vertcat(nTrial_X_Y.fa_binned);
            GLM_X(:,7) = vertcat(nTrial_X_Y.miss_binned);
            GLM_X(:,8) = vertcat(nTrial_X_Y.gocue_binned);
            c = GLM_X;
            %%%%%%%%%%%%%%%%%% GLM %%%%%%%%%%%%%%%%%%%%
            
            mdl=fitglm(c,GLM_Y,'linear','Distribution','poisson');
            
            B0{neu}(bln)=mdl.Coefficients.Estimate(1);
            B{neu}(bln,:)=mdl.Coefficients.Estimate(2:9);  % B{i}: (block个数，event个数)
            p{neu}(bln,:)=mdl.Coefficients.pValue;
            t{neu}(bln,:)=mdl.Coefficients.tStat;
            SE{neu}(bln,:)=mdl.Coefficients.SE;
            disp(['Neu',num2str(neu),'Block',num2str(bln)]);


        end
    end

clear Block_start_trial Block_end_trial c CR_binned fa_binned GLM_X GLM_Y go_binned gocue_binned iTrial;
clear Last_block_xbinned Lick Lick_binned LickRate mdl miss_binned nBin_min nBin_Y nogo_binned nTrial_X_Y;
clear p R_square SE t Trial_dur_ER Trial_dur_txt Trial_end_ER Trial_end_txt Trial_seq Trial_seq_type Trial_start_ER;
clear Trial_start_txt Y Y_binned;

GLM.B0 = B0;GLM.B = B;
clearvars -except GLM Behav_perf PSTH;
% figure
figure('Color','w','Position',[100 500 1500 300]);
X_matrix = [ones(Behav_perf.Block_num, 1) Behav_perf.dprime(1:Behav_perf.Block_num)];       
subplot(1,5,1);hold on;box off;
gd_dprime_id = find(Behav_perf.dprime(1:Behav_perf.Block_num)>1.5);
if mean(GLM.B{1,1}(gd_dprime_id,1))>mean(GLM.B{1,1}(gd_dprime_id,2))
    pref_B = GLM.B{1,1}(:,1);
    nonpref_B = GLM.B{1,1}(:,2);
else
    pref_B = GLM.B{1,1}(:,2);
    nonpref_B = GLM.B{1,1}(:,1);
end
scatter(Behav_perf.dprime(1:Behav_perf.Block_num),pref_B,'b');
scatter(Behav_perf.dprime(1:Behav_perf.Block_num),nonpref_B,'r');
[a, bint, r, rint, stats] = regress(pref_B, X_matrix);
linefit = X_matrix*a;   
plot(Behav_perf.dprime(1:Behav_perf.Block_num),linefit,'-b');
[a, bint, r, rint, stats] = regress(nonpref_B, X_matrix);
linefit = X_matrix*a;   
plot(Behav_perf.dprime(1:Behav_perf.Block_num),linefit,'-r');
set(gca,'YLim',[-0.2 1],'YTick',[0 0.4 0.8]);
xlabel('d''');ylabel('β');
legend('pref cue','nonpref cue','Location','southwest');
subplot(1,5,2);hold on;box off;
CueSI = (pref_B-nonpref_B)/mean([abs(pref_B);abs(nonpref_B)]);    
scatter(Behav_perf.dprime(1:Behav_perf.Block_num),CueSI,'k');
[a, bint, r, rint, stats] = regress(CueSI, X_matrix);
linefit = X_matrix*a;   
plot(Behav_perf.dprime(1:Behav_perf.Block_num),linefit,'-k');    
xlabel('d''');ylabel('Cue SI');  
subplot(1,5,3);hold on;box off;
FA_B = GLM.B{1,1}(:,6);Hit_B = GLM.B{1,1}(:,4);
scatter(Behav_perf.dprime(1:Behav_perf.Block_num),FA_B,'b');
scatter(Behav_perf.dprime(1:Behav_perf.Block_num),Hit_B,'r');
[a, bint, r, rint, stats] = regress(FA_B, X_matrix);
linefit = X_matrix*a;   
plot(Behav_perf.dprime(1:Behav_perf.Block_num),linefit,'-b');
[a, bint, r, rint, stats] = regress(Hit_B, X_matrix);
linefit = X_matrix*a;   
plot(Behav_perf.dprime(1:Behav_perf.Block_num),linefit,'-r');
set(gca,'YLim',[-0.2 1],'YTick',[0 0.4 0.8]);
xlabel('d''');ylabel('β');
legend('FA','Hit','Location','southwest');
subplot(1,5,4);hold on;box off;
OutcomeSI = (FA_B-Hit_B)/mean([abs(FA_B);abs(Hit_B)]);    
scatter(Behav_perf.dprime(1:Behav_perf.Block_num),OutcomeSI,'k');
[a, bint, r, rint, stats] = regress(OutcomeSI, X_matrix);
linefit = X_matrix*a;   
plot(Behav_perf.dprime(1:Behav_perf.Block_num),linefit,'-k');    
xlabel('d''');ylabel('Outcome SI');  
subplot(1,5,5);hold on;box off;
Lick_B = GLM.B{1,1}(:,3);
scatter(Behav_perf.dprime(1:Behav_perf.Block_num),Lick_B,'k');
[a, bint, r, rint, stats] = regress(Lick_B, X_matrix);
linefit = X_matrix*a;   
plot(Behav_perf.dprime(1:Behav_perf.Block_num),linefit,'-k');    
xlabel('d''');ylabel('β'); title('Lick rate');
set(gca,'YLim',[-0.2,1]);




    