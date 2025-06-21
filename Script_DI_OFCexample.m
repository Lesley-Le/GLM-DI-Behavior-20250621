
load('Data_DI_OFCexample.mat');
Fs = 1;
StimuDelay_duration = 1.5;
Baseline_duration = 1;
PSTH_bin_inSec = 0.1;
%% Section 1: calculate PSTH

index_trigger = TrialStart_recording(Go_noLaser);
index_PSTH = [];
index_meanFR_StimuDelay1s = [];
for i = 1:length(index_trigger)
    Spk = SpikeTrain-index_trigger(i);
    Spk = Spk(Spk>=-Baseline_duration & Spk<StimuDelay_duration);
    index_PSTH(i,:) = histcounts(Spk,-Baseline_duration:PSTH_bin_inSec:StimuDelay_duration)/PSTH_bin_inSec;
    index_meanFR_StimuDelay1s(i,1) = length(find(Spk>=0.5&Spk<StimuDelay_duration));
end
Go_noLaser_PSTH = index_PSTH;
Go_noLaser_meanFR = index_meanFR_StimuDelay1s;

index_trigger = TrialStart_recording(Go_Laser);
index_PSTH = [];
index_meanFR_StimuDelay1s = [];
for i = 1:length(index_trigger)
    Spk = SpikeTrain-index_trigger(i);
    Spk = Spk(Spk>=-Baseline_duration & Spk<StimuDelay_duration);
    index_PSTH(i,:) = histcounts(Spk,-Baseline_duration:PSTH_bin_inSec:StimuDelay_duration)/PSTH_bin_inSec;
    index_meanFR_StimuDelay1s(i,1) = length(find(Spk>=0.5&Spk<StimuDelay_duration));
end
Go_Laser_PSTH = index_PSTH;
Go_Laser_meanFR = index_meanFR_StimuDelay1s;

index_trigger = TrialStart_recording(NoGo_noLaser);
index_PSTH = [];
index_meanFR_StimuDelay1s = [];
for i = 1:length(index_trigger)
    Spk = SpikeTrain-index_trigger(i);
    Spk = Spk(Spk>=-Baseline_duration & Spk<StimuDelay_duration);
    index_PSTH(i,:) = histcounts(Spk,-Baseline_duration:PSTH_bin_inSec:StimuDelay_duration)/PSTH_bin_inSec;
    index_meanFR_StimuDelay1s(i,1) = length(find(Spk>=0.5&Spk<StimuDelay_duration));
end
NoGo_noLaser_PSTH = index_PSTH;
NoGo_noLaser_meanFR = index_meanFR_StimuDelay1s;

index_trigger = TrialStart_recording(NoGo_Laser);
index_PSTH = [];
index_meanFR_StimuDelay1s = [];
for i = 1:length(index_trigger)
    Spk = SpikeTrain-index_trigger(i);
    Spk = Spk(Spk>=-Baseline_duration & Spk<StimuDelay_duration);
    index_PSTH(i,:) = histcounts(Spk,-Baseline_duration:PSTH_bin_inSec:StimuDelay_duration)/PSTH_bin_inSec;
    index_meanFR_StimuDelay1s(i,1) = length(find(Spk>=0.5&Spk<StimuDelay_duration));
end
NoGo_Laser_PSTH = index_PSTH;
NoGo_Laser_meanFR = index_meanFR_StimuDelay1s;


%% Section 2: calculate DI
Laseroff_signal = [Go_noLaser_meanFR;NoGo_noLaser_meanFR];
Laseroff_type = [ones(length(Go_noLaser_meanFR),1);2*ones(length(NoGo_noLaser_meanFR),1)];
[actualAUC,significant,significance] = Permutation_ROC_20250602(Laseroff_type,Laseroff_signal);
DI_off = 2*actualAUC-1;

Laseron_signal = [Go_Laser_meanFR;NoGo_Laser_meanFR];
Laseron_type = [ones(length(Go_Laser_meanFR),1);2*ones(length(NoGo_Laser_meanFR),1)];
[actualAUC,significant,significance] = Permutation_ROC_20250602(Laseron_type,Laseron_signal);
DI_on = 2*actualAUC-1;


%% Section 3: plot PSTH
figure('Color','w','Position',[100 100 800 300]);
TimeStamps = -Baseline_duration+PSTH_bin_inSec/2:PSTH_bin_inSec:StimuDelay_duration-PSTH_bin_inSec/2;
subplot(1,2,1);hold on;box off;
psth_Go_mean = mean(Go_noLaser_PSTH,1);
psth_Go_sem = std(Go_noLaser_PSTH,1)/sqrt(size(Go_noLaser_PSTH,1));
plot(TimeStamps,psth_Go_mean,'-b');
fill([TimeStamps fliplr(TimeStamps)],[psth_Go_mean-psth_Go_sem fliplr(psth_Go_mean+psth_Go_sem)],'b','FaceAlpha',0.1,'EdgeAlpha',0);
psth_NoGo_mean = mean(NoGo_noLaser_PSTH,1);
psth_NoGo_sem = std(NoGo_noLaser_PSTH,1)/sqrt(size(NoGo_noLaser_PSTH,1));
plot(TimeStamps,psth_NoGo_mean,'-r');
fill([TimeStamps fliplr(TimeStamps)],[psth_NoGo_mean-psth_NoGo_sem fliplr(psth_NoGo_mean+psth_NoGo_sem)],'r','FaceAlpha',0.1,'EdgeAlpha',0);
set(gca,'XLim',[-1 1.5],'box','off','fontsize',13,'ylim',[0 90]);
xlabel('Time (s)');
ylabel('Spikes/s');
title('Laser off');
legend('Go');
subplot(1,2,2);hold on;box off;
psth_Go_mean = mean(Go_Laser_PSTH,1);
psth_Go_sem = std(Go_Laser_PSTH,1)/sqrt(size(Go_Laser_PSTH,1));
plot(TimeStamps,psth_Go_mean,'-b');
fill([TimeStamps fliplr(TimeStamps)],[psth_Go_mean-psth_Go_sem fliplr(psth_Go_mean+psth_Go_sem)],'b','FaceAlpha',0.1,'EdgeAlpha',0);
psth_NoGo_mean = mean(NoGo_Laser_PSTH,1);
psth_NoGo_sem = std(NoGo_Laser_PSTH,1)/sqrt(size(NoGo_Laser_PSTH,1));
plot(TimeStamps,psth_NoGo_mean,'-r');
fill([TimeStamps fliplr(TimeStamps)],[psth_NoGo_mean-psth_NoGo_sem fliplr(psth_NoGo_mean+psth_NoGo_sem)],'r','FaceAlpha',0.1,'EdgeAlpha',0);
set(gca,'XLim',[-1 1.5],'box','off','fontsize',13,'ylim',[0 90]);
xlabel('Time (s)');
ylabel('Spikes/s');
title('Laser on');




