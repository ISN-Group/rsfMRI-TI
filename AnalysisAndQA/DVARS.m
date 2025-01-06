% Compute DVARS
% Load in all data 
addpath('C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\FinalAnalysis\Inputs\DVARS');
DVARSpath='C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\FinalAnalysis\Inputs\DVARS';
files=dir('C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\FinalAnalysis\Inputs\DVARS\*.txt');
% For ICC, violinplot, and colors
addpath('C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\Functions');
% For outputs
outputpath='C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\FinalAnalysis\Outputs';
addpath(outputpath)

%%
CondName={'TI','TIst','SHAM','pre'};
SubjNum={'TIDES_P2T_S1' ,'TIDES_P2T_S13' ,'TIDES_P2T_S18' ,'TIDES_P2T_S21' ,'TIDES_P2T_S5' ,'TIDES_P2T_S10' ,'TIDES_P2T_S14' ,'TIDES_P2T_S19' ,'TIDES_P2T_S22' ,'TIDES_P2T_S6' ,'TIDES_P2T_S11' ,'TIDES_P2T_S15' ,'TIDES_P2T_S2' ,'TIDES_P2T_S3' ,'TIDES_P2T_S7' ,'TIDES_P2T_S12' ,'TIDES_P2T_S16' ,'TIDES_P2T_S20' ,'TIDES_P2T_S4' ,'TIDES_P2T_S8'};

for task=1:size(CondName,2)
for subj=1:size(SubjNum,2)
SubjDF=strcat(DVARSpath,'\',SubjNum{1,subj},'_',CondName{1,task},'_DVARS.txt');
SubjDF = readtable(SubjDF);
DVARStbl(subj,task)=(sum(table2array(SubjDF),'all')/230)*100;
end
end

% Now Get mean per col for hte table
mean(DVARStbl)
std(DVARStbl)

% Run stats
[DVARsp,~,DVARsStats]=friedman(DVARStbl,1);
Eta2=DVARsStats.n/DVARsStats.sigma;

% Annoying - main effect
for ii=1:size(CondName,2)
for jj=1:ii
[p(ii,jj),Tstat(ii,jj),cohensd(ii,jj),maxT(ii,jj)]=PermTest(DVARStbl(:,ii),DVARStbl(:,jj));
end
end

% Plot
C=linspecer(4);
figure()
hold on
set(gca,'FontSize',20)
set(gca,'FontName','Arial')
violins = violinplot(DVARStbl);
for ii=1:4
violins(1,ii).ShowData = 1;
violins(1,ii).ShowMean = 'yes';
violins(1,ii).ViolinColor = C(2,:);
end
% H=sigstar({[4,2],[3,2]},[0.0001,0.02])
group={[2,1],[3,2],[4,2]};
sigstats=[p(2,1),p(3,2),p(4,2)];
H=sigstar(group,sigstats);
hold off
cd(outputpath)
saveas(gcf,'QA_DVARsViolins.png');


%%

% Compute motion outliers FD

% For Motion outlier data
OutlierPath='C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\MotionOutliers';
addpath(OutlierPath);

MotionDF={};
ShamFiles=dir(strcat(OutlierPath,'\*Sham*'));
TIFiles=dir(strcat(OutlierPath,'\*TI_*'));
TIstFiles=dir(strcat(OutlierPath,'\*TIst*'));
PreFiles=dir(strcat(OutlierPath,'\*pre*'));

for subj=1:size(SubjNum,2)
MotionDF{subj,1}=load(ShamFiles(subj).name);
ShamMotion(:,subj)=load(ShamFiles(subj).name);

MotionDF{subj,2}=load(TIFiles(subj).name);
TIMotion(:,subj)=load(TIFiles(subj).name);

MotionDF{subj,3}=load(TIstFiles(subj).name);
TIstMotion(:,subj)=load(TIstFiles(subj).name);

MotionDF{subj,4}=load(PreFiles(subj).name);
PreMotion(:,subj)=load(PreFiles(subj).name);
end

MeanShamMotion=mean(ShamMotion,1);
MeanTIMotion=mean(TIMotion,1);
MeanTIstMotion=mean(TIstMotion,1);
MeanPreMotion=mean(PreMotion,1);

% Make mean and std as cols 1 & 2 for Appendix Table A.3 
% Rows 1=TI11
AppendixTable(1,1)=mean(MeanTIMotion);
AppendixTable(1,2)=std(MeanTIMotion);
% Row 2 = TI 13
AppendixTable(2,1)=mean(MeanTIstMotion);
AppendixTable(2,2)=std(MeanTIstMotion);
% Post
AppendixTable(3,1)=mean(MeanShamMotion);
AppendixTable(3,2)=std(MeanShamMotion);
% Pre
AppendixTable(4,1)=mean(MeanPreMotion);
AppendixTable(4,2)=std(MeanPreMotion);

% Run main effect
Motion=[MeanTIMotion;MeanTIstMotion;MeanShamMotion;MeanPreMotion]';
[Motionp,~,MotionStats]=friedman(Motion,1);
Eta2=MotionStats.n/MotionStats.sigma;

% % Plot
% C=linspecer(4);
% figure()
% hold on
% set(gca,'FontSize',20)
% set(gca,'FontName','Arial')
% violins = violinplot(Motion);
% for ii=1:4
% violins(1,ii).ShowData = 1;
% violins(1,ii).ShowMean = 'yes';
% violins(1,ii).ViolinColor = C(2,:);
% end
% hold off
% cd(outputpath)
% saveas(gcf,'QA_FDViolins.png');