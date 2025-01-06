%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currloc='C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\FinalAnalysis';

% For ICC, violinplot, and colors
addpath('C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\Functions');

% For LMEM
addpath('C:\Users\dk818\OneDrive - Imperial College London\LilThingsAndMATLABPrograms\emmeans')

% For outputs
str='\Outputs';
outputpath=strcat(currloc,str);
addpath(outputpath);

% For atlas and other data
str='\Inputs';
addpath(strcat(currloc,str));

% For Motion outlier data
OutlierPath='C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\MotionOutliers';
addpath(OutlierPath);

load('ComputeMI_5.mat')
clearvars -except outputpath MI num_rois df n_subj n_sess OutlierPath ResTbl

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the effect of condition on Motion Outliers - FD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute motion outliers FD
MotionDF={};
ShamFiles=dir(strcat(OutlierPath,'\*Sham*'));
TIFiles=dir(strcat(OutlierPath,'\*TI_*'));
TIstFiles=dir(strcat(OutlierPath,'\*TIst*'));
PreFiles=dir(strcat(OutlierPath,'\*pre*'));

for subj=1:n_subj
MotionDF{subj,1}=load(ShamFiles(subj).name);
ShamMotion(:,subj)=load(ShamFiles(subj).name);

MotionDF{subj,2}=load(TIFiles(subj).name);
TIMotion(:,subj)=load(TIFiles(subj).name);

MotionDF{subj,3}=load(TIstFiles(subj).name);
TIstMotion(:,subj)=load(TIstFiles(subj).name);

MotionDF{subj,4}=load(PreFiles(subj).name);
PreMotion(:,subj)=load(PreFiles(subj).name);
end

%%%%%%%%%%%%%
% NEWLY ADDED - see if any participants have more than 20% of their volumes
% with an FD>0.25. Since everyone has 230 volumes, that'd be 46 volumes
for ii =1:n_sess
    for subj=1:n_subj
        
        % Get the subj's FDs per session
        tmp=MotionDF{subj,ii};
        if sum(tmp(:,1)>0.5) >= 46
            disp(strcat('THIS SUBJECT',num2str(subj) ,'AND SESSION SHOULD BE EXCLUDED',num2str(ii)))
        end
        
    end
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
[Motionp,~,MotionStats]=friedman(Motion,1,'off');
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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ICC
% Assess between subject reliability of each region's MI across sessions
% Helpful link: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4913118/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute ICC
for sesh=1:n_sess
for rr=1:num_rois
for subj=1:n_subj
ROIMI(:,subj)=MI{subj,sesh}(rr,:)';
end
[BetSubjRelR(sesh,rr), ~, ~, ~, ~, ~, BetSubjRelP(sesh,rr)] = ICC(ROIMI,'A-1');
end
end 
BetSubjRelR=BetSubjRelR';

% Get mean and std ICC for AppendixTable
for ii=1:n_sess
AppendixTable(ii,3)=mean(BetSubjRelR(:,ii));
AppendixTable(ii,4)=std(BetSubjRelR(:,ii));
end

% Run main effect
[Friedp,~,Friedstats]=friedman(BetSubjRelR);
FriedPartEta2=Friedstats.n/Friedstats.sigma;

% Run post hocs
for ii=1:n_sess
for jj=1:ii
[ICCRes{1,1}(ii,jj),ICCRes{2,1}(ii,jj),ICCRes{3,1}(ii,jj)]=PermTest(BetSubjRelR(:,ii),BetSubjRelR(:,jj));
end
end

% figure()
% hold on
% set(gca,'FontSize',20)
% set(gca,'FontName','Arial')
% violins = violinplot(BetSubjRelR);
% for ii=1:4
% violins(1,ii).ShowData = 1;
% violins(1,ii).ShowMean = 'yes';
% violins(1,ii).ViolinColor = C(2,:);
% end
% % doing this manually since it's a small case
% % TI and TIst, TI to Post, Post to TIst, Pre to TIst, and Pre to Post
% group={[2,1],[3,1],[3,2],[4,2],[4,3]};
% stats=[ICCRes{1,1}(2,1),ICCRes{1,1}(3,1),ICCRes{1,1}(3,2),ICCRes{1,1}(4,2),ICCRes{1,1}(4,3)];
% H=sigstar(group,stats);
% hold off
% saveas(gcf,'QA_ICCScoreViolins.png');
% close all

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TSNR
% tSNR table was computed using a shell script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in the table 
SNRdf=readtable('tSNROutput.csv');

% Delete the subjs we're not using - subj from rest
SNRdf([41:44,65:68],:)=[];

% Now I need to make a list of each condition type - col order will be TI,
% TIst, Post, Pre
Pre=table2array(SNRdf(1:n_sess:end,1));
Sham=table2array(SNRdf(2:n_sess:end,1));
TI=table2array(SNRdf(3:n_sess:end,1));
TIst=table2array(SNRdf(4:n_sess:end,1));

% Excellent, now let's make them cols in a table - 
SNRdf2(:,1:n_sess)=[TI,TIst,Sham,Pre];

% Get mean and std for AppendixTable
for ii=1:n_sess
AppendixTable(ii,5)=mean(SNRdf2(:,ii));
AppendixTable(ii,6)=std(SNRdf2(:,ii));
end

% Run a Friedman's test
[SNRP,SNRTABLE,SNRSTATS] = friedman(SNRdf2,1,'off');
SNREta2=SNRSTATS.n/SNRSTATS.sigma;
% There is a significant, main effect of run on tSNR (p=0.04, Eta2=15.49)

for ii=1:n_sess
for jj=1:ii
[tSNR{1,1}(ii,jj),tSNR{2,1}(ii,jj),tSNR{3,1}(ii,jj)]=PermTest(SNRdf2(:,ii),SNRdf2(:,jj));
end
end
% There is a significant difference between TIst and Pre (p=0.03, tstat=-6.60, cohensD=0.62)

% % Make plots
% C=linspecer(4);
% figure()
% hold on
% set(gca,'FontSize',20)
% set(gca,'FontName','Arial')
% violins = violinplot(SNRdf2);
% for ii=1:4
% violins(1,ii).ShowData = 1;
% violins(1,ii).ShowMean = 'yes';
% violins(1,ii).ViolinColor = C(2,:);
% end
% hold off
% group={[4,2]};
% stats=[ICCRes{1,1}(4,2)];
% H=sigstar(group,stats);
% saveas(gcf,'QA_tSNRViolins.png');
% close all


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the effect of condition on Number of Noise Components
% Components identidied from ICA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tbl=readtable('rsfMRI_ICAInfo.xlsx');

% Remove Subjs we don't need
Tbl([1:8,41:44,],:)=[];

% Take what we need - otherwise it's an annoying table
NumComps=table2array(Tbl(:,7));

% Now I need to make a list of each condition type - col order will be TI,
% TIst, Post, Pre
Pre=NumComps(1:n_sess:end,1);
Sham=NumComps(2:n_sess:end,1);
TI=NumComps(3:n_sess:end,1);
TIst=NumComps(4:n_sess:end,1);

% Excellent, now let's make them cols in a table - 
NumComps2(:,1:n_sess)=[TI,TIst,Sham,Pre];

% Get mean and std for AppendixTable
for ii=1:n_sess
AppendixTable(ii,7)=mean(NumComps2(:,ii));
AppendixTable(ii,8)=std(NumComps2(:,ii));
end

% Run a Friedman's test
[NumCompsP,NumCompsTABLE,NumCompsSTATS] = friedman(NumComps2,1,'off');
NumCompsEta2=NumCompsSTATS.n/NumCompsSTATS.sigma;
% There is not a significant efect of run on the number of noise components
% (p=0.09, Eta=16.15)

% % Plot
% figure()
% hold on
% set(gca,'FontSize',20)
% set(gca,'FontName','Arial')
% violins = violinplot(NumComps2);
% for ii=1:4
% violins(1,ii).ShowData = 1;
% violins(1,ii).ShowMean = 'yes';
% violins(1,ii).ViolinColor = C(2,:);
% end
% hold off
% saveas(gcf,'QA_ICAViolins.png');
% close all

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TSNR AND DEGREE
% Need to see if changes in miFC (degree) are associated with tSNR. 
% So for each contrast (e.g., pre vs post), find the degree (number of
% signficiant edges that contain that number), and put in one column, then,
% in the next, put the tSNR for that region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNRdf=readtable('tSNROutput_AllROIs.csv');
Pre=table2array(SNRdf(1:n_sess:end,3:end));
Sham=table2array(SNRdf(2:n_sess:end,3:end));
TI=table2array(SNRdf(3:n_sess:end,3:end));
TIst=table2array(SNRdf(4:n_sess:end,3:end));

% Now add the tSNR per ROI - 
% HOWEVER, need to account that ROIs 201-218 ARE NOT the same between the
% tSNR and the res table numbering - this is because of how tSNR outputs
% went - so I need to reorder 

% Rows 201-218 need this index: 
Index=[206, 205, 202, 213, 215, 217, 204, 203, 201, 212, 211, 208, 214, 216, 218, 210, 209, 207];

for ii=1:num_rois
if ii<201
% TI11 is currently subj X region - need to average across subj
ROItSNR(ii,1)=mean(TI(:,ii));
% TI13
ROItSNR(ii,2)=mean(TIst(:,ii));
% Post 
ROItSNR(ii,3)=mean(Sham(:,ii));
% Pre
ROItSNR(ii,4)=mean(Pre(:,ii));
else
% TI
ROItSNR(ii,1)=mean(TI(:,Index(ii-200)));
% TIst
ROItSNR(ii,2)=mean(TIst(:,Index(ii-200)));
% Sham
ROItSNR(ii,3)=mean(Sham(:,Index(ii-200)));
% Pre
ROItSNR(ii,4)=mean(Pre(:,Index(ii-200)));
end
end

% Degree per ROI for pre vs post
for ii=1:num_rois
CountForROI(ii,1)=ii;
CountForROI(ii,2)=(sum(ResTbl(:,1)==ii))+(sum(ResTbl(:,2)==ii));
end

% Should get rid of rois wiht no sig change
toDel=CountForROI(:,2)<1;
CountForROI(toDel,:)=[];

% How do I see whether tSNR effected the change in degree? Maybe the
% change in tSNR? So here, for every ROI, compute Pre - Post
for ii=1:size(CountForROI,1)
CountForROI(ii,3)=ROItSNR(CountForROI(ii,1),4)-ROItSNR(CountForROI(ii,1),3);
end

[r,p]=corr(CountForROI(:,2),CountForROI(:,3),'type','Spearman');
disp('Corr stats between Pre vs Post tSNR and Degree')
disp(strcat('(r=',num2str(r),',p=',num2str(p)))

% Now I need to redo this for the degree for Post vs TI
% tSNR should be the same for all ROIs though, so just need to change Count
% for ROI
ResTbl2=ResTbl;
toDel=ResTbl2(:,6)>0.05;
ResTbl2(toDel,:)=[];
clear CountForROI

% Compute Degree
for ii=1:num_rois
CountForROI(ii,1)=ii;
CountForROI(ii,2)=(sum(ResTbl2(:,1)==ii))+(sum(ResTbl2(:,2)==ii));
end

% Should get rid of rois wiht no sig change
toDel=CountForROI(:,2)<1;
CountForROI(toDel,:)=[];

% How do I see whether tSNR effected the change in degree? Maybe the
% change in tSNR? So here, for every ROI, compute Post - TI 
for ii=1:size(CountForROI,1)
CountForROI(ii,3)=ROItSNR(CountForROI(ii,1),3)-ROItSNR(CountForROI(ii,1),1);
end

[r,p]=corr(CountForROI(:,2),CountForROI(:,3),'type','Spearman');
disp('Corr stats between Post vs TI11 tSNR and Degree')
disp(strcat('(r=',num2str(r),',p=',num2str(p)))


% Now I need to redo this for the degree for Post vs TI
% tSNR should be the same for all ROIs though, so just need to change Count
% for ROI
ResTbl2=ResTbl;
toDel=ResTbl2(:,9)>0.05;
ResTbl2(toDel,:)=[];
clear CountForROI

% Compute Degree
for ii=1:num_rois
CountForROI(ii,1)=ii;
CountForROI(ii,2)=(sum(ResTbl2(:,1)==ii))+(sum(ResTbl2(:,2)==ii));
end

% Should get rid of rois wiht no sig change
toDel=CountForROI(:,2)<1;
CountForROI(toDel,:)=[];

% How do I see whether tSNR effected the change in degree? Maybe the
% change in tSNR? So here, for every ROI, compute Post - TI 
for ii=1:size(CountForROI,1)
CountForROI(ii,3)=ROItSNR(CountForROI(ii,1),3)-ROItSNR(CountForROI(ii,1),2);
end

[r,p]=corr(CountForROI(:,2),CountForROI(:,3),'type','Spearman');
disp('Corr stats between Post Vs TI13 tSNR and Degree')
disp(strcat('(r=',num2str(r),',p=',num2str(p)))
%%
% DEPRECEATED
% 
% % Compute number of outliers -  what % are outside 3x SD
% for nn=1:n_sess
%     benchmark=meanICC(nn,1)-(3*meanICC(nn,2));
%     outliercount=0;
%     h=0;
%     for ii=1:size(BetSubjRelR,1)
%     if BetSubjRelR(ii,nn)<benchmark
%         h=h+1;
%         outliercount=outliercount+1;
%         % nameRec(h,nn)=SchaeferNames(ii,2);
%     end
%     end
%     OutlierCount(1,nn)=outliercount;
% end