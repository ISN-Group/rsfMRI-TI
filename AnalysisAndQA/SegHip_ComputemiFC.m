%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add paths
currloc='C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics';

% For outputs
outputpath='C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\FinalAnalysis\Outputs';
addpath(outputpath);

% For miFC
str='\Functions\mi';
addpath(strcat(currloc,str));

% For ICC, violinplot, and colors
str='\Functions';
addpath(strcat(currloc,str));

% For atlas labels
str='\FinalAnalysis\Inputs';
addpath(strcat(currloc,str));
SchaeferNames=readtable('ComboNamesSegHipp.xlsx');

load('df.mat')
[n_subj, n_sess]=size(df);
[~, num_rois]=size(df{1,1});

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADDED SECTION FOR NORMALISATION OF TIMESERIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zscore everything to have centre mean 0 and 1 SD
% Order = TI, TIst, sham, pre
for subj=1:n_subj
for sess=1:n_sess
for rr=1:num_rois
    BigStructNorm{subj,sess}(:,rr)=zscore(df{subj,sess}(:,rr));
end
end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute mi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for subj=1:n_subj
for sesh=1:n_sess

tmp=BigStructNorm{subj,sesh};

for rr1=1:num_rois
    for rr2=1:num_rois
        MI{subj,sesh}(rr1,rr2) = mutualinfo(tmp(:,rr1),tmp(:,rr2));
    end
end  

end
end

% Since this is an expensive computational bit, let's make a checkpoint
cd(outputpath)
save("ComputeMI_4.mat",'-v7.3')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=0;
%%% REMINDER - conds 1-4 = TI, TIst, sham, rsfMRI. Or, TI, TIst, post, pre
% For each ROI pair. . . 
for rr1=1:num_rois
disp(strcat('Running PermTests for rr',num2str(rr1)))
for rr2=1:rr1
if rr1 ~= rr2

    % Make a table for the pairwise miFC for each subj and
    % session
    for subj=1:n_subj
    for sesh=1:n_sess
    ROIMI(subj,sesh)=MI{subj,sesh}(rr1,rr2);
    end
    end

    % Comparison 1: pre vs post:
    h=h+1;
    [PreVsPost_P(rr1,rr2),PreVsPost_Tstat(rr1,rr2),PreVsPost_CohensD(rr1,rr2),PreVsPost_MaxT(h,1)]=PermTest(ROIMI(:,3),ROIMI(:,4));

end
end
end

% Correct for family wise error
% See if Tstat from sig p vals are outside
% the 95% CI for the MaxT distribution 
% Find the 95% CI for the MaxT
maxTDist=fitdist(PreVsPost_MaxT, 'Normal');
CI=paramci(maxTDist);
h=0;
for rr1=1:num_rois
for rr2=1:rr1
if rr1 ~= rr2 && PreVsPost_P(rr1,rr2)<0.05
h=h+1;
if PreVsPost_Tstat(rr1,rr2)<CI(1,1) || PreVsPost_Tstat(rr1,rr2)>CI(2,1)
PreVsPostLongFormat(h,1)=rr1;
PreVsPostLongFormat(h,2)=rr2;
PreVsPostLongFormat(h,3)=PreVsPost_P(rr1,rr2);
end
end
end
end

%%
% For Pre vs Post that survive corrections for multiple comparisons, see if
% there's a sig effect of TI
ResTbl=[];
h=0;

% If PreVsPost is signficiant, store stats, run PostVsTI and PostvsTIst
for ii=1:size(PreVsPostLongFormat)
rr1=PreVsPostLongFormat(ii,1);
rr2=PreVsPostLongFormat(ii,2);

% Store stats for pre vs post
h=h+1;
ResTbl(h,1)=rr1;
ResTbl(h,2)=rr2;
ResTblName(h,1)=SchaeferNames(rr1,1);
ResTblName(h,2)=SchaeferNames(rr2,1);
ResTbl(h,3)=PreVsPost_P(rr1,rr2);
ResTbl(h,4)=PreVsPost_Tstat(rr1,rr2);
ResTbl(h,5)=PreVsPost_CohensD(rr1,rr2);

for subj=1:n_subj
for sesh=1:n_sess
ROIMI(subj,sesh)=MI{subj,sesh}(rr1,rr2);
end
end

% Run and store other tests
[PostVsTI_P(rr1,rr2),PostVsTI_Tstat(rr1,rr2),PostVsTI_CohensD(rr1,rr2),PostVsTI_MaxT(h,1)]=PermTest(ROIMI(:,3),ROIMI(:,1));
ResTbl(h,6)=PostVsTI_P(rr1,rr2);
ResTbl(h,7)=PostVsTI_Tstat(rr1,rr2);
ResTbl(h,8)=PostVsTI_CohensD(rr1,rr2);

[PostVsTIst_P(rr1,rr2),PostVsTIst_Tstat(rr1,rr2),PostVsTIst_CohensD(rr1,rr2),PostVsTIst_MaxT(h,1)]=PermTest(ROIMI(:,3),ROIMI(:,2));
ResTbl(h,9)=PostVsTIst_P(rr1,rr2);
ResTbl(h,10)=PostVsTIst_Tstat(rr1,rr2);
ResTbl(h,11)=PostVsTIst_CohensD(rr1,rr2);

end

save('ComputeMI_5.mat')

clearvars -except ResTbl ResTblName
ResTblName(:,3:11)=array2table(ResTbl(:,3:11));
ResTblName.Properties.VariableNames={'Region 1','Region 2','PreVsPost p-value','PreVsPost Tstat','PreVsPost CohensD','PostVsTI11 p-value','PostVsTI11 Tstat','PostVsTI11 CohensD','PostVsTI13 p-value','PostVsTI13 Tstat','PostVsTI13 CohensD'};
writetable(ResTblName,'AppendixTable4.xlsx')