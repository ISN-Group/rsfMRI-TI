%{
Compute effects of EF on miFC
I need to run corrs between the EF per subj and the miFC between the pairs
significantly effected by stim 
I need to run comparisons between each condition and it's EF
%}
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: only get the hipp fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\AtlasAndOtherInputs')
addpath('C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\FinalAnalysis')
addpath('C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\Functions')
addpath('C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\FinalAnalysis\Inputs')
outputpath='C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\FinalAnalysis\Outputs';
addpath(outputpath)
load('ComputeMI_5.mat')
tbl=readtable("EFsubjects.csv");
ROINames=readtable("ComboNamesSegHipp.xlsx");
% Get rid of electrode rows
rm=[];
aa=[1,1,1,0,0,0]';
tt=height(tbl)/height(aa);
for ii=1:tt
    start=height(rm)+1;
    fin=start+5;
    rm(start:fin,1)=aa;
end
rm=logical(rm);
tbl(rm,:)=[];
% Because we're not using Subj9, remove them
tbl(91:96,:)=[];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: separate EF by stim cond
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rm2=rm(1:height(tbl),1);
TIEF=[];
TIstEF=[];

for ii=1:height(rm2)
    if rm2(ii,1)==logical(1)
        TIEF(height(TIEF)+1,1)=table2array(tbl(ii,7));   % col 7 has the mean EF
    else
        TIstEF(height(TIstEF)+1,1)=table2array(tbl(ii,7));
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: separate EF by hipp region - ant, mid, or post
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TIef=[];
TIstef=[];
leest=[];
aa=[1,2,3]';
tt=height(TIEF)/height(aa);

for ii=1:tt
    start=height(leest)+1;
    fin=start+2;
    leest(start:fin,1)=aa;
end

for ii=1:height(leest)
    place=height(TIef)+1;
    if leest(ii,1)==1
        TIef(place,1)=TIEF(ii,1);   
        TIstef(place,1)=TIstEF(ii,1);
    elseif leest(ii,1)==2
        TIef(place-1,2)=TIEF(ii,1);  
        TIstef(place-1,2)=TIstEF(ii,1);
    else
        TIef(place-1,3)=TIEF(ii,1);
        TIstef(place-1,3)=TIstEF(ii,1);
    end
end

% Now we have two tables, one for TI and one for TIst, where columns 1,2,3
% are their EF for their ant, med, and post hipp

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: redo with normalised by each pars total hipp EF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normTIef=[];
normTIstef=[];

for seg=1:width(TIef)
for ii=1:height(TIef)
    normTIef(ii,seg)=TIef(ii,seg)/(TIef(ii,1)+TIef(ii,2)+TIef(ii,3));
    normTIstef(ii,seg)=TIstef(ii,seg)/(TIstef(ii,1)+TIstef(ii,2)+TIstef(ii,3));
end
end

% Make proportion, same as Ines
normTIef=normTIef*10;
normTIstef=normTIstef*10;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADDING LATE IN THE GAME  - previous script (EF_miFC) in the NewCodeSegHip
% folder was used to compute the effect of condition on EF. Mercifully, the
% results replicate. however, Ines is asking why I went straight to
% Friedman rather than ANOVA, so now I need to test whether they're
% nonnormally distributed. 

for ii=1:size(normTIef,2)
H=kstest(normTIef(:,ii));
if H==1
    disp('Run friedmans for TI!')
end

H=kstest(normTIstef(:,ii));
if H==1
    disp('Run friedmans for TIst!')
end
end

[TIP,~,TIStats]=friedman(normTIef);
[TIstP,~,TIstStats]=friedman(normTIstef);

% Reminderr - 1:3 = ant, mid, post hipp
for ii=1:3
    for jj=1:ii
        
            [PostHocTIpval(ii,jj),~,statsTI]=signrank(normTIef(:,ii),normTIef(:,jj),'method','approximate');
            [PostHocTIstpval(ii,jj),~,statsTIst]=signrank(normTIstef(:,ii),normTIstef(:,jj),"method","approximate");
            if ii ~= jj
            PostHocTIstZval(ii,jj)=statsTIst.zval;
            PostHocTIZval(ii,jj)=statsTI.zval;
            end
            clear stats
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: See if there's a significant relationship between EF and the
% change in miFC due to stim. 
% So run through a sig table, get the deltamiFC, then run corr, then store
% stats
% Only difficulty is that I can only use miFC from the 16 subjs that we
% have EF for
% EF table Subj List = 
% 01 10 11 12 13 14 15 18 19 20 21 22 05 06 08
% After checking with GetReadyFormiFC, that row index = 
% 10 01 02 03 04 05 06 08 09 11 12 13 17 18 20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HippROIs=[204,205,206,213,214,215];
TI11Res=ResTbl(ResTbl(:,6)<0.05,:);
TI13Res=ResTbl(ResTbl(:,9)<0.05,:);
h=0;
% Get delta miFC (DmiFC) for these regions
% MI Cols = TI TIst Sham Pre
for SubjIndex=[10,1,2,3,4,5,6,8,9,11,12,13,17,18,20]
h=h+1;

% Sig TI 
for ii=1:size(TI11Res,1)
Post2TI11(h,ii)=MI{SubjIndex,3}(TI11Res(ii,1),TI11Res(ii,2))-MI{SubjIndex,1}(TI11Res(ii,1),TI11Res(ii,2));
end

for ii=1:size(TI13Res,1)
Post2TI13(h,ii)=MI{SubjIndex,3}(TI13Res(ii,1),TI13Res(ii,2))-MI{SubjIndex,2}(TI13Res(ii,1),TI13Res(ii,2));
end

end

%%
% Since TIef adn TIstef have Cols 1:3 = ant, mid, post, I'm only going to
% run the corr between the mid and ant for TI 1:1 and TI 1:3

% First TI 1:1 : corr btwn each DmiFC and EF in the medial hipp during TI
for ii=1:size(Post2TI11,2)
Post2TI11Res(ii,1:2)=TI11Res(ii,1:2);
[Post2TI11Res(ii,3),Post2TI11Res(ii,4)]=corr(Post2TI11(:,ii),normTIef(:,2),'type','Spearman');
end

% Now TI 1:3
for ii=1:size(Post2TI13,2)
Post2TI13Res(ii,1:2)=TI13Res(ii,1:2);
[Post2TI13Res(ii,3),Post2TI13Res(ii,4)]=corr(Post2TI13(:,ii),normTIstef(:,1),'type','Spearman');
end

% FDR correct
[~,~,Post2TI11Res(:,5)]=fdr(Post2TI11Res(:,4));
[~,~,Post2TI13Res(:,5)]=fdr(Post2TI13Res(:,4));

% Count up the res!
disp('Number of edges with significantly different miFC pre-task sham vs TI 1:1 where the change in miFC relates to EF in the medial hippocampus:')
disp(num2str(sum(Post2TI11Res(:,5)<0.05)))
disp('Number of edges with significantly different miFC pre-task sham vs TI 1:3 where the change in miFC relates to EF in the ant hippocampus:')
disp(num2str(sum(Post2TI13Res(:,5)<0.05)))

%%
% What about only edges that include the ROI that's getting stimulated?

% First TI 1:1 : corr btwn each DmiFC and EF in the medial hipp during TI
Post2TI112=Post2TI11(:,ismember(Post2TI11Res(:,1),HippROIs));
for ii=1:size(Post2TI112,2)
Post2TI11ResHipp(ii,1:2)=TI11Res(ii,1:2);
[Post2TI11ResHipp(ii,3),Post2TI11ResHipp(ii,4)]=corr(Post2TI112(:,ii),normTIef(:,2),'type','Spearman');
end


% Now TI 1:3
for ii=1:size(Post2TI13,2)
Post2TI13Res(ii,1:2)=TI13Res(ii,1:2);
[Post2TI13Res(ii,3),Post2TI13Res(ii,4)]=corr(Post2TI13(:,ii),normTIstef(:,1),'type','Spearman');
end

% FDR correct
[~,~,Post2TI11Res(:,5)]=fdr(Post2TI11Res(:,4));
[~,~,Post2TI13Res(:,5)]=fdr(Post2TI13Res(:,4));

% Count up the res!
disp('Number of edges with significantly different miFC pre-task sham vs TI 1:1 where the change in miFC relates to EF in the medial hippocampus:')
disp(num2str(sum(Post2TI11Res(:,5)<0.05)))
disp('Number of edges with significantly different miFC pre-task sham vs TI 1:3 where the change in miFC relates to EF in the ant hippocampus:')
disp(num2str(sum(Post2TI13Res(:,5)<0.05)))

%%
cd(outputpath)
% Organise results for tables
for ii=1:size(Post2TI11Res,1)
Post2TI11Res_NiceFormat(ii,1)=ROINames(Post2TI11Res(ii,1),1);
Post2TI11Res_NiceFormat(ii,2)=ROINames(Post2TI11Res(ii,2),1);
Post2TI11Res_NiceFormat(ii,3:5)=array2table(Post2TI11Res(ii,3:5));
end
writetable(Post2TI11Res_NiceFormat,'AppendixTable5.csv')

for ii=1:size(Post2TI13Res,1)
Post2TI13Res_NiceFormat(ii,1)=ROINames(Post2TI13Res(ii,1),1);
Post2TI13Res_NiceFormat(ii,2)=ROINames(Post2TI13Res(ii,2),1);
Post2TI13Res_NiceFormat(ii,3:5)=array2table(Post2TI13Res(ii,3:5));
end
writetable(Post2TI13Res_NiceFormat,'AppendixTable6.csv')
