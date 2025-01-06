%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add paths
currloc='C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics';
addpath('C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\FinalAnalysis')

% For outputs
outputpath='C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\FinalAnalysis\Outputs';
addpath(outputpath);

% For ICC, violinplot, and colors
str='\Functions';
addpath(strcat(currloc,str));

% For atlas labels
str='\FinalAnalysis\Inputs';
addpath(strcat(currloc,str));

%%
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORGANISING RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I would like three sets of figs - one for Sham vs rsfMRI, TI vs rsfMRI,
and TIst vs rsfMRI
Order of the Res Table:
Cols 1-2 = ROI nums
Cols 3-5 = p, t, cohens D for pre vs post
Cols 6-8 = p, t, cohens D for post vs TI11
Cols 9-11 = 3-1: p, t, cohens D for post vs TI13
%}

load('ComputeMI_5.mat')
cd(outputpath)
%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Make circle plots for Pre vs Post
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Get degree per network
[WorkingTable,PropRedNetList,RedNetNames] = MakeFigs_FxnToComputeDegreePerNetwork(ResTbl(:,1:2));

RedNumNetworks=9;
A=zeros(RedNumNetworks);
for ii=1:size(WorkingTable,1)
for aa=1:RedNumNetworks
for bb=1:RedNumNetworks

    T=isequal([WorkingTable(ii,5),WorkingTable(ii,6)],[aa,bb]);
    if  T == logical(1)
        A(WorkingTable(ii,5),WorkingTable(ii,6))= A(WorkingTable(ii,5),WorkingTable(ii,6)) + 1; 
    end    
end
end
end

C=linspecer(2);
RedNetNames2={"Vis","SoMat","DorsAttn","SalVentAttn","Limbic","Control","DMN","TempPar","Subcor"}';
CC(1:9,:)=[C(1,:);C(1,:);C(1,:);C(1,:);C(1,:);C(1,:);C(1,:);C(1,:);C(1,:)];
figure()
circularGraph(A,'Label',RedNetNames2,'Colormap',CC)
set(gca,'FontSize',20)
set(gca,'FontName','Arial')
saveas(gcf,'CirclePlot_PreVsPost.png');
close all

%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Make circle plots for post vs TI 1:1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Only limit results to those with sig diff miFC between Post vs TI11
tmpResTbl=ResTbl;
toDel=tmpResTbl(:,6)>0.05;
tmpResTbl(toDel,:)=[];

[WorkingTable,PropRedNetList,RedNetNames] = MakeFigs_FxnToComputeDegreePerNetwork(tmpResTbl(:,1:2));

A=zeros(RedNumNetworks);
for ii=1:size(WorkingTable,1)
for aa=1:RedNumNetworks
for bb=1:RedNumNetworks

    T=isequal([WorkingTable(ii,5),WorkingTable(ii,6)],[aa,bb]);
    if  T == logical(1)
        A(WorkingTable(ii,5),WorkingTable(ii,6))= A(WorkingTable(ii,5),WorkingTable(ii,6)) + 1; 
    end    
end
end
end

CC(1:9,:)=[C(2,:);C(2,:);C(2,:);C(2,:);C(2,:);C(2,:);C(2,:);C(2,:);C(2,:)];
figure()
circularGraph(A,'Label',RedNetNames2,'Colormap',CC)
set(gca,'FontSize',20)
set(gca,'FontName','Arial')
saveas(gcf,'CirclePlot_PostVsTI11.png');
close all

%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Make circle plots for post vs TI 1:3
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Only limit results to those with sig diff miFC between Post vs TI11
tmpResTbl=ResTbl;
toDel=tmpResTbl(:,9)>0.05;
tmpResTbl(toDel,:)=[];

[WorkingTable,PropRedNetList,RedNetNames] = MakeFigs_FxnToComputeDegreePerNetwork(tmpResTbl(:,1:2));

A=zeros(RedNumNetworks);
for ii=1:size(WorkingTable,1)
for aa=1:RedNumNetworks
for bb=1:RedNumNetworks

    T=isequal([WorkingTable(ii,5),WorkingTable(ii,6)],[aa,bb]);
    if  T == logical(1)
        A(WorkingTable(ii,5),WorkingTable(ii,6))= A(WorkingTable(ii,5),WorkingTable(ii,6)) + 1; 
    end    
end
end
end

figure()
circularGraph(A,'Label',RedNetNames2,'Colormap',CC)
set(gca,'FontSize',20)
set(gca,'FontName','Arial')
saveas(gcf,'CirclePlot_PostVsTI13.png');
close all

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make hippocampal violins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Taking results from each portion of the hippocampus
% 204-206 = LH ant, mid, post
% 213-215 = RH ant, mid, post
clear AntHipRes MidHipRes PostHipRes

% Get degree per network
[WorkingTable,PropRedNetList,RedNetNames] = MakeFigs_FxnToComputeDegreePerNetwork(ResTbl(:,1:2));

a=0;
m=0;
p=0;

% Gives ROI1, ROI2, p-val, and z-stat
for ii=1:size(ResTbl,1)
    if ResTbl(ii,1)==204 || ResTbl(ii,1)==213
        a=a+1;
        AntHipRes(a,1:4)=ResTbl(ii,1:4);
        AntHipResName(a,1)=RedNetNames(WorkingTable(ii,6),1);

    elseif ResTbl(ii,1)==205 || ResTbl(ii,1)==214
        m=m+1;
        MidHipRes(m,1:4)=ResTbl(ii,1:4);
        MidHipResName(m,1)=RedNetNames(WorkingTable(ii,6));

    elseif ResTbl(ii,1)==206 || ResTbl(ii,1)==215
        p=p+1;
        PostHipRes(p,1:4)=ResTbl(ii,1:4);
        PostHipResName(p,1)=RedNetNames(WorkingTable(ii,6));
    else
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find ROIs that are a part of the AT or PM networks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ATNetROIs=[51,52,53,54,83,92,93,94,98,144,157,158,161,162,163,164,177,193,196,202,211];
PMNetROIs=[44,68,72,73,74,75,77,78,79,87,97,98,100,108,110,134,139,148,149,164,165,173,181,182,183,184,186,195,196];

for ii=1:size(AntHipRes,1)
    if ismember(AntHipRes(ii,2),PMNetROIs)
        AntHipResName(ii,2)="PM";
        AntHipResName(ii,3)=table2array(SchaeferNames(AntHipRes(ii,2),1));
        AntHipResName(ii,4)=AntHipRes(ii,1);
        AntHipResName(ii,5)=AntHipRes(ii,2);

    elseif ismember(AntHipRes(ii,2),ATNetROIs)
        AntHipResName(ii,2)="AT";
        AntHipResName(ii,3)=table2array(SchaeferNames(AntHipRes(ii,2),1));
        AntHipResName(ii,4)=AntHipRes(ii,1);
        AntHipResName(ii,5)=AntHipRes(ii,2);        
    else
    end
end

for ii=1:size(MidHipRes,1)
    if ismember(MidHipRes(ii,2),PMNetROIs)
        MidHipResName(ii,2)="PM";
        MidHipResName(ii,3)=table2array(SchaeferNames(MidHipRes(ii,2),1));
        MidHipResName(ii,4)=MidHipRes(ii,1);
        MidHipResName(ii,5)=MidHipRes(ii,2);

    elseif ismember(MidHipRes(ii,2),ATNetROIs)
        MidHipResName(ii,2)="AT";
        MidHipResName(ii,3)=table2array(SchaeferNames(MidHipRes(ii,2),1));
        MidHipResName(ii,4)=MidHipRes(ii,1);
        MidHipResName(ii,5)=MidHipRes(ii,2);
    else
    end
end

for ii=1:size(PostHipRes,1)
    if ismember(PostHipRes(ii,2),PMNetROIs)
        PostHipResName(ii,2)="PM";
        PostHipResName(ii,3)=table2array(SchaeferNames(PostHipRes(ii,2),1));
        PostHipResName(ii,4)=PostHipRes(ii,1);
        PostHipResName(ii,5)=PostHipRes(ii,2);

    elseif ismember(PostHipRes(ii,2),ATNetROIs)
        PostHipResName(ii,2)="AT";
        PostHipResName(ii,3)=table2array(SchaeferNames(PostHipRes(ii,2),1));
        PostHipResName(ii,4)=PostHipRes(ii,1);
        PostHipResName(ii,5)=PostHipRes(ii,2);        
    else
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consolidate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AntHipResName(:,6:9)=AntHipRes;
AntHipResName=rmmissing(AntHipResName);
MidHipResName(:,6:9)=MidHipRes;
MidHipResName=rmmissing(MidHipResName);
PostHipResName(:,6:9)=PostHipRes;
PostHipResName=rmmissing(PostHipResName);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW DISTRIB PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gonna try to make 6 plots - one for AT & PM ROIs w/ sig diff miFC to Ant,
% Mid, and Post Hip.
% Gonna start with ant hip to AT and PM ROIs
% Using AntHipResName, I can see rows 11, 20, 24 and 27 are AT ROIs; the rest are
% PM
C=linspecer(4);
AntHipResNameAT=AntHipResName([11,20,24,27],:);
AntHipResNamePM=AntHipResName([1:10,12:19,21:26,28:end],:);

% Let's get all miFC values for the AT ROIs to AntHip
for zz=1:size(AntHipResNameAT,1)
for subj=1:n_subj
for sesh=1:n_sess
ROIMI(subj,sesh)=MI{subj,sesh}(str2num(AntHipResNameAT(zz,4)),str2num(AntHipResNameAT(zz,5)));
end
end
if zz==1
AntHipResAT=ROIMI;
else
AntHipResAT=[AntHipResAT;ROIMI];
end
end

figure()
hold on
set(gca,'FontSize',20)
set(gca,'FontName','Arial')
violins = violinplot(AntHipResAT);
for ii=1:4
violins(1,ii).ShowData = 1;
violins(1,ii).ShowMean = 'yes';
violins(1,ii).ViolinColor = C(ii,:);
end
hold off
saveas(gcf,'AntHip_AT.png');
close all

% Trying PM ROIs now
for zz=1:size(AntHipResNamePM,1)
for subj=1:n_subj
for sesh=1:n_sess
ROIMI(subj,sesh)=MI{subj,sesh}(str2num(AntHipResNamePM(zz,4)),str2num(AntHipResNamePM(zz,5)));
end
end
if zz==1
AntHipResPM=ROIMI;
else
AntHipResPM=[AntHipResPM;ROIMI];
end
end

figure()
hold on
set(gca,'FontSize',20)
set(gca,'FontName','Arial')
violins = violinplot(AntHipResPM);
for ii=1:4
violins(1,ii).ShowData = 1;
violins(1,ii).ShowMean = 'yes';
violins(1,ii).ViolinColor = C(ii,:);
end
hold off
saveas(gcf,'AntHip_PM.png');
close all


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mid Hip plot time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MidHipResNameAT=MidHipResName([1,3,8,9],:);
MidHipResNamePM=MidHipResName([2,4:7,10:end],:);

% Let's get all miFC values for the AT ROIs to AntHip
for zz=1:size(MidHipResNameAT,1)
for subj=1:n_subj
for sesh=1:n_sess
ROIMI(subj,sesh)=MI{subj,sesh}(str2num(MidHipResNameAT(zz,4)),str2num(MidHipResNameAT(zz,5)));
end
end
if zz==1
MidHipResAT=ROIMI;
else
MidHipResAT=[MidHipResAT;ROIMI];
end
end

figure()
hold on
set(gca,'FontSize',20)
set(gca,'FontName','Arial')
violins = violinplot(MidHipResAT);
for ii=1:4
violins(1,ii).ShowData = 1;
violins(1,ii).ShowMean = 'yes';
violins(1,ii).ViolinColor = C(ii,:);
end
hold off
saveas(gcf,'MidHip_AT.png');
close all

% Trying PM ROIs now
for zz=1:size(MidHipResNamePM,1)
for subj=1:n_subj
for sesh=1:n_sess
ROIMI(subj,sesh)=MI{subj,sesh}(str2num(MidHipResNamePM(zz,4)),str2num(MidHipResNamePM(zz,5)));
end
end
if zz==1
MidHipResPM=ROIMI;
else
MidHipResPM=[MidHipResPM;ROIMI];
end
end

figure()
hold on
set(gca,'FontSize',20)
set(gca,'FontName','Arial')
violins = violinplot(MidHipResPM);
for ii=1:4
violins(1,ii).ShowData = 1;
violins(1,ii).ShowMean = 'yes';
violins(1,ii).ViolinColor = C(ii,:);
end
hold off
saveas(gcf,'MidHip_PM.png');
close all


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post Hip plot time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PostHipResNameAT=PostHipResName([2,4,5],:);
PostHipResNamePM=PostHipResName([1,3],:);

% Let's get all miFC values for the AT ROIs to AntHip
for zz=1:size(PostHipResNameAT,1)
for subj=1:n_subj
for sesh=1:n_sess
ROIMI(subj,sesh)=MI{subj,sesh}(str2num(PostHipResNameAT(zz,4)),str2num(PostHipResNameAT(zz,5)));
end
end
if zz==1
PostHipResAT=ROIMI;
else
PostHipResAT=[PostHipResAT;ROIMI];
end
end

figure()
hold on
set(gca,'FontSize',20)
set(gca,'FontName','Arial')
violins = violinplot(PostHipResAT);
for ii=1:4
violins(1,ii).ShowData = 1;
violins(1,ii).ShowMean = 'yes';
violins(1,ii).ViolinColor = C(ii,:);
end
hold off
saveas(gcf,'PostHip_AT.png');
close all

% Trying PM ROIs now
for zz=1:size(PostHipResNamePM,1)
for subj=1:n_subj
for sesh=1:n_sess
ROIMI(subj,sesh)=MI{subj,sesh}(str2num(PostHipResNamePM(zz,4)),str2num(PostHipResNamePM(zz,5)));
end
end
if zz==1
PostHipResPM=ROIMI;
else
PostHipResPM=[PostHipResPM;ROIMI];
end
end

figure()
hold on
set(gca,'FontSize',20)
set(gca,'FontName','Arial')
violins = violinplot(PostHipResPM);
for ii=1:4
violins(1,ii).ShowData = 1;
violins(1,ii).ShowMean = 'yes';
violins(1,ii).ViolinColor = C(ii,:);
end
hold off
saveas(gcf,'PostHip_PM.png');
close all

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See what AT and PM network ROIs we need to keep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ATNetROIs=[51:54,83,92:94,98,144,157,158,161:164,177,193,196,202,211];
PMNetROIs=[44,68,72,73,74,75,77,78,79,87,97,98,100,108,110,134,139,148,149,164,165,173,181,182,183,184,186,195,196];

PostTaskOrTIDecrease=ResTbl(ResTbl(:,4)<0,:);

% How many were subcortical?
HippROIs=[204,205,206,213,214,215];
h=0;

for ii=1:size(PostTaskOrTIDecrease,1)
% If subcortical, were they hippocampal?
if ismember(PostTaskOrTIDecrease(ii,1),HippROIs)
h=h+1;
PostTaskOrTIDecreaseHipp(h,:)=PostTaskOrTIDecrease(ii,:);
end
end


h1=0;
h2=0;
for ii=1:size(PostTaskOrTIDecreaseHipp,1)
if ismember(PostTaskOrTIDecreaseHipp(ii,2),ATNetROIs)
h1=h1+1;
ATNetRegion(h1,1)=PostTaskOrTIDecreaseHipp(ii,2);

elseif ismember(PostTaskOrTIDecreaseHipp(ii,2),PMNetROIs)
h2=h2+1;
PMNetRegion(h2,1)=PostTaskOrTIDecreaseHipp(ii,2);

else
end
end
