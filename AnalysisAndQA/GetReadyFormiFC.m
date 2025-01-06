%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add paths
addpath('C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\FinalAnalysis');
% For outputs
outputpath='C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\FinalAnalysis\Outputs';
addpath(outputpath);
% For data
TimeseriesDir='C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\Timeseries';
addpath(TimeseriesDir);
df={};
num_rois=218;
CondName={'TI_','TIst','SHAM','Pre'};
n_subj=20;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and organise data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:size(CondName,2)
disp(strcat('Loading files for Cond #',num2str(ii)))
Files_Cor=dir(strcat(TimeseriesDir,'\*',CondName{1,ii},'*Cortex*'));
Files_l_acc=dir(strcat(TimeseriesDir,'\*',CondName{1,ii},'*l_acc*'));
Files_l_amy=dir(strcat(TimeseriesDir,'\*',CondName{1,ii},'*l_amy*'));
Files_l_cau=dir(strcat(TimeseriesDir,'\*',CondName{1,ii},'*l_cau*'));
Files_l_anthip=dir(strcat(TimeseriesDir,'\*',CondName{1,ii},'*Ant_LH*'));
Files_l_midhip=dir(strcat(TimeseriesDir,'\*',CondName{1,ii},'*Mid_LH*'));
Files_l_posthip=dir(strcat(TimeseriesDir,'\*',CondName{1,ii},'*Post_LH*'));
Files_l_pal=dir(strcat(TimeseriesDir,'\*',CondName{1,ii},'*l_pal*'));
Files_l_put=dir(strcat(TimeseriesDir,'\*',CondName{1,ii},'*l_put*'));
Files_l_thl=dir(strcat(TimeseriesDir,'\*',CondName{1,ii},'*l_thl*'));
Files_r_nac=dir(strcat(TimeseriesDir,'\*',CondName{1,ii},'*r_nac*'));
Files_r_amy=dir(strcat(TimeseriesDir,'\*',CondName{1,ii},'*r_amy*'));
Files_r_cau=dir(strcat(TimeseriesDir,'\*',CondName{1,ii},'*r_cau*'));
Files_r_anthip=dir(strcat(TimeseriesDir,'\*',CondName{1,ii},'*Ant_RH*'));
Files_r_midhip=dir(strcat(TimeseriesDir,'\*',CondName{1,ii},'*Mid_RH*'));
Files_r_posthip=dir(strcat(TimeseriesDir,'\*',CondName{1,ii},'*Post_RH*'));
Files_r_pal=dir(strcat(TimeseriesDir,'\*',CondName{1,ii},'*r_pal*'));
Files_r_put=dir(strcat(TimeseriesDir,'\*',CondName{1,ii},'*r_put*'));
Files_r_thl=dir(strcat(TimeseriesDir,'\*',CondName{1,ii},'*r_thl*'));

for subj=1:n_subj
% Load and transpose cortex 
Cor = load(Files_Cor(subj).name);
 
% Now the l_acc
l_acc=load(Files_l_acc(subj).name);

% Now the l_amy
l_amy=load(Files_l_amy(subj).name);

% Now the l_cau
l_cau=load(Files_l_cau(subj).name);

% Now the ant, med, and post l_hip
l_anthip=load(Files_l_anthip(subj).name);
l_midhip=load(Files_l_midhip(subj).name);
l_posthip=load(Files_l_posthip(subj).name);

% Now the l_pal
l_pal=load(Files_l_pal(subj).name);

% Now the l_put
l_put=load(Files_l_put(subj).name);

% Now the l_thl
l_thl=load(Files_l_thl(subj).name);

% Now the r_nac
r_nac=load(Files_r_nac(subj).name);

% Now the r_amy
r_amy=load(Files_r_amy(subj).name);

% Now the r_cau
r_cau=load(Files_r_cau(subj).name);

% Now the r_hip
r_anthip=load(Files_r_anthip(subj).name);
r_midhip=load(Files_r_midhip(subj).name);
r_posthip=load(Files_r_posthip(subj).name);

% Now the r_pal
r_pal=load(Files_r_pal(subj).name);

% Now the r_put
r_put=load(Files_r_put(subj).name);

% Now the r_thl
r_thl=load(Files_r_thl(subj).name);

% % Stack 'em
stack = [Cor,l_acc,l_amy,l_cau,l_anthip,l_midhip,l_posthip,l_pal,l_put,l_thl,r_nac,r_amy,r_cau,r_anthip,r_midhip,r_posthip,r_pal,r_put,r_thl] ;
df{subj,ii}=stack;

end
end

%%
cd(outputpath)
clearvars -except df
save('df.mat')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make shuffled data for null models
dfRand={};
[n_subj,n_sess]=size(df);
[Tmax,num_rois]=size(df{1,1});
TR=2;
