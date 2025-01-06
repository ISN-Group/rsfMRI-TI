% Add paths
currloc='C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics';
% For outputs
outputpath='C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\FinalAnalysis\Outputs';
addpath(outputpath);
% For Brain Connectivity Toolbox
addpath('C:\Users\dk818\OneDrive - Imperial College London\LilThingsAndMATLABPrograms\BCT-main\BCT\2019_03_03_BCT')
% For functions
addpath('C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\Functions')
load('ComputeMI_5.mat')
clearvars -except MI n_sess n_subj 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Global Connectivity, Modularity, Route efficiency, adn
% Betweenness Centrality
%%%%%%%%%%%%%%%%%%%%%%%%%
for subj=1:n_subj
for sesh=1:n_sess
    GlobalConn(subj,sesh)=mean(MI{subj,sesh},'all');
    [~,Q(subj,sesh)]=modularity_und(MI{subj,sesh});
    ErOut(subj,sesh)=rout_efficiency(MI{subj,sesh});
end
end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for effect of session on modularity, global conn, and route eff
%%%%%%%%%%%%%%%%%%%%%%%%%
% MI cols 1-4 = TI11, TI13, Post, Pre
% Gonna make rows 1:3 PrevsPost, PostvsTI11, and PostvsTI13
% Cols = p, tstat, cohensD

% PrevsPost
[ResQ(1,1),ResQ(1,2),ResQ(1,3),~]=PermTest(Q(:,3),Q(:,4));
[ResGlobalConn(1,1),ResGlobalConn(1,2),ResGlobalConn(1,3),~]=PermTest(GlobalConn(:,3),GlobalConn(:,4));
[ResErOut(1,1),ResErOut(1,2),ResErOut(1,3),~]=PermTest(ErOut(:,3),ErOut(:,4));

% PostvsTI11 
[ResQ(2,1),ResQ(2,2),ResQ(2,3),~]=PermTest(Q(:,3),Q(:,1));
[ResGlobalConn(2,1),ResGlobalConn(2,2),ResGlobalConn(2,3),~]=PermTest(GlobalConn(:,3),GlobalConn(:,1));
[ResErOut(2,1),ResErOut(2,2),ResErOut(2,3),~]=PermTest(ErOut(:,3),ErOut(:,1));

% PostvsTI13
[ResQ(3,1),ResQ(3,2),ResQ(3,3),~]=PermTest(Q(:,3),Q(:,2));
[ResGlobalConn(3,1),ResGlobalConn(3,2),ResGlobalConn(3,3),~]=PermTest(GlobalConn(:,3),GlobalConn(:,2));
[ResErOut(3,1),ResErOut(3,2),ResErOut(3,3),~]=PermTest(ErOut(:,3),ErOut(:,2));
