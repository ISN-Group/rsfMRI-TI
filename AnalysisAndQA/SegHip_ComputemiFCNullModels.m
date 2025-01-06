%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add paths
addpath('C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\FinalAnalysis');
% For outputs
outputpath='C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\FinalAnalysis\Outputs';
addpath(outputpath);
% For atlas and other data
addpath('C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\AtlasAndOtherInputs');
% For data
addpath('C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\Timeseries');
homeDir='C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics';
% For miFC
addpath('C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\Functions');
addpath('C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\Functions\mi');
load('ComputeMI_5.mat')
clearvars -except df ResTbl outputpath
[n_subj,n_sess]=size(df);
[Tmax,num_rois]=size(df{1,1});
TR=2;
fnq=1/(2*TR);                 % Nyquist frequency.
flp = .02;                    % lowpass frequency of filter (Hz). This allows anything below 50 Hz. 
fhi = 0.1;                    % highpass- we've already applied one, though it is less inclusive (allowing anything over 0.01 to pass). I will keep this for now. 
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency. 
k=2;                          % 2nd order butterworth filter. This determines the steepness of the gain function (and 2 is pretty smooth). 
[bfilt,afilt]=butter(k,Wn);   % "Butter" is a MATLAB function than constructs the butterworth filter using defined cutoff frequencies.
clear df Wn fhi flp fnq k

%%
for zz=1:1000
disp(strcat('Currently on permutation ',num2str(zz)))
% Start randomisation of timeseries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
load('df.mat')
dfRand={};
for s=1:n_subj 
    for task=1:n_sess 
        
        % Get the BOLD signals from this subject in this task
        BOLD = df{s,task};
        Phase_BOLD=zeros(Tmax,num_rois); 
        Rand_BOLD=zeros(Tmax,num_rois);
        ShuffSig=zeros(Tmax,num_rois);
        % Get the BOLD phase using the Hilbert transform
        for seed=1:num_rois
            BOLD(:,seed)=BOLD(:,seed)-mean(BOLD(:,seed)); % demean the timecourse
            signal_filt =filtfilt(bfilt,afilt,BOLD(:,seed)); % filter                   
            Phase_BOLD(:,seed) = angle(fft(signal_filt)); % extract phase angle
        end

        % Now shuffle all phase coefficients per region
        % Preserve the row indices
        rowIndex = repmat((1:Tmax)',[1 num_rois]);
        % Get randomized column indices by sorting a second random array
        [~,randomizedColIndex] = sort(rand(Tmax,num_rois),2);
        % Need to use linear indexing to create B
        newLinearIndex = sub2ind([Tmax,num_rois],rowIndex,randomizedColIndex);
        ShuffSig = Phase_BOLD(newLinearIndex);
        
        for seed=1:num_rois
            % Fourier transform of original signal
            z =fft(filtfilt(bfilt,afilt,BOLD(:,seed))); 

            % Use the shuffled phase angles to get the permuted signal
            % The angles in theta are such that z = abs(z).*exp(i*theta).
            ZZ=abs(z).*exp(i*ShuffSig(:,seed));

            % Do the inverse transform to get the real, but permuted, signal
            Rand_BOLD(:,seed)=ifft(ZZ,'symmetric');
            
        end
        dfRand{s,task}=Rand_BOLD;
    end
end
df=dfRand;

% Zscore everything to have centre mean 0 and 1 SD
BigStructNorm={};
for subj=1:n_subj
for sess=1:n_sess
for rr=1:num_rois
    BigStructNorm{subj,sess}(:,rr)=zscore(df{subj,sess}(:,rr));
end
end
end

%
% Compute mi for edges w/ sig diff miFC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=0;
for ii=1:size(ResTbl,1)
% Get the first pair of ROIs
rr1=ResTbl(ii,1);
rr2=ResTbl(ii,2);

% Compute MI
for subj=1:n_subj
for sesh=1:n_sess
tmp=BigStructNorm{subj,sesh};
ROIMI(subj,sesh)= mutualinfo(tmp(:,rr1),tmp(:,rr2));
end
end

% Compute pre vs post diff
h=h+1;
[PreVsPost_P(rr1,rr2,zz),PreVsPost_Tstat(rr1,rr2,zz),PreVsPost_CohensD(rr1,rr2,zz),PreVsPost_MaxT(h,zz)]=PermTest(ROIMI(:,3),ROIMI(:,4));
end

% Correct for family wise error
% See if Tstat from sig p vals are outside
% the 95% CI for the MaxT distribution 
% Find the 95% CI for the MaxT
maxTDist=fitdist(PreVsPost_MaxT(:,zz), 'Normal');
CI=paramci(maxTDist);
h=0;
for ii=1:size(ResTbl,1)
rr1=ResTbl(ii,1);
rr2=ResTbl(ii,2);
if rr1 ~= rr2 && PreVsPost_P(rr1,rr2,zz)<0.05
h=h+1;
if PreVsPost_Tstat(rr1,rr2,zz)<CI(1,1) || PreVsPost_Tstat(rr1,rr2,zz)>CI(2,1)
PreVsPostLongFormat(h,1)=rr1;
PreVsPostLongFormat(h,2)=rr2;
PreVsPostLongFormat(h,3)=PreVsPost_P(rr1,rr2,zz);
end
end
end

% If PreVsPost is signficiant, run PostVsTI and PostvsTIst
h=0;
for ii=1:size(PreVsPostLongFormat)
rr1=PreVsPostLongFormat(ii,1);
rr2=PreVsPostLongFormat(ii,2);

for subj=1:n_subj
for sesh=1:n_sess
tmp=BigStructNorm{subj,sesh};
ROIMI(subj,sesh)=mutualinfo(tmp(:,rr1),tmp(:,rr2));
end
end

% Run and store pairwise tests
[PostVsTI11_P(rr1,rr2,zz),PostVsTI11_Tstat(rr1,rr2,zz),PostVsTI11_CohensD(rr1,rr2,zz),~]=PermTest(ROIMI(:,3),ROIMI(:,1));
[PostVsTI13_P(rr1,rr2,zz),PostVsTI13_Tstat(rr1,rr2,zz),PostVsTI13_CohensD(rr1,rr2,zz),~]=PermTest(ROIMI(:,3),ROIMI(:,2));
end

end

%%
% NOTE - there's no way to run this for 10,000 permutations, since it takes
% ~1 min per permutation, and 10,000 perms = 166 hrs. 1,000 perms = 16 hrs.
% I will run this off and on, and continue adding to a database of the pval
% distributions. 
cd(outputpath)

% Clean up
clearvars -except zz PreVsPost* PostVs* outputpath
if size(PreVsPost_P,3)==zz
PreVsPost_P(:,:,zz)=[];
PreVsPost_Tstat(:,:,zz)=[];
PreVsPost_CohensD(:,:,zz)=[];
end

if size(PostVsTI11_P,3)==zz
PostVsTI11_P(:,:,zz)=[];
PostVsTI11_Tstat(:,:,zz)=[];
PostVsTI11_CohensD(:,:,zz)=[];
end

if size(PostVsTI13_P,3)==zz
PostVsTI13_P(:,:,zz)=[];
PostVsTI13_Tstat(:,:,zz)=[];
PostVsTI13_CohensD(:,:,zz)=[];
end



% Load and add to it
if isfile('NullPvalDistrib.mat')

load('NullPvalDistrib.mat')
start=size(NullPreVsPost_P,3)+1;
End=start+(zz-2);
NullPreVsPost_P(:,:,start:End)=PreVsPost_P(:,:,:);
NullPreVsPost_Tstat(:,:,start:End)=PreVsPost_Tstat(:,:,:);
NullPreVsPost_CohensD(:,:,start:End)=PreVsPost_CohensD(:,:,:);

start=size(NullPostVsTI11_P,3)+1;
End=start+(zz-2);
NullPostVsTI11_P(:,:,start:End)=PostVsTI11_P(:,:,:);
NullPostVsTI11_Tstat(:,:,start:End)=PostVsTI11_Tstat(:,:,:);
NullPostVsTI11_CohensD(:,:,start:End)=PostVsTI11_CohensD(:,:,:);

start=size(NullPostVsTI13_P,3)+1;
End=start+(zz-2);
NullPostVsTI13_P(:,:,start:End)=PostVsTI13_P(:,:,:);
NullPostVsTI13_Tstat(:,:,start:End)=PostVsTI13_Tstat(:,:,:);
NullPostVsTI13_CohensD(:,:,start:End)=PostVsTI13_CohensD(:,:,:);

clearvars -except Null*
save('NullPvalDistrib.mat','-append')

else
NullPreVsPost_P=PreVsPost_P;
NullPreVsPost_CohensD=PreVsPost_CohensD;
NullPreVsPost_Tstat=PreVsPost_Tstat;

NullPostVsTI11_P=PostVsTI11_P;
NullPostVsTI11_Tstat=PostVsTI11_Tstat;
NullPostVsTI11_CohensD=PostVsTI11_CohensD;

NullPostVsTI13_P=PostVsTI13_P;
NullPostVsTI13_Tstat=PostVsTI13_Tstat;
NullPostVsTI13_CohensD=PostVsTI13_CohensD;

save('NullPvalDistrib.mat','append')
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute significance after permutation testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The objective of permutation testing is to see how many more
% values are more extreme than our empirical value, or test
% statistic. So our test statistic is the Pvalue for each pair of
% regions, and we need to see what proportion of permutated pvalues
% are SMALLER than our pvalue. That % is our alpha, or pvalue. 

load('ComputeMI_5.mat') 
clearvars -except ResTbl
load('NullPvalDistrib.mat')
sigCount=size(ResTbl,1);
SurviveNullSigCount=0;

% For each Pvalue we've deemed signficant. . . 
for ii=1:size(ResTbl,1)
rr1=ResTbl(ii,1);
rr2=ResTbl(ii,2);

% . . .get the distribution of pvals from permutation testing. . . 
arr=squeeze(NullPreVsPost_P(rr1,rr2,:));

% find what percent of the total number of permutation pvalues are more extreme
% (i.e., smaller) than the test statistic (i.e., our empirical
% pvalue)
ResTblSurvival(ii,1)=sum(ResTbl(ii,4)>arr)/(size(arr,1)+1);

% If that % is smaller than 0.05, then this result survives
% permutation testing.
if ResTblSurvival(ii,1)<0.05
    SurviveNullSigCount=SurviveNullSigCount+1;
end

end

PropSaved=(SurviveNullSigCount/sigCount)*100;
disp(strcat(num2str(PropSaved),'% of results survive null testing'))


%%
% See what AT and PM regions are saved! 
ResTblSurvived=ResTbl(ResTblSurvival(:,1)<0.05,:);

