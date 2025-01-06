
load('C:\Users\dk818\OneDrive - Imperial College London\INS\TIrsfMRI_BrainStateTopologyAndDynamics\FinalAnalysis\Inputs\SegHip_UseToComputeDegPerNet.mat')

Cont=[ContA;ContB;ContC];
DMN=[DMNA;DMNB;DMNC];
DorsAttn=[DorsAttnA;DorsAttnB];
Lim=[LimA;LimB];
SalVentAttn=[SalVentAttnA;SalVentAttnB];
SoMat=[SoMatA;SoMatB];
SubCor=SubCor;
TmpPar=TmpPar;
Vis=[VisCent;VisPeri];

ContCount=0;
DMNCount=0;
DorsAttnCount=0;
LimCount=0;
SalVentAttnCount=0;
SoMatCount=0;
SubCorCount=0;
TmpParCount=0;
VisCount=0;
ATCount=0;
PMCount=0;

MassiveNetworkTbl=[Cont;DMN;DorsAttn;Lim;SalVentAttn;SoMat;SubCor;TmpPar;Vis];
MassiveNetworkTbl(:,2)=[ones(size(Cont,1),1)*1;ones(size(DMN,1),1)*2;ones(size(DorsAttn,1),1)*3;ones(size(Lim,1),1)*4;ones(size(SalVentAttn,1),1)*5;ones(size(SoMat,1),1)*6;ones(size(SubCor,1),1)*7;ones(size(TmpPar,1),1)*8;ones(size(Vis,1),1)*9];

% Now add AT or PM ROIs as networks 10 and 11, respectively
ATNetROIs=[51:54,83,92:94,98,144,157,158,161:164,177,193,196,202,211];
PMNetROIs=[44,68,72,73,74,75,77,78,79,87,97,98,100,108,110,134,139,148,149,164,165,173,181,182,183,184,186,195,196];

for ii=1:size(MassiveNetworkTbl,1)
if ismember(MassiveNetworkTbl(ii,1),ATNetROIs)
    MassiveNetworkTbl(ii,2)=10;
elseif ismember(MassiveNetworkTbl(ii,1),PMNetROIs)
    MassiveNetworkTbl(ii,2)=11;
else
end
end

clearvars -except MassiveNetworkTbl
saveas('MegaNetworkTable.mat')



