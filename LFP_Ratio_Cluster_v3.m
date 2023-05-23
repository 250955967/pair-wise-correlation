function LFP_Ratio_Cluster_v3(stExpFN,nCh,stOption_Sp,strSavePath)
%stExpFN.Pre.LogFN
%stExpFN.Pre.NevFN
%stExpFN.Cond.LogFN
%stExpFN.Cond.NevFN
%stExpFN.Post.LogFN
%stExpFN.Post.NevFN
%stExpFN.CCMatFN
%stExpFN.RTableFN

switch nargin
    case 0
       error('SG:LFP_RATIO_CLUSTER:NotEnoughInputs','Not enough input arguments!');

    case 1
        nCh = 1;
        stOption_Sp.tWnd = 2;
        stOption_Sp.tShift = 0.2;
        stOption_Sp.Freq = 0:0.1:60;
        strSavePath = [];

    case 2
        stOption_Sp.tWnd = 2;
        stOption_Sp.tShift = 0.2;
        stOption_Sp.Freq = 0:0.1:60;
        strSavePath = [];

    case 3
        strSavePath = [];

    otherwise
end

if isempty(stOption_Sp)
    stOption_Sp.tWnd = 2;
    stOption_Sp.tShift = 0.2;
    stOption_Sp.Freq = 0:0.1:60;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nsResult Pre_nsData Pre_logData] = GetExpData(stExpFN.Pre,'Induce',0,true,'all');
[nsResult Cond_nsData Cond_logData] = GetExpData(stExpFN.Cond,'Line Real',0,true,'all');
if(~isequalwithequalnans(Pre_logData.NeuroProp,Cond_logData.NeuroProp))
    error('SG:LFP_RATIO_CLUSTER:NeuroPropChanged','The Cond NeuroProp is changed!');
end
[nsResult Post_nsData Post_logData] = GetExpData(stExpFN.Post,'Induce',0,true,'all');
if(~isequalwithequalnans(Pre_logData.NeuroProp,Post_logData.NeuroProp))
    error('SG:LFP_RATIO_CLUSTER:NeuroPropChanged','The Post NeuroProp is changed!');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(strSavePath)
    strSaveDir = [strSavePath '\' num2str(Pre_nsData.FileInfo.CreateTime.Year) '-' num2str(Pre_nsData.FileInfo.CreateTime.Month) '-' num2str(Pre_nsData.FileInfo.CreateTime.Day)];
    if ~exist(strSaveDir,'dir')
        mkdir(strSaveDir);
    end
    dirlist = dir(strSaveDir);
    nextdirID = length(find([dirlist.isdir]))-1;
    strNextID = sprintf('%02d',nextdirID);
    strSaveDir = [strSaveDir '\' strNextID];
    if ~exist(strSaveDir,'dir')
        mkdir(strSaveDir);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load stimulus time and sequence
nsData =[Pre_nsData,Cond_nsData,Post_nsData];
stLFP(3) = nsData(3).LFPData(nCh);
Marker(3)=struct('S',[],'E',[]);
logData = [Pre_logData,Cond_logData,Post_logData];
Sti_Seq=cell(1,3);
for n=1:3
    if n~=3
        stLFP(n)=nsData(n).LFPData(nCh);
    end
    StartInd = [nsData(n).SpikeData.nChannelID]==130;
    Marker(n).S = nsData(n).SpikeData(StartInd).dTimestamps;
    EndInd = [nsData(n).SpikeData.nChannelID]==131;
    Marker(n).E = nsData(n).SpikeData(EndInd).dTimestamps;
    Sti_Seq{n}=logData(n).StimulusInfo.Test_Sequence;
end
clear('logData','Pre_logData','Cond_logData','Post_logData');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load CCs
stCCMat = load(stExpFN.CCMatFN);
stCCVal(3) = struct('C',[],'A',[]);
CCFields={'PreRho','CondRho','PostRho'};
for nField = 1:3
    stCCVal(nField).C = stCCMat.(CCFields{nField}).C.SVal;
    stCCVal(nField).A = stCCMat.(CCFields{nField}).A.SVal;
end
clear('stCCMat','CCFields');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FieldNames={'Pre','Cond','Post'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate Ratios
Time = struct('Pre',[],'Cond',[],'Post',[]);
PSD= struct('Pre',[],'Cond',[],'Post',[]);
Freq = struct('Pre',[],'Cond',[],'Post',[]);
P_Ratio_1 = struct('Pre',[],'Cond',[],'Post',[]);
P_Ratio_2 = struct('Pre',[],'Cond',[],'Post',[]);

hFig = figure();
for nField = 1:length(FieldNames)
    FieldName = FieldNames{nField};
    Smp_Wnd = 2^nextpow2(stOption_Sp.tWnd *stLFP(nField).dSampleRate);
    Smp_nOverlap = ceil((stOption_Sp.tWnd- stOption_Sp.tShift)*stLFP(nField).dSampleRate);
    [Y,F,T,P]=spectrogram(stLFP(nField).dVal,Smp_Wnd,Smp_nOverlap,stOption_Sp.Freq,stLFP(nField).dSampleRate);
    bF_1_L = F>=1&F<10;
    PSum_1_L = sum(P(bF_1_L,:),1);
    bF_1_H = F>=1&F<25;
    PSum_1_H = sum(P(bF_1_H,:),1);

    bF_2_L = F>=15&F<30;
    PSum_2_L = sum(P(bF_2_L,:),1);
    bF_2_H = F>=15&F<60;
    PSum_2_H = sum(P(bF_2_H,:),1);
    P_Ratio_1.(FieldNames{nField}) = PSum_1_L./PSum_1_H;
    P_Ratio_2.(FieldNames{nField}) = PSum_2_L./PSum_2_H;
    subplot(2,2,nField)
    scatter(P_Ratio_1.(FieldNames{nField}),P_Ratio_2.(FieldNames{nField}),1,'b');
    xlim([0 1]);
    ylim([0 1]);
    title(FieldName);
    PSD.(FieldName) = P;
    Freq.(FieldName) = F;
    Time.(FieldName) = T;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_Ratio_1_All = [P_Ratio_1.Pre P_Ratio_1.Cond P_Ratio_1.Post];
P_Ratio_2_All = [P_Ratio_2.Pre P_Ratio_2.Cond P_Ratio_2.Post];
subplot(2,2,4)
scatter(P_Ratio_1_All,P_Ratio_2_All,1,'b');
xlim([0 1]);
ylim([0 1]);
title('All');

if ~isempty(strSavePath)
    save([strSaveDir '\Spt.mat'],'PSD','Time','Freq','stLFP');
    save([strSaveDir '\P_Ratio.mat'],'P_Ratio_1','P_Ratio_2');
    hgsave(hFig,[strSaveDir '\Ratio_Scatter.fig']);
end
close(hFig);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%smooth
nWinLen = 17;
w=hann(nWinLen);
w=w/sum(w);

hFig = figure('Name','Ratio_Scatter');
hCurve_1 = figure('Name','Curve_1');
hCurve_2 = figure('Name','Curve_2');
nT=floor(nWinLen/2);
for nField = 1:length(FieldNames)
    FieldName = FieldNames{nField};
    tmpH= repmat(P_Ratio_1.(FieldName)(1),[1 nT]);
    tmpR= repmat(P_Ratio_1.(FieldName)(end),[1 nT]);
    P_Ratio_S1.(FieldName) = conv([tmpH P_Ratio_1.(FieldName) tmpR],w);
    P_Ratio_S1.(FieldName)([1:2*nT end-2*nT+1:end])=[];

    tmpH= repmat(P_Ratio_2.(FieldName)(1),[1 nT]);
    tmpR= repmat(P_Ratio_2.(FieldName)(end),[1 nT]);
    P_Ratio_S2.(FieldName) = conv([tmpH P_Ratio_2.(FieldName) tmpR],w);
    P_Ratio_S2.(FieldName)([1:2*nT end-2*nT+1:end])=[];

    figure(hCurve_1);
    subplot(3,1,nField);
    plot(Time.(FieldName),P_Ratio_1.(FieldName),'r');
    hold on;
    plot(Time.(FieldName),P_Ratio_S1.(FieldName),'b');
    xlim([Time.(FieldName)(1) Time.(FieldName)(end)]);
    title(FieldName);

    figure(hCurve_2);
    subplot(3,1,nField);
    plot(Time.(FieldName),P_Ratio_2.(FieldName),'r');
    hold on;
    plot(Time.(FieldName),P_Ratio_S2.(FieldName),'b');
    xlim([Time.(FieldName)(1) Time.(FieldName)(end)]);
    title(FieldName);

    figure(hFig);
    subplot(2,2,nField)
    scatter(P_Ratio_S1.(FieldName),P_Ratio_S2.(FieldName),1,'b');
    xlim([0 1]);
    ylim([0 1]);
    title(FieldName);
end

figure(hFig);
P_Ratio_S1_All = [P_Ratio_S1.Pre P_Ratio_S1.Cond P_Ratio_S1.Post];
P_Ratio_S2_All = [P_Ratio_S2.Pre P_Ratio_S2.Cond P_Ratio_S2.Post];
subplot(2,2,4)
scatter(P_Ratio_S1_All,P_Ratio_S2_All,1,'r');
xlim([0 1]);
ylim([0 1]);
title('All');

if ~isempty(strSavePath)
    save([strSaveDir '\P_Ratio_S.mat'],'P_Ratio_S1','P_Ratio_S2');
    hgsave(hFig,[strSaveDir '\Ratio_Scatter_S.fig']);
    hgsave(hCurve_1,[strSaveDir '\Ratio1_Curve.fig']);
    hgsave(hCurve_2,[strSaveDir '\Ratio2_Curve.fig']);
end
close([hFig hCurve_1 hCurve_2]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clustering
X=[P_Ratio_S1.Pre P_Ratio_S1.Cond P_Ratio_S1.Post;P_Ratio_S2.Pre P_Ratio_S2.Cond P_Ratio_S2.Post]';
[center,U,obj_fcn] = fcm(X, 2,[NaN NaN NaN 0]); %#ok<NASGU>
maxU = max(U);
idxL = find(U(1,:) == maxU);
idxH = find(U(2, :) == maxU);
if center(1,1) < center(2,1)%make [idxL idxH] [low high]
    tmp = idxL;
    idxL=idxH;
    idxH=tmp;
    center = flipud(center);
    clear('tmp');
end
Seg_Len = [ length(P_Ratio_S1.Pre) length(P_Ratio_S1.Cond) length(P_Ratio_S1.Post)];
Seg_Sep_N = [0 cumsum(Seg_Len)];
Seg_idx_L = struct('Pre',[],'Cond',[],'Post',[]);
Seg_idx_H= struct('Pre',[],'Cond',[],'Post',[]);

hFig=figure('Name','Ratio_Scatter_C');
for nField = 1:length(FieldNames)
    FieldName = FieldNames{nField};
    Seg_idx_L.(FieldName)= idxL((Seg_Sep_N(nField)<idxL) &( idxL<=Seg_Sep_N(nField+1)))-Seg_Sep_N(nField);
    Seg_idx_H.(FieldName)= idxH((Seg_Sep_N(nField)<idxH) &( idxH<=Seg_Sep_N(nField+1)))-Seg_Sep_N(nField);

    subplot(2,2,nField);
    scatter(P_Ratio_S1.(FieldName)(Seg_idx_L.(FieldName)),P_Ratio_S2.(FieldName)(Seg_idx_L.(FieldName)),2,'b');
    hold on;
    scatter(P_Ratio_S1.(FieldName)(Seg_idx_H.(FieldName)),P_Ratio_S2.(FieldName)(Seg_idx_H.(FieldName)),2,'r');
    xlim([0 1]);
    ylim([0 1]);
    title(FieldName);
end

subplot(2,2,4);
scatter(X(idxL,1),X(idxL,2),2,'b');
hold on;
scatter(X(idxH,1),X(idxH,2),2,'r');
plot(center(1,1),center(1,2),'ko','markersize',7,'LineWidth',2);
plot(center(2,1),center(2,2),'k*','markersize',7,'LineWidth',2);
xlim([0 1]);
ylim([0 1]);
title('All');

if ~isempty(strSavePath)
    save([strSaveDir '\Seg_idx.mat'],'Seg_idx_L','Seg_idx_H');
    hgsave(hFig,[strSaveDir '\Ratio_Scatter_C.fig']);
end
close(hFig);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the distribution of CCs in different clusters

stCCDist = struct('Pre',[],'Cond',[],'Post',[]);
for nField = 1:length(FieldNames)
    FieldName = FieldNames{nField};
    if nField ~=2
        bC = Sti_Seq{nField}=='C';
        tCS = Marker(nField).S(bC);
        bA = Sti_Seq{nField}=='A';
        tAS = Marker(nField).S(bA);
    else
        tCS = Marker(nField).S(1:2:end);
        tAS=[];
    end

    Time_H =Time.(FieldName)((Seg_idx_H.(FieldName)));
    idxH_Diff=diff(Seg_idx_H.(FieldName));
    idxH_Diff_Ind = find((idxH_Diff ~= 1));
    idxH_S = [0 idxH_Diff_Ind]+1;
    idxH_E = [idxH_Diff_Ind  length(Seg_idx_H.(FieldName))];

    H_C_Ind = zeros(size(tCS));
    nH_C_Count = 0;
    H_A_Ind = zeros(size(tAS));
    nH_A_Count = 0;
    for nSegH = 1:length(idxH_S)
        tSegH_S = Time_H(idxH_S(nSegH));
        tSegH_E = Time_H(idxH_E(nSegH));
        H_C_Ind_tmp=find(tCS>=tSegH_S&tCS<tSegH_E);
        tmpCount= length(H_C_Ind_tmp);
        H_C_Ind(nH_C_Count+1:nH_C_Count+tmpCount)=H_C_Ind_tmp;
        nH_C_Count = nH_C_Count+tmpCount;
        
        H_A_Ind_tmp=find(tAS>=tSegH_S&tAS<tSegH_E);
        tmpCount= length(H_A_Ind_tmp);
        H_A_Ind(nH_A_Count+1:nH_A_Count+tmpCount)=H_A_Ind_tmp;
        nH_A_Count = nH_A_Count+tmpCount;
    end
    if(nH_C_Count<length(H_C_Ind))
        H_C_Ind(nH_C_Count+1:end) = [];
    end
    if(nH_A_Count<length(H_A_Ind))
        H_A_Ind(nH_A_Count+1:end) = [];
    end
    L_C_Ind = setdiff(1:length(tCS),H_C_Ind);
    stCCDist.(FieldName).HInd = H_C_Ind';
    stCCDist.(FieldName).HCCVal = stCCVal(nField).C(H_C_Ind);
    stCCDist.(FieldName).LInd = L_C_Ind;
    stCCDist.(FieldName).LCCVal = stCCVal(nField).C(L_C_Ind);
    
    L_A_Ind = setdiff(1:length(tAS),H_A_Ind);
    stCCDist.(FieldName).HInd_A = H_A_Ind';
    stCCDist.(FieldName).HCCVal_A = stCCVal(nField).A(H_A_Ind);
    stCCDist.(FieldName).LInd_A = L_A_Ind;
    stCCDist.(FieldName).LCCVal_A = stCCVal(nField).A(L_A_Ind);
end

if ~isempty(strSavePath)
    save([strSaveDir '\stCCDist.mat'],'stCCDist');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot clustered LFP and spectrogram
% GapTime = 1;
hFig = figure('Name','LFP_Spectrogram');
hAxes= zeros(1,2);
for nField = 1:length(FieldNames)
    FieldName = FieldNames{nField};
    clf(hFig);
    hAxes(1)=subplot(2,1,1);
    LFP_T=(0:(stLFP(nField).dwItemCount-1))/stLFP(nField).dSampleRate;
    plot(LFP_T,stLFP(nField).dVal,'b');
    hold on;
    Time_H =Time.(FieldName)((Seg_idx_H.(FieldName)));
    idxH_Diff=diff(Seg_idx_H.(FieldName));
    idxH_Diff_Ind = find((idxH_Diff ~= 1));
    idxH_S = [0 idxH_Diff_Ind]+1;
    idxH_E = [idxH_Diff_Ind  length(Seg_idx_H.(FieldName))];
    idxH_Smp_S = round(Time_H(idxH_S)*stLFP(nField).dSampleRate);
    idxH_Smp_E = round(Time_H(idxH_E)*stLFP(nField).dSampleRate);
    for nSeg = 1:length(idxH_S)
         line((idxH_Smp_S(nSeg):idxH_Smp_E(nSeg))/stLFP(nField).dSampleRate,stLFP(nField).dVal(idxH_Smp_S(nSeg):idxH_Smp_E(nSeg)),'color',[1 0 0],'linewidth',1);
    end
    xlim([LFP_T(1) LFP_T(end)]);
    title(FieldName);

    hAxes(2)=subplot(2,1,2);
    imagesc(Time.(FieldName),Freq.(FieldName),10*log10(abs(PSD.(FieldName))));
    set(hAxes(2),'YDir','normal');
    xlabel('Time');
    ylabel('Frequency (Hz)');

    linkaxes(hAxes,'x');
    hzoom = zoom(hFig);
    setAxesZoomMotion(hzoom,hAxes,'horizontal');
    hpan= pan(hFig);
    setAxesPanMotion(hpan,hAxes,'horizontal');
    
    if ~isempty(strSavePath)
        hgsave(hFig,[strSaveDir '\' FieldName '_LFP_Spectrogram.fig']);
    end
end
close(hFig);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hFig= figure();
HistEdge= -1:0.1:1;
for nField = 1:length(FieldNames)
    clf();
    FieldName = FieldNames{nField};
    set(hFig,'Name', [FieldName ' C CC Distribution']);
    Hist_H = histc(stCCDist.(FieldName).HCCVal,HistEdge);
    Hist_L = histc(stCCDist.(FieldName).LCCVal,HistEdge);
    h=bar(HistEdge(:),[Hist_H(:),Hist_L(:)],'grouped');
    set(h(1),'EdgeColor','r','FaceColor','r');
    set(h(2),'EdgeColor','b','FaceColor','b');
    xlim([-1 1]);
    legend('High','Low');

    if ~isempty(strSavePath)
        hgsave(hFig,[strSaveDir '\' FieldName ' C CC Dist.fig']);
    end
end
close(hFig);

hFig= figure();
HistEdge= -1:0.1:1;
for nField = 1:2:length(FieldNames)
    clf();
    FieldName = FieldNames{nField};
    set(hFig,'Name', [FieldName ' A CC Distribution']);
    Hist_H = histc(stCCDist.(FieldName).HCCVal_A,HistEdge);
    Hist_L = histc(stCCDist.(FieldName).LCCVal_A,HistEdge);
    h=bar(HistEdge(:),[Hist_H(:),Hist_L(:)],'grouped');
    set(h(1),'EdgeColor','r','FaceColor','r');
    set(h(2),'EdgeColor','b','FaceColor','b');
    xlim([-1 1]);
    legend('High','Low');

    if ~isempty(strSavePath)
        hgsave(hFig,[strSaveDir '\' FieldName ' A CC Dist.fig']);
    end
end
close(hFig);