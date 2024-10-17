function fmr1CircTrack_x_decodeByThetaCycle_v2(group)
% function fmr1CircTrack_x_thetaDecodeByTime(group)
%
% PURPOSE:
%   Decode across the theta cycle by theta phase to better understand theta
%   sequence activity in WT and FXS rats.
%
% MMD (Based in part on code written by JBT and CZ for the Colgin Lab)
% 08/2024
% Colgin Lab

%% OPTIONS

saveOrNot = 1;

plotEachDay = 0;

%% INITALIZE

minFr = 1;

phiEdges = 0:20:360; %for estimating group 0 phase
% phiBins = 0:30:360;

% Decoding parameters
Fs = 1000; % spike raster
bayesWin = 40/1000;
bayesStep = 10/1000;

%spk parameters
minCell = 3; %min number of cells that need to participate - from Zheng et al. 2016
minSpks = 3; %min number of total spikes in a cycle


%theta cycle parameters - Wang et al. 2020
thMin = 80/1000;
thMax = 180/1000;

runThresh = 5; %cm/s

degBinCtrs = group(2).rat(1).day(1).binCtrs; %doesn't change across days/rats
radBinCtrs = deg2rad(degBinCtrs)';

saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\THETA_SEQUENCES\thetaCycleDecodingByTime';
dataDir = 'E:\FMR1_CIRCTRACK\RAW_DATA';

groupNames = {'WT', 'FXS'};
cols = {'Blue', 'Red'};

radCmConv = 50;

if saveOrNot == 1 && plotEachDay == 1
    cd(saveDir)
    cd('plotEachDay')
end

%% CALCULATE DATA

fprintf('CALCULATING SEQUENCES PROPERTIES\n');

for g = 1:2
    fprintf('Group %d\n', g);
    for r = 1:length(group(g).rat)
        ratIn = 0;
        fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);

        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day))

            %             seqPxnExt = [];
            %             seqPxn = {};

            %             seqByDay = [];
            %             lfpByDay = [];

            if group(g).rat(r).day(d).decThresh == 0
                continue %o next day
            end %check decoding threshold
            ratIn = 1;

            %get ratemaps for day
            rateMaps = zeros(length(group(g).rat(r).day(d).xBegUnitInfo),...
                length(group(g).rat(r).day(d).xBegUnitInfo(1).smRateMap));

            badU = [];
            uIDs = zeros(length(group(g).rat(r).day(d).xBegUnitInfo),2);
            for u = 1:length(group(g).rat(r).day(d).xBegUnitInfo)
                if max(group(g).rat(r).day(d).xBegUnitInfo(u).smRateMap)>= minFr
                    rateMaps(u,:) = group(g).rat(r).day(d).xBegUnitInfo(u).smRateMap; %Smoothed ratemap
                    uIDs(u,:) = group(g).rat(r).day(d).xBegUnitInfo(u).ID;
                else
                    badU = [badU u]; %#ok
                end
            end
            rateMaps(badU,:) = [];
            rateMaps(rateMaps==0) = 0.0001; %get rid of zeros because our Bayesian decoder can't handle 'em.
            uIDs(badU,:) = [];

            % need to get theta global phase 0 - do it Zheng et al. 2021 way, cutting cycles at phase w least # of spikes
            tetNum = group(g).rat(r).day(d).thetaTet;

            if isfile([dataDir '\' group(g).name '\' group(g).rat(r).name '\' group(g).rat(r).day(d).name '\CSC' num2str(tetNum) 'ZeroPhase_v2.mat'])
                load([dataDir '\' group(g).name '\' group(g).rat(r).name '\' group(g).rat(r).day(d).name '\CSC' num2str(tetNum) 'ZeroPhase_v2.mat'],...
                    'phiHist', 'phiEdges', 'zeroPhase');
            else
                spkPhases = [];
                for b = 1:length(group(g).rat(r).day(d).begin)
                    begSpks = [];

                    radPos = group(g).rat(r).day(d).begin(b).radPos;
                    instRs = get_runspeed_circtrack(radPos);
                    smRs = smooth_runspeed(instRs);

                    for u = 1:length(group(g).rat(r).day(d).xBegUnitInfo)
                        if max(group(g).rat(r).day(d).xBegUnitInfo(u).smRateMap)>= minFr

                            velSpkTms = velocity_filter_spikes(group(g).rat(r).day(d).begin(b).unit(u).spkTms, smRs, runThresh, inf);
                            %                         begSpks = [begSpks; group(g).rat(r).day(d).begin(b).unit(u).spkTms];
                            begSpks = [begSpks; velSpkTms];
                        end %min fr
                    end %unit
                    lfpStruct = read_in_lfp([group(g).rat(r).day(d).begin(b).dir '\CSC' num2str(tetNum) '.ncs']); %need ts
                    load([group(g).rat(r).day(d).begin(b).dir '\CSC' num2str(tetNum) '_narrowThetaLfp.mat']);  %#ok
                    phiVctr = angle(hilbert(-filtLfp))*180/pi+180; %theta peaks = 0 and 360 - from CZ code

                    spkInds = match(begSpks, lfpStruct.ts);
                    spkPhases = [spkPhases; phiVctr(spkInds)];
                end %begin

                phiHist = histcounts(spkPhases, phiEdges);
                [~, minInd] = min(phiHist);
                zeroPhase = phiEdges(minInd);

                save([dataDir '\' group(g).name '\' group(g).rat(r).name '\' group(g).rat(r).day(d).name '\CSC' num2str(tetNum) 'ZeroPhase_v2.mat'],...
                    'zeroPhase', 'phiHist', 'phiEdges');
            end %get zero phase or not

            for b = 1:length(group(g).rat(r).day(d).begin)
                fprintf('\t\t\tBegin %d\n', b)

                if  isfile([group(g).rat(r).day(d).begin(b).dir '\seqData_decodeByThetaCycle_com_v4.mat'])
                    continue
                end
seqData = [];

                tmpInd = 0;
                %                 tmpInd = length(seqData);
                %                 allTms = vertcat(seqData(:).tms);

                radPos = group(g).rat(r).day(d).begin(b).radPos;
                instRs = get_runspeed_circtrack(radPos);
                smRs = smooth_runspeed(instRs);

                lfpStruct = read_in_lfp([group(g).rat(r).day(d).begin(b).dir '\CSC' num2str(tetNum) '.ncs']);
                load([group(g).rat(r).day(d).begin(b).dir '\CSC' num2str(tetNum) '_narrowThetaLfp.mat']);  %okay

                phiVctr = angle(hilbert(filtLfp))*180/pi+180; %troughs are 0 for this one

                phiVctrGZ = wrapTo360(phiVctr - zeroPhase);
                % phiVctrGZ = phiVctr;

                %get peak times of the pahse vector - 0 = 360 degrees so peak = trough!
                [~, trLocs] = findpeaks(phiVctrGZ, 'MinPeakHeight', 357); %troughs at 0/360
                trTms = lfpStruct.ts(trLocs);

                load([group(g).rat(r).day(d).begin(b).dir '\seqInfo_allBeg.mat'], 'spkRstr', 'timeMap')
                w = waitbar(0, 'Running through cycles...', 'Name', 'thetaCycles');
                for i = 1:length(trTms)
                    waitbar(i/length(trTms), w, sprintf('Cycle %d of %d', i, length(trTms)));

                    if trTms(i) > timeMap(end)
                        continue
                    end %all fits

                    seqStart = find(phiVctrGZ(1:trLocs(i)) < 180, 1, 'Last'); %inds
                    seqEnd = find(phiVctrGZ(trLocs(i)+1:end) > 180, 1, 'First') + trLocs(i);

                    seqStartTm = lfpStruct.ts(seqStart);
                    seqEndTm = lfpStruct.ts(seqEnd);
                    if isempty(seqStartTm) || isempty (seqEndTm)
                        continue
                    end

                    if (seqEndTm - seqStartTm) < thMin || (seqEndTm - seqStartTm) > thMax
                        continue
                    end

                    %check speed
                    tmpSpd = mean(smRs(match(seqStartTm, smRs(:,1)):match(seqEndTm, smRs(:,1)),2));
                    tmpPos = radPos(match(seqStartTm, radPos(:,1)):match(seqEndTm, radPos(:,1)),2);
                    actPos = wrapTo360(rad2deg(circ_mean(deg2rad(tmpPos))));

                    if tmpSpd < runThresh
                        continue %to next cycle
                    end %doesn't meet speed

                    %check spks/cells

                    startInd = match(seqStartTm, timeMap);
                    endInd = match(seqEndTm, timeMap)-1;
                    if endInd > size(spkRstr,2) %sampling rate used to get the spike raster is lower, so some are cut off
                        continue
                    end %doesn't fit

                    pullRstr = full(spkRstr(:,startInd:endInd)); %across entire cycle
                    spkSum = sum(pullRstr,2);

                    if sum(spkSum) < minSpks || nnz(spkSum) < minCell
                        continue %to next sequence
                    end %doesn't meet cell/spk

                    %                     if ismember(seqStartTm, allTms(:,1))
                    %                         continue
                    %                     end
                    pxn = BayesianDecoder(pullRstr, rateMaps, bayesWin, bayesStep, Fs);

                    [nWin, winStartInds] = find_num_windows(size(pullRstr,2), bayesWin*Fs, bayesStep*Fs);
                    if winStartInds(end)+bayesWin*Fs < size(pullRstr,2)
                        nWin = nWin + 1;
                        %                         winStartInds(end+1) = winStartInds(end)+bayesStep*Fs; %okay
                    end
                    %                     winStartTms = seqStartTm-timeBnd + winStartInds/Fs - 1/Fs;
                    %                     winEndTms = winStartTms + bayesWin;
                    pxn = pxn(:,1:nWin);

                    tmpBin = ~isnan((pxn(1,:)));
                    goodBins = find(tmpBin);
                    con = bwconncomp(tmpBin);
                    if con.NumObjects == 1
                        bins2use = con.PixelIdxList{1}';
                        tSpan = bayesWin + bayesStep*length(bins2use);
                    else
                        % lenObj = cellfun(@length, con.PixelIdxList);
                        % useInd = max(lenObj);
                        tmpFunc = @(x)(x(1));
                        findFir = cellfun(tmpFunc, con.PixelIdxList);
                        minOne = findFir - 1;
                        minOne = minOne(minOne>0);
                        tmpBin(minOne) = 1;

                        con = bwconncomp(tmpBin);
                        lenObj = cellfun(@length, con.PixelIdxList);
                        [~, useInd] = max(lenObj);
                        bins2use = con.PixelIdxList{useInd}';
                        if ~ismember(bins2use(1), goodBins)
                            bins2use = bins2use(2:end);
                        end
                        tSpan = bayesWin + bayesStep*length(bins2use);
                    end %whether we have on or it's more complicated

                    %                     bins2use = find(~isnan(pxn(1,:)));
                    if length(bins2use) == 1
                        continue
                    end

                    tmpInd = tmpInd + 1;
                    seqData(tmpInd).tms = [seqStartTm seqEndTm];
                    seqData(tmpInd).actCell =  uIDs(find(spkSum),:);
                    seqData(tmpInd).spkRstr =  pullRstr;
                    seqData(tmpInd).actPos = actPos;
                    seqData(tmpInd).pxn = pxn;
                    seqData(tmpInd).bins2use = bins2use;
                    seqData(tmpInd).tSpan = tSpan;

                    tBins = 0:bayesStep:bayesStep*size(pxn,2);
                    [seqData(tmpInd).p_resid, seqData(tmpInd).p_r2, seqData(tmpInd).r2, seqData(tmpInd).slope,...
                        seqData(tmpInd).xSpan, seqData(tmpInd).calphase, seqData(tmpInd).propPostClose, seqData(tmpInd).com] =...
                        get_seq_pVal_com(pxn, radBinCtrs, tBins, bins2use);

                end %i

                save([group(g).rat(r).day(d).begin(b).dir '\seqData_decodeByThetaCycle_com_v4.mat'], 'seqData')
                close(w)
            end %begin

            %             save([dataDir '\' group(g).name '\' group(g).rat(r).name '\' group(g).rat(r).day(d).name...
            %                 '\seqData_decodeByThetaCycle_v2.mat'], 'seqData')
            %             keyboard
        end %day

    end %rat
end %group
fprintf('\n')
cd(saveDir)

%% ADD PROP TO ACT

% for g = 1:2
%     fprintf('Group %d\n', g);
%     for r = 1:length(group(g).rat)
%         fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);
%         for d = 1:length(group(g).rat(r).day)
%             if group(g).rat(r).day(d).decThresh == 0
%                 continue %o next day
%             end %check decoding threshold
%             for b = 1:4
%                 load([group(g).rat(r).day(d).begin(b).dir '\seqData_decodeByThetaCycle_com_v4.mat'], 'seqData')
% 
%                 for i = 1:length(seqData)
%                     
%                     seqData(i).minDist = min(wrapTo2Pi(abs(circ_dist(deg2rad(seqData(i).actPos), seqData(i).calphase))));
% if seqData(i).minDist > pi
%     keyboard
% end
%                 end
% 
%   save([group(g).rat(r).day(d).begin(b).dir '\seqData_decodeByThetaCycle_com_v4.mat'], 'seqData')
%             end %begin
%         end %day
%     end %rat
% end %group

%% GET DATA AND PLOT

fprintf('GATHERING DATA FOR PLOTTING\n')

trajSlopes = cell(2,1);
trajX = cell(2,1);
trajT = cell(2,1);
trajR = cell(2,1);

trajSlopesRat = cell(2,1);
trajSlopesRat(:) = {cell(1,6)};
trajXRat = cell(2,1);
trajXRat(:) = {cell(1,6)};
trajTRat = cell(2,1);
trajTRat(:) = {cell(1,6)};
trajRRat = cell(2,1);
trajRRat(:) = {cell(1,6)};


flatSlopes = cell(2,1);
flatX = cell(2,1);
flatT = cell(2,1);
flatR = cell(2,1);
flatThresh = 0.01; %r2 thresh for flat slopes

dataForGLM = [];
errorAxis = 0:.05:pi;
radThresh = 0.35;
radThreshInd = match(radThresh, errorAxis);

rValPull = zeros(3,2,2); %pulling random examples for trajectory and flat slopes for plot
pxnPull = cell(3,2,2);
slopePull = zeros(3,2,2);
ratNamePull = cell(3,2,2);

seqPerRat = zeros(6,2); %for creating table w number of events per rat

for g = 1:2
    fprintf('Group %d\n', g);
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);
        for d = 1:length(group(g).rat(r).day)
            if group(g).rat(r).day(d).decThresh == 0
                continue %o next day
            end %check decoding threshold
            for b = 1:4
                load([group(g).rat(r).day(d).begin(b).dir '\seqData_decodeByThetaCycle_com_v4.mat'], 'seqData')

                load([dataDir '\' group(g).name '\' group(g).rat(r).name '\' group(g).rat(r).day(d).name '\dayDecodeErr.mat'], 'cumErr')

                forInds = vertcat(seqData(:).slope)>0;
                forSeq = seqData(forInds); %forward sequences

                % get sig traj slopes
                sigTrajInds = find(vertcat(forSeq(:).p_r2)<0.05&vertcat(forSeq(:).propPostClose)>=0.6&vertcat(forSeq(:).minDist)<=0.35);
                seqPerRat(r,g) = seqPerRat(r,g) + length(sigTrajInds);
                if isempty(sigTrajInds)
                    continue
                end

%                 storeInd = []; %check from randomly sampled examples indices
%                 if g == 1 && r == 1 && d == 2 && b == 4
%                     storeInd = 1;
%                     it = 33;
%                     fi = 803;
%                 end
%                 if g == 1 && r == 5 && d == 2 && b == 3
%                     storeInd = 2;
%                     it = 441;
%                     fi = 264;
%                 end
%                 if g == 1 && r == 6 && d == 1 && b == 1
%                     storeInd = 3;
%                     it = 155;
%                     fi = 356;
%                 end
%                 if g == 2 && r == 2 && d == 1 && b == 2
%                     storeInd = 1;
%                     it = 873; %CHANGE
%                     fi = 1856;
%                 end
%                 if g == 2 && r == 4 && d == 3 && b == 3
%                     storeInd = 2;
%                     it = 1051;
%                     fi = 856;
%                 end
%                 if g == 2 && r == 6 && d == 1 && b == 1
%                     storeInd = 3;
%                     it = 1197;
%                     fi = 884;
%                 end
% 
%                 if ~isempty(storeInd)
% 
%                     rValPull(storeInd,g,1) = forSeq(it).r2;
%                     pxnPull{storeInd,g,1} = shift_pxn(forSeq(it).pxn, forSeq(it).actPos, spatBinSz);
%                     slopePull(storeInd,g,1) = forSeq(it).slope.*radCmConv;
%                     ratNamePull{storeInd,g,1} = group(g).rat(r).name;
% 
%                     rValPull(storeInd,g,2) = forSeq(fi).r2;
%                     pxnPull{storeInd,g,2} = shift_pxn(forSeq(fi).pxn, forSeq(fi).actPos, spatBinSz);
%                     slopePull(storeInd,g,2) = forSeq(fi).slope.*radCmConv;
%                     ratNamePull{storeInd,g,2} = group(g).rat(r).name;
%                 end %store trajectory example


                trajSlopes{g} = [trajSlopes{g}; vertcat(forSeq(sigTrajInds).slope).*radCmConv];
                trajSlopesRat{g}{r} = [trajSlopesRat{g}{r}; vertcat(forSeq(sigTrajInds).slope).*radCmConv];
                trajX{g} = [trajX{g}; abs(vertcat(forSeq(sigTrajInds).xSpan)).*radCmConv];
                trajXRat{g}{r} = [trajXRat{g}{r}; abs(vertcat(forSeq(sigTrajInds).xSpan)).*radCmConv];
                trajR{g} = [trajR{g}; vertcat(forSeq(sigTrajInds).r2)];
                trajRRat{g}{r} = [trajRRat{g}{r}; vertcat(forSeq(sigTrajInds).r2)];

                tmpT = (cellfun(@length, vertcat({forSeq(sigTrajInds).bins2use}))' * bayesStep) + bayesWin;
                trajT{g} = [trajT{g}; tmpT];
                trajTRat{g}{r} = [trajTRat{g}{r}; tmpT];

                tmpGLMData = cat(2, repmat(g, length(sigTrajInds), 1), repmat(cumErr(radThreshInd), length(sigTrajInds), 1),...
                    vertcat(forSeq(sigTrajInds).slope).*radCmConv,  vertcat(forSeq(sigTrajInds).xSpan).*radCmConv, tmpT);
                dataForGLM = cat(1, dataForGLM, tmpGLMData);

                sigFlatInds =  find(vertcat(forSeq(:).r2)>flatThresh & vertcat(forSeq(:).p_resid)<0.5 & vertcat(forSeq(:).propPostClose)>=0.5);

                flatSlopes{g} = [flatSlopes{g}; vertcat(forSeq(sigFlatInds).slope).*radCmConv];
                flatX{g} = [flatX{g}; abs(vertcat(forSeq(sigFlatInds).xSpan)).*radCmConv];
                flatR{g} = [flatR{g}; vertcat(forSeq(sigFlatInds).r2)];

                tmpT = (cellfun(@length, vertcat({forSeq(sigFlatInds).bins2use}))' * bayesStep) + bayesWin;
                flatT{g} = [trajT{g}; tmpT];
            end %begin
        end %day
    end %rat
end %group
keyboard
%% MAKE EXAMPLES PLOTS

figtitle = 'ThetaSequences_examples_trajectory_v2';
%     figtitle = [figtitle '_' methodNames{detectMethod}];
figure('Name', figtitle, 'Position', [248 448 1349 532])
spMap = [1:3; 4:6];

for g = 1:2

    tmpPxnForPlot = pxnPull(:,g,1);
    tmpR2ForPlot = rValPull(:,g,1);
    tmpSlopeForPlot = slopePull(:,g,1);
    rNamesForPlot = ratNamePull(:,g,1);

    [~, sortOrd] = sort(tmpR2ForPlot, 'descend');

    for i = 1:3
        subplot(2,3,spMap(g,find(sortOrd == i)))

        maxTm = size(tmpPxnForPlot{i},2) * bayesStep + bayesWin/2;

        imagesc(0:1, -180:180, tmpPxnForPlot{i})
        axis xy
        axis square
        colormap jet
        caxis([0 0.3])
        hold on
        plot([0 1], [0 0], 'Color', 'w', 'LineStyle', '--')

        xticks([0 1])
        xticklabels({'', num2str(maxTm)})
        if g == 2
            xlabel('Time (s)')
        end

        ylim([-90 90])
        yticks(-180:90:180)
        if find(sortOrd == i) == 1
            ylabel({groupNames{g}, 'Angular position (deg)'})
        end

        title({rNamesForPlot{i}; ['r^2 = ' num2str(round(tmpR2ForPlot(i),2)) ' | slope = ' num2str(round(tmpSlopeForPlot(i))) ' cm/s']})

    end %which map

end %group
cbr = colorbar;
set(cbr, 'Position', [.92 .59 .01 .3])
ylabel(cbr, 'Probability')

same_caxis(0);

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end


%% FIG - SLOPE DISTRIBUTION

figtitle = 'SlopeDistribution_trajectory';
figure('Name', figtitle)

minData = 0;
maxData = 3000;
delta = 125;
nboot = 1000;

leg = cell(2,1);
for g = 1:2
    % subplot(1,2,g)
    hold on
    axis square
    % histogram(trajSlopes{g}, 'DisplayStyle', 'stairs')
    [wDist, xAx] = WeightedProportion(trajSlopes{g}, minData, maxData, delta);
    %     xAx = xAx.*radCmConv;
    tmpFunc = @(x)(WeightedProportion(x, minData, maxData, delta));
%     plot(xAx, wDist, 'Color', rgb(cols{g}))

    CI = bootci(nboot, {tmpFunc, trajSlopes{g}}, 'type', 'per');

    lh(g) = plot_filled_ci(xAx, wDist, CI, rgb(cols{g}));

    leg{g} = [groupNames{g} ' ( n = ' num2str(length( trajSlopes{g})) ')'];
end %group

     xlim([0 2000])
xlabel('Slope (cm/s)')
ylabel('Proportion')
ylim([0 0.05])
legend(lh, leg)

[p, observeddifference, effectsize] = permutationTest(trajSlopes{2}, trajSlopes{1}, 5000, 'showprogress', 1, 'plotresult', 0);

title(['p = ' num2str(p)])


if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end
% 
% tbl = table(dataForGLM(:,1), dataForGLM(:,2), dataForGLM(:,3), 'VariableNames', {'genotype', 'error', 'slope'}); %from Ernie
% tbl.genotype = categorical(tbl.genotype);
% modelspec = 'slope ~ genotype*error';
% mdl =  fitglm(tbl, modelspec, 'Distribution', 'gamma', 'link', 'log');

%% SLOPE BY RAT

figtitle = 'SlopeByRat_trajectory';
figure('Name', figtitle)

ratCntr = 0;
nboot = 1000;
tmpTicks = [];
tickLabs = {};

for g = 1:2
    for r = 1:length(trajSlopesRat{g})
        if isempty(trajSlopesRat{g}{r})
            continue
        end
        hold on
        ratCntr = ratCntr + 0.1;
        xData = ratCntr;
        tmpTicks = [tmpTicks xData];
        tickLabs = cat(1, tickLabs, group(g).rat(r).name);

        tmpCI = bootci(nboot, {@mean, trajSlopesRat{g}{r}}, 'type', 'per');
        yData = mean(trajSlopesRat{g}{r});

        lowErr = yData - min(tmpCI);
        highErr = max(tmpCI) - yData;
        er = errorbar(xData, yData, lowErr, highErr, 'o', 'Color', rgb(cols{g}));
        er.CapSize = 0;

    end %rat
    ratCntr = ratCntr + 0.25;
end %group

xlim([tmpTicks(1)-0.1 tmpTicks(end)+0.1])
xticks(tmpTicks)
xticklabels(tickLabs)
xlabel('Rat #')


%% ALL X

figtitle = 'xDistribution_trajectory';
figure('Name', figtitle)

minData = 0;
maxData = 300;
delta = 10;
nboot = 1000;

leg = cell(2,1);
for g = 1:2
    % subplot(1,2,g)
    hold on
    axis square

    [wDist, xAx] = WeightedProportion(trajX{g}, minData, maxData, delta);
    %     xAx = xAx.*radCmConv;
    tmpFunc = @(x)(WeightedProportion(x, minData, maxData, delta));
%      plot(xAx, wDist, 'Color', rgb(cols{g}))

    CI = bootci(nboot, {tmpFunc, trajX{g}}, 'type', 'per');

    lh(g) = plot_filled_ci(xAx, wDist, CI, rgb(cols{g}));

    leg{g} = [groupNames{g} ' ( n = ' num2str(length( trajX{g})) ')'];
end %group

xlim([0 125])
xticks([0:25:125])
xlabel('x-span (cm)')
ylabel('Proportion')
ylim([0 0.05])
legend(lh, leg)

[p, observeddifference, effectsize] = permutationTest(trajX{2}, trajX{1}, 5000, 'showprogress', 1, 'plotresult', 0);
title(['p = ' num2str(p)])

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end


%% TRAJ T

figtitle = 'tDistribution_trajectory';
figure('Name', figtitle)

minData = 0;
maxData = 0.250;
delta = 0.01;
nboot = 1000;

leg = cell(2,1);
for g = 1:2
    % subplot(1,2,g)
    hold on
    axis square

    [wDist, xAx] = WeightedProportion(trajT{g}, minData, maxData, delta);
    %     xAx = xAx.*radCmConv;
    tmpFunc = @(x)(WeightedProportion(x, minData, maxData, delta));
%      plot(xAx, wDist, 'Color', rgb(cols{g}))

    CI = bootci(nboot, {tmpFunc, trajT{g}}, 'type', 'per');

    lh(g) = plot_filled_ci(xAx, wDist, CI, rgb(cols{g}));

    leg{g} = [groupNames{g} ' ( n = ' num2str(length( trajT{g})) ')'];
end %group

xlim([0 0.25])
xlabel('t-span (s)')
ylabel('Proportion')
% ylim([0 0.05])
legend(lh, leg, 'Location', 'northwest')

[p, observeddifference, effectsize] = permutationTest(trajT{2}, trajT{1}, 5000, 'showprogress', 1, 'plotresult', 0);
title(['p = ' num2str(p)])

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end

%% R2

figtitle = 'r2Distribution_trajectory';
figure('Name', figtitle)

minData = 0;
maxData = 1;
delta = 0.1;
nboot = 1000;

leg = cell(2,1);
for g = 1:2
    % subplot(1,2,g)
    hold on
    axis square

    [wDist, xAx] = WeightedProportion(trajR{g}, minData, maxData, delta);
    %     xAx = xAx.*radCmConv;
    tmpFunc = @(x)(WeightedProportion(x, minData, maxData, delta));
%      plot(xAx, wDist, 'Color', rgb(cols{g}))

    CI = bootci(nboot, {tmpFunc, trajR{g}}, 'type', 'per');

    lh(g) = plot_filled_ci(xAx, wDist, CI, rgb(cols{g}));

    leg{g} = [groupNames{g} ' ( n = ' num2str(length( trajR{g})) ')'];
end %group

%      xlim([0 0.25])
xlabel('r2')
ylabel('Proportion')
ylim([0 0.035])
legend(lh, leg, 'Location', 'northwest')
[p, observeddifference, effectsize] = permutationTest(trajR{2}, trajR{1}, 5000, 'showprogress', 1, 'plotresult', 0);
title(['p = ' num2str(p)])

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end
% 
% %% DOWN SAMPLE
% 
% nByG = cellfun(@length, trajSlopes);
% [minVal, minG] = min(nByG);
% [maxVal, maxG] = max(nByG);
% s = RandStream('mt19937ar', 'Seed', 1); 
% keepInds = randsample(s, maxVal, minVal, false);
% 
% p = permutationTest(trajSlopes{maxG}(keepInds), trajSlopes{minG}, 5000, 'showprogress', 0, 'plotresult', 0);
% fprintf('Random down sample permutation test - effect of genotypes on slope: p = %d\n', p)
% p = permutationTest(trajX{maxG}(keepInds), trajX{minG}, 5000, 'showprogress', 0, 'plotresult', 0);
% fprintf('Random down sample permutation test - effect of genotypes on x-span: p = %d\n', p)
% p = permutationTest(trajT{maxG}(keepInds), trajT{minG}, 5000, 'showprogress', 0, 'plotresult', 0);
% fprintf('Random down sample permutation test - effect of genotypes on t-span: p = %d\n', p)
% p = permutationTest(trajR{maxG}(keepInds), trajR{minG}, 5000, 'showprogress', 0, 'plotresult', 0);
% fprintf('Random down sample permutation test - effect of genotypes on r2: p = %d\n', p)


end %function