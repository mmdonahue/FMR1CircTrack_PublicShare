function fmr1CircTrack_x_plotSingleCell_replay(group)
% function fmr1CircTrack_x_plotSingleCell_replay(group)
%
% PURPOSE:
%   To plot single cell analyses for replay events from WT and KO rats.
%
% INPUTS:
%   group: data struct
%
% OPTIONS:
%   saveOrNot: save figs (1) or don't (0)
%   See function for other options.
%
% MMD
% 7/2021
% Colgin Lab

%% OPTIONS

saveOrNot = 1;

prepForStats = 1;

downSampEvents = 0; %so number of events equal between groups

% Binned firing rate options
binSz = 0.001; %seconds, for binning firing rate
gWinStd = 10/1000; %5 ms in Ernie's paper
% gWinDur = gWinStd;
gWinDur = gWinStd * 6;
nboot = 1000; %for calculating confidence intervals, from Hwaun & Colgin 2019

%% INITIALIZE

% Binned firing rate stuff
totTime = 1; %seconds, to plot binned firing rate (with event in middle)
preTime = [-0.4 -0.1]; %time frame before ripple to get normalizing firing rate; 0.4-0.1s in Hwaun & Colgin 2017
numBins = totTime/binSz;
binTms = -totTime/2:binSz:totTime/2-binSz;
onsetBin = ceil(numBins/2);
preBins = match(preTime, binTms)';

binFr = cell(2,1); %by group
binFrPreNorm = cell(2,1); %by group - firing rate normalized by firing avg firing rate 0.4-0.1 s before replay event

binFrAct = cell(2,1); %by group
binFrNormAct = cell(2,1); %by group - firing rate normalized by firing avg firing rate 0.4-0.1 s before replay event

gWinDur = gWinDur/binSz; %convert based on bin size
gWinStd = gWinStd/binSz; %convert based on bin size
gKrnl = gausskernel(gWinDur, gWinStd);

% FR and SPK/CELL stuff
frInEv = cell(2,1);
frInEv(:) = {cell(1,6)};
spkPerCell = cell(2,1);
spkPerCell(:) = {cell(1,6)};

% Misc
cols = {'Blue', 'Red'};
groupNames = {'WT', 'FXS'};
markers = {'o', 'square', '*', '+', 'x', 'diamond', '^', '+'};
fCols = cat(1, [1 1 1], rgb('Blue'), [1 1 1], [1 1 1], [1 1 1], rgb('Red'), [1 1 1], [1 1 1]);

evCntr = zeros(1,2);
r2Thresh = 0.5;
propCloseThresh = 0;
jumpThresh = 0.25*2*pi; %one quarter of the track

saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\REPLAY\replayEvents_singleCell';

if prepForStats == 1
    statAll = [];
end %stats

%% PREP FOR DOWN SAMPLE
%
% if downSampEvents == 1
%     rsE = RandStream('mt19937ar');
%
%     gNum = zeros(1,2);
%
%     for g = 1:2
%         for r = 1:length(group(g).rat)
%             for d = 1:length(group(g).rat(r).day)
%                 for s = 2:5
%                     if ~isempty(group(g).rat(r).day(d).sleep(s).unit)
%                         if detectMethod == 1
%                             events = group(g).rat(r).day(d).sleep(s).popEv;
%                         else
%                             events = group(g).rat(r).day(d).sleep(s).rip;
%                         end %detcet method
%                         for i = 1:length(events)
%                             if ~isnan(events(i).r2) && events(i).r2 >= r2Thresh
%                                 gNum(g) = gNum(g) + 1;
%                             end %meets thresh
%                         end %i - pop events
%                     end %if there are events
%                 end %sleep
%             end %day
%         end %rat
%     end %group
%
%     [minEv, minG] = min(gNum);
%     keepEvs = datasample(rsE, 1:max(gNum), minEv, 'Replace', false);
%     keepEvs = sort(keepEvs);
%
%     allEvCntr = 0; %initialize
%
% end %if down sampling events

%% GET DATA

for g = 1:2
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
            if group(g).rat(r).day(d).decThresh == 0 %if meets decoding threshold
                continue
            end

            binFrByU = [];
            spkPerByU = [];
            frByU = [];
            %             numEvents = 0;
            %             if downSampEvents == 1 && g ~= minG %find out which events from this day to keep
            %
            %                 for s = 2:5 %getting getting total number of events for day
            %                     if ~isempty(group(g).rat(r).day(d).sleep(s).unit)
            %                         for i = 1:length(group(g).rat(r).day(d).sleep(s).popEv)
            %                             if group(g).rat(r).day(d).sleep(s).popEv(i).r2 >= r2Thresh
            %                                 numEvents = numEvents + 1;
            %                             end %meets threshold
            %                         end %i
            %                     end %if there is anything in popEvents
            %                 end %sleep
            %
            %                 dayEvs = 1:numEvents; %initialize
            %                 keepInds = find(keepEvs > allEvCntr & keepEvs <= allEvCntr+numEvents); %which from keepEvs are in this day (no = in first since allEvCntr starts at 0)
            %                 dayEvInds = keepEvs(keepInds) - allEvCntr;
            %                 dayEvs = dayEvs(dayEvInds); %the events from this day that get to be included
            %
            %                 allEvCntr = allEvCntr + numEvents;
            %             end %down sample

            for s = 2:5
                if isempty(group(g).rat(r).day(d).sleep(s).unit) %if there is anything in the sleep
                    continue
                end

                load([group(g).rat(r).day(d).sleep(s).dir '\PopEventOut_CA1'], 'spkRstr', 'timeMap')

                events = group(g).rat(r).day(d).sleep(s).popEv;
                for i = 1:length(events)
                    r2 = events(i).r2;
                    if r2 < r2Thresh || isnan(r2) || events(i).maxJump > jumpThresh || events(i).propClose < propCloseThresh %it doesn't meet threshold
                        continue %to next event
                    end %check r2

                    startInd = match(events(i).tms(1)-totTime/2, timeMap);
                    endInd = match(events(i).tms(1)+totTime/2, timeMap);

                    pullRstr = full(spkRstr(:, startInd:endInd));
                    tmpBinFr = pullRstr ./ binSz; %convert to Hz
                    binFrByU = cat(3, binFrByU, tmpBinFr); %events in third dimension

                    matchEnd = match(events(i).tms(2), timeMap);

                    sumSpks = sum(full(spkRstr(:,startInd+onsetBin:matchEnd)),2);
                    spkPerByU = cat(2, spkPerByU, sumSpks);
                    frByU = cat(2, frByU, sumSpks ./ (diff(events(i).tms)));

                end %i - events
            end %sleep

            tmpBinFr = mean(binFrByU,3); %mean for each unit
            smBinFr = convn(tmpBinFr, gKrnl', 'same'); %smooth across bin fr for each cell individually
            binFr{g} = [binFr{g}; smBinFr];

            pull100bin = smBinFr(:,onsetBin:(0.1/binSz)+onsetBin);
            mean100bin = mean(pull100bin,2);

            pullNorms = mean(smBinFr(:,preBins(1):preBins(2)),2);
            pullNorms(pullNorms<0.05) = 0.05; %effectively 0
            tmpNormFr = (smBinFr ./ pullNorms) * 100;
            binFrPreNorm{g} = [binFrPreNorm{g}; tmpNormFr];

            pull100norm = tmpNormFr(:,onsetBin:(0.1/binSz)+onsetBin);
            mean100norm = mean(pull100norm,2);

            spkPerByU(spkPerByU==0) = nan; %active only
            tmpSpk = nanmean(spkPerByU,2);
            spkPerCell{g}{r} = [spkPerCell{g}{r}; tmpSpk(~isnan(tmpSpk))];

            frByU(frByU==0) = nan; %active only
            tmpFr = nanmean(frByU,2);
            frInEv{g}{r} = [frInEv{g}{r}; tmpFr(~isnan(tmpFr))];

            store100Act = [];
            store100ActNorm = [];
            store200Act = [];
            storePeak = [];
            for u = 1:size(smBinFr,1)
                tmpActBinFr = binFrByU(u,:,~isnan(frByU(u,:)));
                if isnan(sum(tmpActBinFr(:))) || isempty(tmpActBinFr)
                    continue
                end
                smActBinFr = convn(mean(tmpActBinFr,3), gKrnl', 'same'); %smooth across bin fr for each cell individually
                binFrAct{g} = [binFrAct{g}; smActBinFr];

                pull100binAct = smActBinFr(:,onsetBin:(0.1/binSz)+onsetBin);
                mean100ActNorm = mean(pull100binAct,2);
                store100Act = [store100Act; mean100ActNorm];

                 pull200binAct = smActBinFr(:,(0.1/binSz)+onsetBin:(0.2/binSz)+onsetBin);
                mean200ActNorm = mean(pull200binAct,2);
                store200Act = [store200Act; mean200ActNorm];

                tmpActNormFr =  (smActBinFr ./ pullNorms(u)) * 100;
                binFrNormAct{g} = [binFrNormAct{g}; tmpActNormFr];

                pull100ActNorm = tmpActNormFr(:,onsetBin:(0.1/binSz)+onsetBin);
                mean100ActNorm = mean(pull100ActNorm,2);
                store100ActNorm = [store100ActNorm; mean100ActNorm];

                pullPk = max(smActBinFr(onsetBin:end));
                storePeak = [storePeak; pullPk];
            end %unit
            store200Act(store200Act<0.0001) = 0; % spps gets mad
            tmpUstat = (1:length(tmpFr(~isnan(tmpFr))))' + size(statAll,1);
            tmpStat = cat(2, repmat(g,sum(~isnan(tmpFr)),1), repmat(r,sum(~isnan(tmpFr)),1), repmat(d,sum(~isnan(tmpFr)),1),...
                tmpUstat, mean100bin(~isnan(tmpFr)), mean100norm(~isnan(tmpFr)), store100Act, store100ActNorm,...
                tmpFr(~isnan(tmpFr)), tmpSpk(~isnan(tmpSpk)), store200Act, storePeak);
            statAll = [statAll; tmpStat];

        end %day
    end %rat
end %group

cd(saveDir)
keyboard

%% FIG 1 - BINNED FR AROUND RIP

figtitle = 'BinnedFr_inReplayEvent';

if downSampEvents == 1
    figtitle = [figtitle '_downSampled_events'];
end %down samp
figure('Position', [333 501 1142 388], 'Name', figtitle)

yLabs = {'Firing rate (Hz)', 'Normalized firing rate (%)'};

for sp = 1:2
    subplot(1,2,sp)
    lh = zeros(1,2);
    leg = cell(1,2);

    for g = 1:2
        if sp == 1
            pullData = binFr{g};
        else
            pullData = binFrPreNorm{g};
        end %which subplot
        cellNum = size(pullData,1);

        mubin = mean(pullData);

        binCut = mubin(onsetBin-(totTime/2/binSz)+1:onsetBin+(totTime/2/binSz));

        CI = bootci(nboot, {@mean, pullData}, 'type', 'per');
        CIcut = CI(:,onsetBin-(totTime/2/binSz)+1:onsetBin+(totTime/2/binSz));

        hold on;

        lh(g) =  plot_filled_ci(binTms, binCut, CIcut, rgb(cols{g}));

        leg{g} = [groupNames{g} ' n = ' num2str(cellNum) ' cells'];
    end %group

    ylabel(yLabs{sp})
    xlabel('Time from replay event onset (s)')
    xlim([0 numBins])
    xlim([-0.4 0.4])

    legend(lh, leg, 'Location', 'northwest')
end %subplot

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end

%% FIG 1 - BINNED FR AROUND RIP - ACTIVE ONLY

figtitle = 'ActBinnedFr_inReplayEvent';

if downSampEvents == 1
    figtitle = [figtitle '_downSampled_events'];
end %down samp
figure('Position', [333 501 1142 388], 'Name', figtitle)

yLabs = {'Firing rate (Hz)', 'Normalized firing rate (%)'};

for sp = 1:2
    subplot(1,2,sp)
    lh = zeros(1,2);
    leg = cell(1,2);

    for g = 1:2
        if sp == 1
            pullData = binFrAct{g};
        else
            pullData = binFrNormAct{g};
        end %which subplot
        cellNum = size(pullData,1);

        mubin = mean(pullData);
        binCut = mubin(onsetBin-(totTime/2/binSz)+1:onsetBin+(totTime/2/binSz));

        CI = bootci(nboot, {@mean, pullData}, 'type', 'per');
        CIcut = CI(:,onsetBin-(totTime/2/binSz)+1:onsetBin+(totTime/2/binSz));

        hold on;
        lh(g) =  plot_filled_ci(binTms, binCut, CIcut, rgb(cols{g}));

        leg{g} = [groupNames{g} ' n = ' num2str(cellNum) ' cells'];
    end %group

    ylabel(yLabs{sp})
    xlabel('Time from replay event onset (s)')
    xlim([-0.25 0.5])

    legend(lh, leg, 'Location', 'northeast')
end %subplot

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end

%% FIG 2 - FR IN RIP

figtitle = ('Fr_byRat_CI');
if downSampEvents == 1
    figtitle = [figtitle '_downSampEvents'];
end %down samp

figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = {};
for g = 1:2
    for r = 1:length(frInEv{g})
        if isempty(frInEv{g}{r})
            continue
        end
        hold on
        ratCntr = ratCntr + 1;
        h = dotplot(ratCntr, {frInEv{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([ratCntr-0.5 ratCntr-0.5], [0 100], 'Color', 'k')
        line([ratCntr+0.5 ratCntr+0.5], [0 100], 'Color', 'k')
        tickLabs = cat(1, tickLabs, group(g).rat(r).name);

        quants = quantile(frInEv{g}{r}, 3);
        tmpCI = bootci(nboot, {@mean, frInEv{g}{r}}, 'type', 'per');
        meanData = mean(frInEv{g}{r});
        %         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
        pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
        line([ratCntr-0.25 ratCntr+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 9])
xticks(1:8)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0 80])
ylabel('Firing rate (Hz)')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end

%% FIG - FR IN RIP SORTED

figtitle = ('Fr_byRat_CI_sorted');
if downSampEvents == 1
    figtitle = [figtitle '_downSampEvents'];
end %down samp

figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = cell(8,1);
for g = 1:2
    tmpMeans = cellfun(@mean, frInEv{g});
    [~, sortOrd] = sort(tmpMeans);
    %     sortOrd = sortOrd - sum(isnan(tmpMeans));
    for r = 1:length(frInEv{g})
        if isempty(frInEv{g}{r})
            continue
        end
        hold on
        %         ratCntr = ratCntr + 1;
        xVal = find (r==sortOrd) + 4*(g-1);
        h = dotplot(xVal, {frInEv{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([xVal-0.5 xVal-0.5], [0 100], 'Color', 'k')
        line([xVal+0.5 xVal+0.5], [0 100], 'Color', 'k')
        tickLabs{xVal} = group(g).rat(r).name;

        quants = quantile(frInEv{g}{r}, 3);
        tmpCI = bootci(nboot, {@mean, frInEv{g}{r}}, 'type', 'per');
        meanData = mean(frInEv{g}{r});
        %         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
        pt = patch([xVal-0.25 xVal-0.25 xVal+0.25 xVal+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
        line([xVal-0.25 xVal+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 9])
xticks(1:8)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0 75])
ylabel('Firing rate (Hz)')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end

%% FIGS - SPK PER BY RAT CI SORTED

figtitle = ('SpkPer_byRat_CI_sorted');
if downSampEvents == 1
    figtitle = [figtitle '_downSampEvents'];
end %down samp

figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = cell(8,1);
for g = 1:2
    tmpMeans = cellfun(@mean, spkPerCell{g});
    [~, sortOrd] = sort(tmpMeans);
    %     sortOrd = sortOrd - sum(isnan(tmpMeans));
    for r = 1:length(spkPerCell{g})
        if isempty(spkPerCell{g}{r})
            continue
        end
        hold on
        %         ratCntr = ratCntr + 1;
        xVal = find (r==sortOrd) + 4*(g-1);
        h = dotplot(xVal, {spkPerCell{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([xVal-0.5 xVal-0.5], [0 10], 'Color', 'k')
        line([xVal+0.5 xVal+0.5], [0 10], 'Color', 'k')
        tickLabs{xVal} = group(g).rat(r).name;

        quants = quantile(spkPerCell{g}{r}, 3);
        tmpCI = bootci(nboot, {@mean, spkPerCell{g}{r}}, 'type', 'per');
        meanData = mean(spkPerCell{g}{r});
        %         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
        pt = patch([xVal-0.25 xVal-0.25 xVal+0.25 xVal+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
        line([xVal-0.25 xVal+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 9])
xticks(1:8)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0.5 7])
ylabel('Spikes/event, per cell')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end


end %function