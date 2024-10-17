function fmr1CircTrack_x_replayEvent_TFR(group)
% function fmr1CircTrack_x_replayEvent_TFR(group)
%
% PURPOSE:
%   This function gets the TFR from replay events from both genotypes as
%   well as the power x frequency plot, peak ripple requency. 
%
% INPUT:
%   group = data struct
%
%
% OPTIONS:
%   See function for options, including the pre- and post-event times,
%   whether or not to down sample events, and which event detection method
%   to use (among others).
%
% MMD
% 12/2021
% Colgin Lab

%% OPTIONS

saveOrNot = 1;

downSampEvents = 0;

nboot = 5000;

preEvTm = 0.4; %time pre-event
postEvTm = 0.5;

%% INITIALIZE

r2Thresh = 0.5;
propCloseThresh = 0;
jumpThresh = 0.25*2*pi; %one quarter of the track

fRange = [2 250];
fRangeCont = fRange(1):fRange(2);
% gamRange = [25 95]; %gamma
width = 7; %as in Mably et al. 2017

ripRange = [150 250];
ripInds = [ripRange(1)-(fRange(1)-1) ripRange(2)-(fRange(1)-1)];

gamRange = [25 55];
gamInds = [gamRange(1)-(fRange(1)-1) gamRange(2)-(fRange(1)-1)];

Fs = 2000; %for lfp
order = 5; %prewhiten
useOrigNorm = 0; %see get_wavelet_power code for more information
dbConv = 1; %don't convert to db (?)

normTms = [0.45 0.4]; %100-400 ms before event onset for getting the normalized power
tmBnds = -preEvTm:1/Fs:postEvTm;

allTFR = cell(1,2); %event TFR, normalized time by group
TFRbyRat = cell(1,2);
% allTFRbyEv = cell(2,2); %group x for/rev

PxF = cell(1,2); %power x freqeuncy from 1-100 ms after event detection (Guillespie et al. 2016)
PxFRat = cell(1,2);
PxFRat(:) = {cell(1,6)};
PxFbyEv = cell(2,2); %group x ev type
pkRipFreq = cell(1,2); %for storing peak ripple frequency across event types
pkRipFreqRat = cell(2,1);
pkRipFreqRat(:) = {cell(1,6)};
freqxSlope = cell(1,2);
gammaPowRat = cell(2,1);
gammaPowRat(:) = {cell(1,6)};

cols = {'Blue', 'Red'};
groupNames = {'WT', 'FXS'};

preWhiten = 1; %whether or not to pre-whiten the LFP (using prewhitening function in Matlab)

saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\REPLAY\replayEventTFR';

statAll = [];

%% PREP FOR DOWN SAMPLE

if downSampEvents == 1
    rsE = RandStream('mt19937ar');
    gNum = zeros(1,2);
    for g = 1:2
        for r = 1:length(group(g).rat)
            for d = 1:length(group(g).rat(r).day)
                for s = 2:5
                    events = group(g).rat(r).day(d).sleep(s).popEv; %shorten variable name

                    for i = 1:length(events)
                        if ~isnan(events(i).r2) && events(i).r2 >= r2Thresh
                            gNum(g) = gNum(g) + 1;
                        end %meets thresh
                    end %i - pop events
                end %sleep
            end %day
        end %rat
    end %group

    [minEv, minG] = min(gNum);
    keepEvs = datasample(rsE, 1:max(gNum), minEv, 'Replace', false);
    keepEvs = sort(keepEvs);

    allEvCntr = 0; %initialize

end %if down sampling events

%% INITIALIZE RAT FIG

figtitle = 'TFR_allRats';
figure('Name', figtitle, 'Position', [155 197 1659 771])

%% GET DATA

for g = 1:2
    fprintf('%s\n', group(g).name)
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d\n', r, length(group(g).rat))
        subplot(2, 6, (g-1)*6+r)
        ratTFR = [];
        for d = 1:length(group(g).rat(r).day)
            numEvents = 0;
            for s = 2:5 %getting getting total number of events for day
                if ~isempty(group(g).rat(r).day(d).sleep(s).unit)

                    events = group(g).rat(r).day(d).sleep(s).popEv; %shorten variable name

                    for i = 1:length(events)
                        if ~isnan(events(i).r2) && events(i).r2 >= r2Thresh
                            numEvents = numEvents + 1;
                        end %meets thresh
                    end %i - events
                end %if there is anything in popEvents
            end %sleep

            if downSampEvents == 1 && g ~= minG %find out which events from this day to keep
                dayEvs = 1:numEvents; %initialize
                keepInds = find(keepEvs > allEvCntr & keepEvs <= allEvCntr+numEvents); %which from keepEvs are in this day (no = in first since allEvCntr starts at 0)
                dayEvInds = keepEvs(keepInds) - allEvCntr;
                dayEvs = dayEvs(dayEvInds); %the events from this day that get to be included

                allEvCntr = allEvCntr + numEvents;
            end %down sample events

            if downSampEvents == 0 || g == minG
                fprintf('\t\tDay %d/%d - %d events\n', d, length(group(g).rat(r).day), numEvents)
            else
                fprintf('\t\tDay %d/%d - %d events (downsampled)\n', d, length(group(g).rat(r).day), length(dayEvs))
            end
            if numEvents == 0
                continue %to next day
            end % no events

            dayTFR = [];
            dayTFRbyEv = cell(1,2);
            dayPxF = [];
            dayPxFByEv = cell(1,2); %for or rev
            dayPkFreq = [];
            dayPkFreqByEv = cell(1,2);
            meanRipPow = [];
            meanGamPow = [];

%             tetNums = group(g).rat(r).day(d).tetNums;
            allIDs = vertcat(group(g).rat(r).day(d).xBegUnitInfo(:).ID);
            tetNums = mode(allIDs(:,1)); %use tetrode w the msot cells

            for tt = 1:length(tetNums)
                dayEvCntr = 0; %reset before we go through each tetrode

                fprintf('\t\t\tTetrode %d (%d/%d)\n', tetNums(tt), tt, length(tetNums))
                tmpTet = tetNums(tt);
                tetLFP = [];
                tetLFPbyEv = cell(2,1);
                tmpSlopes = [];


                for s = 2:5 %after experience
                    fprintf('\t\t\t\tSleep %d\n', s)
                    if isempty(group(g).rat(r).day(d).sleep(s).coords)
                        continue
                    end % no sleep data
                    cd(group(g).rat(r).day(d).sleep(s).dir)

                    lfpStruct = read_in_lfp(['CSC' num2str(tmpTet) '.ncs']);
                    if preWhiten == 1
                        cscData = prewhitening(lfpStruct.data, order); %prewhiten
                    else
                        cscData = lfpStruct.data;
                    end %whether to prewhiten

                    events = group(g).rat(r).day(d).sleep(s).popEv; %shorten variable name
                    for i = 1:length(events)
                        if isnan(events(i).r2) || events(i).r2 < r2Thresh || events(i).propClose < propCloseThresh || events(i).maxJump > jumpThresh
                            continue %to next event
                        end %r2 check

                        dayEvCntr = dayEvCntr + 1;
                        if downSampEvents == 1 && g ~= minG && ~ismember(dayEvCntr, dayEvs)
                            continue %to next events
                        end %we are discarding this event in the downsample

                        evStart = events(i).tms(1);
                        evEnd = events(i).tms(2);
                        slope = events(i).slope;
                        tmpSlopes = [tmpSlopes slope];


                        startTm = evStart - preEvTm;
                        startInd = match(startTm, lfpStruct.ts);
                        endInd = startInd + (preEvTm+postEvTm)*Fs;

                        evLFP = cscData(startInd:endInd);
                        tetLFP = cat(2, tetLFP, evLFP);

                        tmpEvLfp = lfpStruct.data (match(evStart, lfpStruct.ts):match(evEnd, lfpStruct.ts)); %NOT pre-whitened data
                        ripTFR = get_wavelet_power(tmpEvLfp, Fs, ripRange, width, useOrigNorm, dbConv); %get freq x time power for whole window
                        meanRipPow = cat(2, meanRipPow, mean(ripTFR,2));

                        gamTFR = get_wavelet_power(tmpEvLfp, Fs, gamRange, width, useOrigNorm, dbConv); %get freq x time power for whole window
                          meanGamPow = cat(2, meanGamPow, mean(gamTFR,2));

                    end %ripples - i
                end %sleep

                tetTFR = get_wavelet_power(tetLFP, Fs, fRange, width, useOrigNorm, dbConv);
                pullNorm = tetTFR(:, match(-max(normTms), tmBnds):match(-min(normTms), tmBnds),:);
                pullNorm = mean(pullNorm,2);
                zTFR = zscore(tetTFR, 0, 2);

                dayTFR = cat(3, dayTFR, zTFR);
%                  dayTFR = cat(3, dayTFR, tetTFR);

                pullEv = zTFR(:,preEvTm*Fs:(preEvTm+0.100)*Fs); %up to 100 ms after event detection
                tmpPxF = mean(pullEv, 2)';
                dayPxF = [dayPxF; tmpPxF];

                [~, pkInd] = max(tmpPxF(ripInds(1):ripInds(2)));
                tmpPkFreq = pkInd + ripRange(1) - 1;
                dayPkFreq = [dayPkFreq tmpPkFreq];

                pullRip = tetTFR(ripInds(1):ripInds(2),:);
                zRipPow = (meanRipPow - mean(pullRip,2)) ./ std(pullRip, 0, 2);

 [~, maxInd] = max(meanRipPow);
                maxFreqs = ripRange(1) + maxInd - 1;

                pkRipFreqRat{g}{r} = [pkRipFreqRat{g}{r}; maxFreqs'];

freqxSlope{g} = [freqxSlope{g}; cat(2, abs(tmpSlopes)', maxFreqs')];
                 pullGam = tetTFR(gamInds(1):gamInds(2),:);
                 zGamPow = (meanGamPow - mean(pullGam,2)) ./ std(pullGam, 0, 2);


meanPows = mean(meanGamPow,1);
                gammaPowRat{g}{r} = [gammaPowRat{g}{r}; meanPows'];

                evStat = (1:length(maxFreqs))+size(statAll,1);
                tmpStat = cat(2, repmat(g, length(maxFreqs),1), repmat(r, length(maxFreqs),1),...
                    repmat(d, length(maxFreqs),1), evStat', maxFreqs', meanPows');
                statAll = [statAll; tmpStat];
             
            end %tet

            if isempty(dayPxF)
                continue
            end %continue onto the next day

            allTFR{g} = cat(3, allTFR{g}, mean(dayTFR,3));
            PxF{g} = [PxF{g}; mean(dayPxF,1)];
            PxFRat{g}{r} = [PxFRat{g}{r};  mean(dayPxF,1)];
            pkRipFreq{g} = [pkRipFreq{g} mean(dayPkFreq)];

            for sInd = 1:2
                if ~isempty(dayTFRbyEv{sInd})
                    allTFRbyEv{g,sInd} = cat(3, allTFRbyEv{g,sInd}, mean(dayTFRbyEv{sInd},3));
                    PxFbyEv{g,sInd} = [PxFbyEv{g,sInd}; mean(dayPxFByEv{sInd},1)];
                    pkRipFreqByEv{g,sInd} = [pkRipFreqByEv{g,sInd} mean(dayPkFreqByEv{sInd})];
                end %there is any data to add
            end %sInd
            ratTFR = cat(3, ratTFR, dayTFR);
        end %day

        ratMean = mean(ratTFR, 3);
        imagesc(-preEvTm:1/lfpStruct.Fs:postEvTm, fRange(1):fRange(2), ratMean)
        axis xy
        axis square

        ylim(fRange)

        xticks(-0.3:0.1:postEvTm)
        caxis([-1 4])

        colormap(jet)
        caxis([-1.5 5])
        if r == 1
            ylabel({groupNames{g}; 'Frequency (Hz)'})
        end
        ylim(fRange)
        if g == 2
            xlabel('Time (s)')
        end
            xlim([-0.2 0.3])
        title(group(g).rat(r).name)

        TFRbyRat{g} = cat(3, TFRbyRat{g}, ratMean);
    end %rat
end %group


cd(saveDir)

cbr = colorbar;
cbr.Position = [0.92 0.18 0.0069 0.2];
ylabel(cbr, 'Power (z-score)')
same_caxis(1);

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

keyboard


%% FIG 1 - TFR BY RAT

tmpFunc = @(x)(mean(x,3));
meanByGroup = cellfun(tmpFunc, TFRbyRat, 'UniformOutput', false);

figtitle = 'TFR_byRat';
if downSampEvents == 1
    figtitle = [figtitle '_downSampEvents'];
end %down samp
% if preWhiten == 1
%     figtitle = [figtitle '_preWhiten'];
% end %add prewhiten to title
figure('Name', figtitle, 'Position', [722 174 518 804])

for g = 1:2
    subplot(2,1,g)

    imagesc(-preEvTm:1/Fs:postEvTm, fRange(1):fRange(2), meanByGroup{g})
    axis xy
    axis square
    colormap(jet)

    ylabel('Frequency (Hz)')
    ylim(fRange)
    xlabel('Time (s)')
    xticks(-0.3:0.1:postEvTm)
    title(group(g).name)
    caxis([-1 3])

    xlim([-0.2 0.3])
end %group

same_caxis(gcf, 1)
cbr = colorbar;
cbr.Position = [0.85 0.58 0.0261 0.35];
ylabel(cbr, 'z-scored power')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% PxF BY RAT

figtitle = 'PxF_allRats';
figure('Name', figtitle, 'Position', [23 551 1849 417])

for g = 1:2
    for r = 1:6
        subplot(2, 6, (g-1)*6+r)
        title(group(g).rat(r).name)
        if isempty(PxFRat{g}{r})
            continue
        end
        pullData = PxFRat{g}{r};
        meanData = mean(pullData,1);
        hold on
        plot(fRange(1):fRange(2), meanData, 'Color', rgb(cols{g}))

        xlim([fRange(1) fRange(2)]);
        %         ylim([-0.5 4])
        if g == 2
            xlabel('Frequency (Hz)')
        end
        if r == 1
            ylabel('z-scored power')
        end

        [~, pkInd] = max(meanData(ripInds(1):ripInds(2)));
        pkFreq = fRangeCont(pkInd+ripInds(1)-1);
        line([pkFreq pkFreq], [-0.5 4], 'LineStyle', '--', 'Color', rgb(cols{g}))
        text(10, 3.75, ['Peak freq = ' num2str(pkFreq) ' Hz'])

    end %rat
end %group

same_axes;

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIG 3 - PxF BY RAT

figtitle = 'PxF_byRat';
if downSampEvents == 1
    figtitle = [figtitle '_downSampEvents'];
end %down samp
% if preWhiten == 1
%     figtitle = [figtitle '_preWhiten'];
% end %add prewhiten to title
figure('Name', figtitle, 'Position', [194 558 1204 420])

leg = cell(2,1);
for g = 1:2
    tmpFunc = @(x)(mean(x,1));
    tmpMean = cellfun(tmpFunc, PxFRat{g}, 'UniformOutput', false)';
    tmpGMean = vertcat(tmpMean{:});

    meanData = mean(tmpGMean,1);
    CI = bootci(nboot, {tmpFunc, tmpGMean}, 'type', 'per');

    lh(g) =  plot_filled_ci(fRange(1):fRange(2), meanData, CI, rgb(cols{g}));
    leg{g} = [groupNames{g} ' (n = ' num2str(size(tmpGMean,1)) ' rats)'];
end %group

xlim([fRange(1) fRange(2)]);
xlabel('Frequency (Hz)')
ylabel('z-scored power')
yticks([0 1:5])
legend(lh, leg, 'Location', 'southeast')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIG 3 - PxF BY GROUP - ALL

figtitle = 'PxF_byDay';
if downSampEvents == 1
    figtitle = [figtitle '_downSampEvents'];
end %down samp
% if preWhiten == 1
%     figtitle = [figtitle '_preWhiten'];
% end %add prewhiten to title
figure('Name', figtitle, 'Position', [194 558 1204 420])

for g = 1:2

    pullData = PxF{g};
    meanData = mean(pullData,1);
    tmpFunc = @(x)(mean(x,1));
    CI = bootci(nboot, {tmpFunc, pullData}, 'type', 'per');

    lh(g) =  plot_filled_ci(fRange(1):fRange(2), meanData, CI, rgb(cols{g}));
end %group

xlim([fRange(1) fRange(2)]);
xlabel('Frequency (Hz)')
ylabel('z-scored power')
legend(lh, {group(1).name, group(2).name}, 'Location', 'southeast')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIG 5 - PEAK RIP FREQ SORTED

figtitle = ('PeakRipFreq_byRat_CI_sorted');
figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = cell(8,1);
for g = 1:2
    tmpMeans = cellfun(@mean, pkRipFreqRat{g});
    [~, sortOrd] = sort(tmpMeans);
%     sortOrd = sortOrd - sum(isnan(tmpMeans));
    for r = 1:length(pkRipFreqRat{g})
        if isempty(pkRipFreqRat{g}{r})
            continue
        end
        hold on
%         ratCntr = ratCntr + 1;
xVal = find (r==sortOrd) + 4*(g-1);
        h = dotplot(xVal, {pkRipFreqRat{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([xVal-0.5 xVal-0.5], [0 350], 'Color', 'k')
        line([xVal+0.5 xVal+0.5], [0 350], 'Color', 'k')
        tickLabs{xVal} = group(g).rat(r).name;

          quants = quantile(pkRipFreqRat{g}{r}, 3);
            tmpCI = bootci(nboot, {@mean, pkRipFreqRat{g}{r}}, 'type', 'per');
    meanData = mean(pkRipFreqRat{g}{r});
%         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
pt = patch([xVal-0.25 xVal-0.25 xVal+0.25 xVal+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
    line([xVal-0.25 xVal+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 9])
xticks(1:8)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([125 275])
yticks(100:50:350)
ylabel('Peak ripple frequency (Hz)')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end


%% FIG 5 - GAMMA POWER SORTED

figtitle = ('GammaPower_byRat_CI_sorted');
figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = cell(8,1);
for g = 1:2
    tmpMeans = cellfun(@mean, gammaPowRat{g});
    [~, sortOrd] = sort(tmpMeans);
%     sortOrd = sortOrd - sum(isnan(tmpMeans));
    for r = 1:length(gammaPowRat{g})
        if isempty(gammaPowRat{g}{r})
            continue
        end
        hold on
%         ratCntr = ratCntr + 1;
xVal = find (r==sortOrd) + 4*(g-1);
        h = dotplot(xVal, {gammaPowRat{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([xVal-0.5 xVal-0.5], [-15 350], 'Color', 'k')
        line([xVal+0.5 xVal+0.5], [-15 350], 'Color', 'k')
        tickLabs{xVal} = group(g).rat(r).name;

          quants = quantile(gammaPowRat{g}{r}, 3);
            tmpCI = bootci(nboot, {@mean, gammaPowRat{g}{r}}, 'type', 'per');
    meanData = mean(gammaPowRat{g}{r});
%         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
pt = patch([xVal-0.25 xVal-0.25 xVal+0.25 xVal+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
    line([xVal-0.25 xVal+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 9])
xticks(1:8)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0 50])
yticks(-10:10:40)
ylabel('Gamma power (dB)')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end


end %function