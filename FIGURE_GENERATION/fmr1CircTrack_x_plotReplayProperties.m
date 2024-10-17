function fmr1CircTrack_x_plotReplayProperties(group)
% function fmr1CircTrack_x_plotReplayProperties(group)
%
% PURPOSE:
%   Plots the properties of replay events between genotypes.
%
% MMD
% 12/2021 - but based on code originally written by MMD 06/2020
% Colgin Lab

%% OPTIONS

saveOrNot = 1; %to save figs

downSampEvents = 0; %to randomly down sample events to to same as lowest group

saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\REPLAY\replayProperties';
cd(saveDir)

cols = {'Blue', 'Red'};
groupNames = {'WT', 'FXS'};

r2Thresh = 0.5;
propCloseThresh = 0;
jumpThresh = 0.25*2*pi; %one quarter of the track

%% INITIALIZE

durAll = cell(2,1);
durRat = cell(2,1);
slopesAll = cell(2,1); %absolute values of all slopes
slopesRat = cell(2,1); %by rat
fitDistAll = cell(2,1); %rescaled so time of SWR event is 0-1
avgJumpAll = cell(2,1); % average spatial jump for decoded positions in replay event
pathDistAll = cell(2,1); %just dist between end and beginning of decoded path
pathDistRat = cell(2,1);
rateRat = cell(2,1);
rateRat(:) = {cell(1,6)};
r2Rat = cell(2,1);
r2Rat(:) = {cell(1,6)};

radBinCtrs = group(2).rat(1).day(1).binCtrs; %doesn't change across days/rats
radBinCtrs = deg2rad(radBinCtrs);

radCmConv = (pi*100) / (2*pi); %convert all track measurements from deg to cm before plotting

statAll = [];
statRate = [];
ratCntr = 0;

%% PREP FOR DOWN SAMPLE

if downSampEvents == 1
    rsE = RandStream('mt19937ar');

    gNum = zeros(1,2);
    for g = 1:2
        for r = 1:length(group(g).rat)
            for d = 1:length(group(g).rat(r).day)
                for s = 2:5
                    if ~isempty(group(g).rat(r).day(d).sleep(s).unit)
                        for i = 1:length(group(g).rat(r).day(d).sleep(s).popEv)
                            if group(g).rat(r).day(d).sleep(s).popEv(i).r2 >= r2Thresh
                                gNum(g) = gNum(g) + 1;
                            end %meets threshold
                        end %i
                    end %if there are events
                end %sleep
            end %day
        end %rat
    end %group

    [minEv, minG] = min(gNum);
    keepEvs = datasample(rsE, 1:max(gNum), minEv, 'Replace', false);
    keepEvs = sort(keepEvs);

    allEvCntr = 0; %initialize

end %if down sampling events


%% GET DATA

for g = 1:2
    for r = 1:length(group(g).rat)
        ratCntr = ratCntr + 1;
        ratSlope = [];
        ratDur  = [];
        ratDist = [];
        for d = 1:length(group(g).rat(r).day)

            numEvents = 0;
            for s = 2:5 %getting getting total number of events for day
                if ~isempty(group(g).rat(r).day(d).sleep(s).unit)

                    for i = 1:length(group(g).rat(r).day(d).sleep(s).popEv)
                        if group(g).rat(r).day(d).sleep(s).popEv(i).r2 >= r2Thresh
                            numEvents = numEvents + 1;
                        end %meets threshold
                    end %i
                end %if there is anything in popEvents
            end %sleep

            if numEvents == 0
                continue %to next day
            end %if there are no events

            if downSampEvents == 1 && g ~= minG %find out which events from this day to keep
                dayEvs = 1:numEvents; %initialize
                keepInds = find(keepEvs > allEvCntr & keepEvs <= allEvCntr+numEvents); %which from keepEvs are in this day (no = in first since allEvCntr starts at 0)
                dayEvInds = keepEvs(keepInds) - allEvCntr;
                dayEvs = dayEvs(dayEvInds); %the events from this day that get to be included

                allEvCntr = allEvCntr + numEvents;
            end %down sample events

            dayEvCntr = 0;

            for s = 2:length(group(g).rat(r).day(d).sleep) %after experience, so just sleep 2 on
                if isempty(group(g).rat(r).day(d).sleep(s).unit) %if there is anything in the sleep
                    continue % to next sleep
                end %if empty
                events = group(g).rat(r).day(d).sleep(s).popEv; %shorten variable name

                tmpRate = length(events) ./ (group(g).rat(r).day(d).sleep(s).coords(end,1)-group(g).rat(r).day(d).sleep(s).coords(1,1));
                rateRat{g}{r} = [rateRat{g}{r}; tmpRate];
                tmpRateStat = [g ratCntr d s tmpRate];
                statRate = [statRate; tmpRateStat];

                for i = 1:length(events)

                    r2 = events(i).r2;
                    if isnan(events(i).r2) || events(i).r2 < r2Thresh || events(i).propClose < propCloseThresh || events(i).maxJump > jumpThresh %it doesn't meet threshold
                        continue %to next event
                    end %r2 thresh
                    dayEvCntr = dayEvCntr + 1;
                    if downSampEvents == 1 && g ~= minG && ~ismember(dayEvCntr, dayEvs)
                        continue %to next events
                    end %we are discarding this event in the downsample

                    r2Rat{g}{r} = [r2Rat{g}{r}; r2];

                    startTm = events(i).tms(1);
                    endTm = events(i).tms(2);

                    slope = events(i).slope;
                    if slope > 0 %forward event
                        sInd = 1; %positive slope
                    else %negative event
                        sInd = 2;
                    end %which direction is slope
                    pxn = events(i).pxn;
                    eventDur = endTm-startTm;

                    calphase = events(i).calphase;
                    calDist = circ_dist(calphase(end), calphase(1));
                    if sInd == 1 && calDist < 0  %circ_dist just takes the shortest abs val amount, + or -
                        calDist = 2*pi + calDist;
                    elseif sInd == 2 && calDist > 0
                        calDist = 2*pi - calDist;
                    end %check that distance matches slope direction

                    durAll{g} = [durAll{g} eventDur];
                    ratDur = [ratDur eventDur];
                    slopesAll{g} = [slopesAll{g} abs(slope)*radCmConv];
                    ratSlope = [ratSlope abs(slope)*radCmConv];
                    fitDistAll{g} = [fitDistAll{g} abs(calDist)*radCmConv];

                    %                     if abs(slope)*radCmConv > 3500 && g == 2 && r == 6
                    %                         keyboard
                    %                     end

                    com = events(i).com;
                    tmpPathDist = circ_dist(com(end), com(1)); %distance between two points
                    if sInd == 1 && tmpPathDist < 0  %circ_dist just takes the shortest abs val amount, + or -
                        tmpPathDist = 2*pi + tmpPathDist;
                    elseif sInd == 2 && tmpPathDist > 0
                        tmpPathDist = 2*pi - tmpPathDist;
                    end %check that distance matches slope direction
                    pathDistAll{g} = [pathDistAll{g} abs(tmpPathDist)*radCmConv];
                    ratDist = [ratDist; abs(tmpPathDist)*radCmConv];

                    jumpDist = nan(1,size(pxn,2)-1);
                    for ji = 1:length(jumpDist)
                        if ismember(ji, events(i).bins2use) && ismember(ji+1, events(i).bins2use)
                            comInd = find(events(i).bins2use==ji);
                            tmpJump = circ_dist(com(comInd+1), com(comInd));
                            jumpDist(ji) = abs(tmpJump);
                        end
                    end %jumpInd
                    avgJumpAll{g} = [avgJumpAll{g} nanmean(jumpDist)*radCmConv];

                    tmpStat = [g ratCntr size(statAll,1)+1 d s eventDur abs(slope)*radCmConv abs(calDist)*radCmConv abs(tmpPathDist)*radCmConv r2 i]; %group rat event sleep dur slope fitDist pathDist
                    statAll = [statAll; tmpStat];

                end %i - popevents


            end %sleep

        end %day
        slopesRat{g} = cat(1, slopesRat{g}, {ratSlope});
        durRat{g} = cat(1, durRat{g}, {ratDur});
        pathDistRat{g} = cat(1, pathDistRat{g}, {ratDist});
    end %rat
end %group
keyboard

%% FIGS - DUR BY RAT CI SORTED

figtitle = ('Dur_byRat_CI_sorted');
if downSampEvents == 1
    figtitle = [figtitle '_downSampEvents'];
end %down samp

figure('Name', figtitle, 'Position',  [680 558 448 420])

tmpJit = 0.25;
nboot = 1000;
tickLabs = cell(8,1);
for g = 1:2
    tmpMeans = cellfun(@mean, durRat{g});
    [~, sortOrd] = sort(tmpMeans);
    %     sortOrd = sortOrd - sum(isnan(tmpMeans));
    for r = 1:length(durRat{g})
        if isempty(durRat{g}{r})
            continue
        end
        hold on
        %         ratCntr = ratCntr + 1;
        xVal = find (r==sortOrd) + 4*(g-1);
        h = dotplot(xVal, {durRat{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([xVal-0.5 xVal-0.5], [0 0.75], 'Color', 'k')
        line([xVal+0.5 xVal+0.5], [0 0.75], 'Color', 'k')
        tickLabs{xVal} = group(g).rat(r).name;

        quants = quantile(durRat{g}{r}, 3);
        tmpCI = bootci(nboot, {@mean, durRat{g}{r}}, 'type', 'per');
        meanData = mean(durRat{g}{r});
        %         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
        pt = patch([xVal-0.25 xVal-0.25 xVal+0.25 xVal+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
        line([xVal-0.25 xVal+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 9])
xticks(1:8)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0 0.5])
ylabel('Duration (s)')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end

%% FIGS

figtitle = ('Slope_byRat_CI_sorted');
if downSampEvents == 1
    figtitle = [figtitle '_downSampEvents'];
end %down samp

figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = cell(8,1);
for g = 1:2
    tmpMeans = cellfun(@mean, slopesRat{g});
    [~, sortOrd] = sort(tmpMeans);
    %     sortOrd = sortOrd - sum(isnan(tmpMeans));
    for r = 1:length(slopesRat{g})
        if isempty(slopesRat{g}{r})
            continue
        end
        hold on
        %         ratCntr = ratCntr + 1;
        xVal = find (r==sortOrd) + 4*(g-1);
        h = dotplot(xVal, {slopesRat{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([xVal-0.5 xVal-0.5], [0 5000], 'Color', 'k')
        line([xVal+0.5 xVal+0.5], [0 5000], 'Color', 'k')
        tickLabs{xVal} = group(g).rat(r).name;

        quants = quantile(slopesRat{g}{r}, 3);
        tmpCI = bootci(nboot, {@mean, slopesRat{g}{r}}, 'type', 'per');
        meanData = mean(slopesRat{g}{r});
        %         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
        pt = patch([xVal-0.25 xVal-0.25 xVal+0.25 xVal+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
        line([xVal-0.25 xVal+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 9])
xticks(1:8)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0 5000])
yticks(0:1000:5000)
ylabel('Slope (cm/s)')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end


%% FIGS - R2 BY RAT SORTED

figtitle = ('r2_byRat_CI_sorted');
if downSampEvents == 1
    figtitle = [figtitle '_downSampEvents'];
end %down samp

figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = cell(8,1);
for g = 1:2
    tmpMeans = cellfun(@mean, r2Rat{g});
    [~, sortOrd] = sort(tmpMeans);
    %     sortOrd = sortOrd - sum(isnan(tmpMeans));
    for r = 1:length(r2Rat{g})
        if isempty(r2Rat{g}{r})
            continue
        end
        hold on
        %         ratCntr = ratCntr + 1;
        xVal = find (r==sortOrd) + 4*(g-1);
        h = dotplot(xVal, {r2Rat{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([xVal-0.5 xVal-0.5], [0 5000], 'Color', 'k')
        line([xVal+0.5 xVal+0.5], [0 5000], 'Color', 'k')
        tickLabs{xVal} = group(g).rat(r).name;

        quants = quantile(r2Rat{g}{r}, 3);
        tmpCI = bootci(nboot, {@mean, r2Rat{g}{r}}, 'type', 'per');
        meanData = mean(r2Rat{g}{r});
        %         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
        pt = patch([xVal-0.25 xVal-0.25 xVal+0.25 xVal+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
        line([xVal-0.25 xVal+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 9])
xticks(1:8)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0.45 1])
yticks(0:0.1:1)
ylabel('r^2')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end


%% FIGS - PATH DIST BY RAT CI SORTED

figtitle = ('PathDist_byRat_CI_sorted');
if downSampEvents == 1
    figtitle = [figtitle '_downSampEvents'];
end %down samp

figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = cell(8,1);
for g = 1:2
    tmpMeans = cellfun(@mean, pathDistRat{g});
    [~, sortOrd] = sort(tmpMeans);
    %     sortOrd = sortOrd - sum(isnan(tmpMeans));
    for r = 1:length(pathDistRat{g})
        if isempty(pathDistRat{g}{r})
            continue
        end
        hold on
        %         ratCntr = ratCntr + 1;
        xVal = find (r==sortOrd) + 4*(g-1);
        h = dotplot(xVal, {pathDistRat{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([xVal-0.5 xVal-0.5], [0 400], 'Color', 'k')
        line([xVal+0.5 xVal+0.5], [0 400], 'Color', 'k')
        tickLabs{xVal} = group(g).rat(r).name;

        quants = quantile(pathDistRat{g}{r}, 3);
        tmpCI = bootci(nboot, {@mean, pathDistRat{g}{r}}, 'type', 'per');
        meanData = mean(pathDistRat{g}{r});
        %         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
        pt = patch([xVal-0.25 xVal-0.25 xVal+0.25 xVal+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
        line([xVal-0.25 xVal+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 9])
xticks(1:8)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0 350])
ylabel('Path distance (cm)')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end

%% FIGS - RATE BY RAT CI SORTED

figtitle = ('Rate_byRat_CI_sorted');
if downSampEvents == 1
    figtitle = [figtitle '_downSampEvents'];
end %down samp

figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = cell(8,1);
for g = 1:2
    tmpMeans = cellfun(@mean, rateRat{g});
    [~, sortOrd] = sort(tmpMeans);
    %     sortOrd = sortOrd - sum(isnan(tmpMeans));
    for r = 1:length(rateRat{g})
        if isempty(rateRat{g}{r})
            continue
        end
        hold on
        %         ratCntr = ratCntr + 1;
        xVal = find (r==sortOrd) + 4*(g-1);
        h = dotplot(xVal, {rateRat{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([xVal-0.5 xVal-0.5], [0 0.75], 'Color', 'k')
        line([xVal+0.5 xVal+0.5], [0 0.75], 'Color', 'k')
        tickLabs{xVal} = group(g).rat(r).name;

        quants = quantile(rateRat{g}{r}, 3);
        tmpCI = bootci(nboot, {@mean, rateRat{g}{r}}, 'type', 'per');
        meanData = mean(rateRat{g}{r});
        %         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
        pt = patch([xVal-0.25 xVal-0.25 xVal+0.25 xVal+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
        line([xVal-0.25 xVal+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 9])
xticks(1:8)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0 0.5])
ylabel('Population event rate (Hz)')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end


end %function