function fmr1CircTrack_x_replayEventISI(group)
% function fmr1CircTrack_x_replayEventISI(group)
%
% PURPOSE:
%   Plot the inter-spike interval during replay events in the WT and KO
%   groups. This function will plot both the ISI of all spikes in the event
%   and the first spike from each cells that fires in the event.
%
%
% MMD
% 12/2021
% Colgin Lab

%% OPTIONS

saveOrNot = 1; %to save figs

detectMethod = 1; %1 = from population events, 2 = from SWRs in the LFP

prepForStats = 1; %1 to prepare output for SPSS and keyboard in code

downSampEvents = 0; %so number of events equal between groups
rsE = RandStream('mt19937ar');

%% INITIALIZE

saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\REPLAY\ISI';
curDir = pwd;

popISI = cell(2,1); %by group
popISIrat = cell(2,1);
popISIrat(:) = {cell(1,6)};
firstSpkPopISI = cell(2,1); %only between the first spike from each cell in the event
fsISIrat = cell(2,1);
fsISIrat(:) = {cell(1,6)};

percParU = cell(2,1);

allISIs = cell(2,1);
firstSpkAllISI = cell(2,1);

nCell = cell(2,1); %controls
nSpks = cell(2,1);

allDurs = cell(2,1);
r2Thresh = 0.5;
maxJumpDist = 0.25*2*pi; %rad
propCloseThresh = 0;

cols = {'Blue', 'Red'};
groupNames = {'WT', 'FXS'};
markers = {'o', 'square', '*', '+', 'x', 'diamond', '^', '+'};
fCols = cat(1, [1 1 1], rgb('Blue'), [1 1 1], [1 1 1], [1 1 1], rgb('Red'), [1 1 1], [1 1 1]);
alpha = 1;

if prepForStats == 1
    statAll = [];
end %stats

ratCntr = 0;

%% GET DATA

for g = 1:2
    for r = 1:length(group(g).rat)
        ratCntr = ratCntr + 1;
        for d = 1:length(group(g).rat(r).day)
            if group(g).rat(r).day(d).decThresh == 0
                continue
            end %not enough to decode

            for s = 2:5 %start from 2, after experience
                if isempty(group(g).rat(r).day(d).sleep(s).unit) %if there is anything in the sleep
                    continue %to next event
                end %if there's any data
                if detectMethod == 1
                    events = group(g).rat(r).day(d).sleep(s).popEv;
                else
                    events = group(g).rat(r).day(d).sleep(s).rip;
                end %which method

                load([group(g).rat(r).day(d).sleep(s).dir '\PopEventOut_allCells.mat'])

                for i = 1:length(events)
                    if events(i).r2 < r2Thresh || isnan(events(i).r2) || events(i).propClose < propCloseThresh || events(i).maxJump > maxJumpDist %it doesn't meet threshold
                        continue
                    end %doesn't meet threshold

                    startTm = events(i).tms(1);
                    endTm = events(i).tms(2);
                    evMat = timeMap(events(i).inds(1):events(i).inds(2));

                    allDurs{g} = [allDurs{g} endTm-startTm];

                    allSpkTms = [];
                    allFirstSpkTms = [];

                    partU = events(i).nCell; %number of participating units
                    nCell{g} = [nCell{g}; events(i).nCell];
%                     if events(i).nCell > 30
%                         keyboard
%                     end
                    %                     [~, uInds] = ismember(events(i).actCell, uIDsAll, 'row'); %get inds of cells that fired in this event
                    for uInd = 1:partU
                        [~, u] = ismember(events(i).actCell(uInd,:), uIDs, 'row');
                        spkInds = full(spkRstr(u, events(i).inds(1):events(i).inds(2)));
                        try
                            evSpks = evMat(find(spkInds))';
                        catch; keyboard; end
                        allSpkTms = [allSpkTms; evSpks];
                        try
                        allFirstSpkTms = [allFirstSpkTms; evSpks(1)];
                        catch; keyboard; end
                    end %units active in this event

                    nSpks{g} = [nSpks{g}; length(allSpkTms)];
                    sortAllSpks = sort(allSpkTms); %ascending order
                    spksISI = diff(sortAllSpks); %difference between spike times - ISI
                    meanISI = mean(spksISI); %get mean for this event
                    popISI{g} = [popISI{g}; meanISI];
                    popISIrat{g}{r} = [popISIrat{g}{r}; meanISI];

                    sortFirstSpks = sort(allFirstSpkTms); %ascending order
                    firstSpksISI = diff(sortFirstSpks); %difference between spike times - ISI
                    meanFirstISI = mean(firstSpksISI); %get mean for this event
                    if isnan(meanFirstISI)
                        keyboard
                    end
                    firstSpkPopISI{g} = [firstSpkPopISI{g}; meanFirstISI];
                    fsISIrat{g}{r} = [fsISIrat{g}{r}; meanFirstISI];

                    percParU{g} = [percParU{g}; partU/length(uIDs)*100];

                    allISIs{g} = [allISIs{g}; spksISI];
                    firstSpkAllISI{g} = [firstSpkAllISI{g}; firstSpksISI];

                    if prepForStats == 1
                        % g ratCntr d s evNum (data point)
                        tmpEvID = size(statAll,1) + 1;

                        if abs(events(i).slope) > 0 %forward event
                            sInd = 1; %positive slope
                        else %negative event
                            sInd = 2;
                        end %which direction is slope

                        lineStart = [g ratCntr d s sInd tmpEvID];

                        statAll = [statAll; lineStart meanISI meanFirstISI];
                    end %stats
                end %i - popEvents
            end %sleep
        end %day
    end %rat
end %group

cd(saveDir)

keyboard

%% FIGS  POP ISI BY RAT CI SORTED

figtitle = ('PopISI_byRat_CI_sorted');
if downSampEvents == 1
    figtitle = [figtitle '_downSampEvents'];
end %down samp

figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = cell(8,1);
for g = 1:2
    tmpMeans = cellfun(@mean, popISIrat{g});
    [~, sortOrd] = sort(tmpMeans);
%     sortOrd = sortOrd - sum(isnan(tmpMeans));
    for r = 1:length(popISIrat{g})
        if isempty(popISIrat{g}{r})
            continue
        end
        hold on
%         ratCntr = ratCntr + 1;
xVal = find (r==sortOrd) + 4*(g-1);
        h = dotplot(xVal, {popISIrat{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([xVal-0.5 xVal-0.5], [0 5000], 'Color', 'k')
        line([xVal+0.5 xVal+0.5], [0 5000], 'Color', 'k')
        tickLabs{xVal} = group(g).rat(r).name;

          quants = quantile(popISIrat{g}{r}, 3);
            tmpCI = bootci(nboot, {@mean, popISIrat{g}{r}}, 'type', 'per');
    meanData = mean(popISIrat{g}{r});
%         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
pt = patch([xVal-0.25 xVal-0.25 xVal+0.25 xVal+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
    line([xVal-0.25 xVal+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 9])
xticks(1:8)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0 0.02])
yticks(0:0.005:0.02)
ylabel('Population ISI (s)')


if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
        saveas(gcf, figtitle, 'epsc');
end


%% FIGS  FIRST SPIKE ISI BY RAT CI SORTED

figtitle = ('FirstSpikeISI_byRat_CI_sorted');

figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = cell(8,1);
for g = 1:2
    tmpMeans = cellfun(@mean, fsISIrat{g});
    [~, sortOrd] = sort(tmpMeans);
%     sortOrd = sortOrd - sum(isnan(tmpMeans));
    for r = 1:length(fsISIrat{g})
        if isempty(fsISIrat{g}{r})
            continue
        end
        hold on
%         ratCntr = ratCntr + 1;
xVal = find (r==sortOrd) + 4*(g-1);
        h = dotplot(xVal, {fsISIrat{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([xVal-0.5 xVal-0.5], [0 5000], 'Color', 'k')
        line([xVal+0.5 xVal+0.5], [0 5000], 'Color', 'k')
        tickLabs{xVal} = group(g).rat(r).name;

          quants = quantile(fsISIrat{g}{r}, 3);
            tmpCI = bootci(nboot, {@mean, fsISIrat{g}{r}}, 'type', 'per');
    meanData = mean(fsISIrat{g}{r});
%         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
pt = patch([xVal-0.25 xVal-0.25 xVal+0.25 xVal+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
    line([xVal-0.25 xVal+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 9])
xticks(1:8)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0 0.04])
yticks(0:0.01:0.04)
ylabel('First spike ISI (s)')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
        saveas(gcf, figtitle, 'epsc');
end


cd(curDir)

end %function