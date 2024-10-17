function fmr1CircTrack_x_plotPlaceCellProperties(group)
% function fmr1CircTrack_x_plotPlaceCellProp(group)
%
% PURPOSE:
%   Plot place cell firing properties as in Figure 3 of Mably et al. 2016
%   for data collected for the WT/FXS circle track data.
%
% INPUT:
%   group data struct
%
% OUTPUT:
%   F1: Spatial correlation.
%   F2: Rate overlap.
%   F3: Spatial information.
%   F4: Spatial sparsity.
%   F5: Place field size.
%   F6: Place fields per cell.
%   F7: Mean firing rate across begin sessions.
%   F8: Peak in-field firing rate.
%   F9: Average in-field firing rate.
%
% MM Donahue
% 5/2021
% Colgin Lab

%% OPTIONS

saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\PLACE_CELLS\placeCellProperties';
saveOrNot = 1;

prepForStats = 0;

cols = {'Blue', 'Red'};
alpha = 0.5;
groupNames = {'WT', 'FXS'};

spatBinSz = 4;

minFr = 1; %for units to be included

%% INITIALIZE

% Include all CA1 units, as in Mably et al. 2016 (Fig 3)
spatCorr = cell(2,5); %group x combo
spatCorr(:) = {cell(6,1)};
rateOverlap = cell(2,5); %same
rateOverlap(:) = {cell(6,1)};

combos = [1 2; 2 3; 3 4; 1 4]; %combos as listed above
combosTxt = {'1-2', '2-3', '3-4', '1-4'};

% Cells with identified place fields
spatInfo = cell(2,1); %by group
spatInfo(:) = {cell(6,1)};
sparsity = cell(2,1);
sparsity(:) = {cell(6,1)};
meanFirRate = cell(2,1);
meanFirRate(:) = {cell(6,1)};
peakFieldFirRate = cell(2,1);
peakFieldFirRate(:) = {cell(6,1)};
meanFieldFirRate = cell(2,1);
meanFieldFirRate(:) = {cell(6,1)};
fieldSize = cell(2,1);
fieldSize(:) = {cell(6,1)};
fieldPerCell = cell(2,1);
fieldPerCell(:) = {cell(6,1)};

convFact = (2*pi*50)/(360/spatBinSz); %cm per bin, because track is 100 cm in diameter and ratemap has 72 bins

uCntr = 0; %initialize for stats (ANOVAs)

statStab = [];
statCell = [];
statField = [];

%% GET DATA

for g = 1:2
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
%             if group(g).rat(r).day(d).decThresh == 0 
%                 continue
%             end
            for u = 1:length(group(g).rat(r).day(d).xBegUnitInfo)
                if group(g).rat(r).day(d).xBegUnitInfo(u).meetMin == 0
                    continue %to next unit
                end %doens't reach max firing rate
                uCntr = uCntr + 1;
                uRateMaps = zeros(4,360/spatBinSz); %for storing rate maps from all begins for this unit

                for b = 1:length(group(g).rat(r).day(d).begin)
                    uRateMaps(b,:) = group(g).rat(r).day(d).begin(b).unit(u).smRateMap; %get all of the ratemaps for the spat corr comparisons
                end %begin

                for c = 1:length(combosTxt)
                    b1 = combos(c,1); %begins to compare
                    b2 = combos(c,2);

                    rm1 = uRateMaps(b1,:); %ratemaps to compare
                    rm2 = uRateMaps(b2,:);

                    scMat = corrcoef(rm1,rm2);
                    if ~isnan(scMat(2))
                        spatCorr{g,c}{r} = [spatCorr{g,c}{r} scMat(2)];
                    end

                    fr1 = mean(uRateMaps(b1,:));
                    fr2 = mean(uRateMaps(b2,:));

                    rateRatio = min([fr1 fr2]) / max([fr1 fr2]); % Calculate the rate ratio with the lower FR as numerator
                    if ~isnan(rateRatio)
                        rateOverlap{g,c}{r} = [rateOverlap{g,c}{r} rateRatio];
                    end

                    tmpStat = [g r d uCntr c scMat(2) rateRatio];
                    statStab = [statStab; tmpStat];
                end %combos
                
                if ~isempty(group(g).rat(r).day(d).xBegUnitInfo(u).pf) %only cells with place fields from avg ratemap

   rateMap = group(g).rat(r).day(d).xBegUnitInfo(u).smRateMap; %get AVG ratemap across all begins
                    timePerBin = group(g).rat(r).day(d).xAllBeginTpb; %same, all begins

                    tmpSpatInfo = get_spatial_info(rateMap, timePerBin);
                    spatInfo{g}{r} = [spatInfo{g}{r} tmpSpatInfo];

                    tmpSparse = get_spatial_sparsity(rateMap, timePerBin);
                    sparsity{g}{r} = [sparsity{g}{r} tmpSparse];

                    fieldPerCell{g}{r} = [fieldPerCell{g}{r} length(group(g).rat(r).day(d).xBegUnitInfo(u).pf)];

                    meanFirRate{g}{r} = [meanFirRate{g}{r} mean(uRateMaps(:))];

                    tmpStat = [g r d uCntr tmpSpatInfo tmpSparse length(group(g).rat(r).day(d).xBegUnitInfo(u).pf)  mean(uRateMaps(:))];
                    statCell = [statCell; tmpStat];
              

                    for pf = 1:length(group(g).rat(r).day(d).xBegUnitInfo(u).pf)
                        peakFieldFirRate{g}{r} = [ peakFieldFirRate{g}{r} group(g).rat(r).day(d).xBegUnitInfo(u).pf(pf).pkFr];

                        pfInds = group(g).rat(r).day(d).xBegUnitInfo(u).pf(pf).inds;
                        pullBins = zeros(1,length(pfInds)); %initialize so I know it works right

                        if pfInds(1) < pfInds(end) %it doesn't wrap around 0 degrees
                            pullBins = rateMap(1,pfInds(1):pfInds(end));
                        else %pf does wrap around 0 degrees
                            findWrap = find(pfInds == 1); %find where it crosses 0 degree boundary
                            pullBins(1,1:findWrap-1) = rateMap(1,pfInds(1):pfInds(findWrap-1));
                            pullBins(1,findWrap:end) = rateMap(1,pfInds(findWrap):pfInds(end));
                        end %whether or not place field wraps around 0

                        meanFieldFirRate{g}{r}  = [meanFieldFirRate{g}{r}  mean(pullBins)];
                        fieldSize{g}{r} = [fieldSize{g}{r} length(pfInds)*convFact];

                        tmpStat = [g r d size(statField,1)+1 group(g).rat(r).day(d).xBegUnitInfo(u).pf(pf).pkFr  mean(pullBins)  length(pfInds)*convFact];
                        statField = [statField; tmpStat];
                    end %pf
                end %if this unit has a place field
            end %unit
        end %day
    end %rat
end %group

cd(saveDir)
keyboard

%% SPATIAL CORRELATION

figtitle = 'SpatialCorr_byRat_CI';
figure('Name', figtitle, 'Position',  [233 479 1448 420])

for c = 1:length(combosTxt)
    subplot(1,length(combosTxt),c)
    ratCntr = 0;
    tmpJit = 0.25;
    nboot = 1000;
    tickLabs = {};
    for g = 1:2
        for r = 1:length(spatCorr{g,c})
            if isempty(spatCorr{g,c}{r})
                continue
            end
            hold on
            ratCntr = ratCntr + 1;
            h = dotplot(ratCntr, {spatCorr{g,c}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

            line([ratCntr-0.5 ratCntr-0.5], [-0.5 1], 'Color', 'k')
            line([ratCntr+0.5 ratCntr+0.5], [-0.5 1], 'Color', 'k')
            tickLabs = cat(1, tickLabs, group(g).rat(r).name);

            quants = quantile(spatCorr{g,c}{r}, 3);
            tmpCI = bootci(nboot, {@mean, spatCorr{g,c}{r}}, 'type', 'per');
            meanData = mean(spatCorr{g,c}{r});
            %         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
            pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
            line([ratCntr-0.25 ratCntr+0.25], [meanData meanData], 'Color', 'k')
        end %rat
    end %group

    title(combosTxt{c})
    axis square
    xlim([0 13])
    xticks(1:12)
    xticklabels(tickLabs)
    xlabel('Rat #')
    ylim([0 1])
    yticks(0:0.5:1)
    if c == 1
        ylabel('Spatial correlation')
    end
end %combo

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end


%% RATE OVERLAP

figtitle = 'RateOverlap_byRat_CI';
figure('Name', figtitle, 'Position',  [233 479 1448 420])

for c = 1:4
    subplot(1,4,c)
    ratCntr = 0;
    tmpJit = 0.25;
    nboot = 1000;
    tickLabs = {};
    for g = 1:2
        for r = 1:length(rateOverlap{g,c})
            if isempty(rateOverlap{g,c}{r})
                continue
            end
            hold on
            ratCntr = ratCntr + 1;
            h = dotplot(ratCntr, {rateOverlap{g,c}{r}}, tmpJit, rgb(cols{g}), rgb('White'));
            line([ratCntr-0.5 ratCntr-0.5], [-0.5 1], 'Color', 'k')
            line([ratCntr+0.5 ratCntr+0.5], [-0.5 1], 'Color', 'k')
            tickLabs = cat(1, tickLabs, group(g).rat(r).name);

            quants = quantile(rateOverlap{g,c}{r}, 3);
            tmpCI = bootci(nboot, {@mean, rateOverlap{g,c}{r}}, 'type', 'per');
            meanData = mean(rateOverlap{g,c}{r});
            %         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
            pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
            line([ratCntr-0.25 ratCntr+0.25], [meanData meanData], 'Color', 'k')
        end %rat
    end %group

    title(combosTxt{c})
        axis square
    xlim([0 13])
    xticks(1:12)
    xticklabels(tickLabs)
    xlabel('Rat #')
    ylim([0 1])
    yticks(0:0.5:1)
    if c == 1
        ylabel('Rate overlap')
    end
end %combo

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end

%% SPAT INFO

figtitle = 'SpatialInfo_byRat_CI';
figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = {};
for g = 1:2
    for r = 1:length(spatInfo{g})
        if isempty(spatInfo{g}{r})
            continue
        end
        hold on
        ratCntr = ratCntr + 1;
        h = dotplot(ratCntr, {spatInfo{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([ratCntr-0.5 ratCntr-0.5], [0 5000], 'Color', 'k')
        line([ratCntr+0.5 ratCntr+0.5], [0 5000], 'Color', 'k')
        tickLabs = cat(1, tickLabs, group(g).rat(r).name);

        quants = quantile(spatInfo{g}{r}, 3);
        tmpCI = bootci(nboot, {@mean, spatInfo{g}{r}}, 'type', 'per');
        meanData = mean(spatInfo{g}{r});
        %         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
        pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
        line([ratCntr-0.25 ratCntr+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 13])
xticks(1:12)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0 6])
yticks(0:1:6)
ylabel('Spatial information (bits/spike)')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end

%% FIGS - SPAT INFO SORTED

figtitle = 'SpatialInfo_byRat_CI_sorted';
figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = cell(12,1);
for g = 1:2
    tmpMeans = cellfun(@mean, spatInfo{g});
    [~, sortOrd] = sort(tmpMeans);
    %     sortOrd = sortOrd - sum(isnan(tmpMeans));
    for r = 1:length(spatInfo{g})
        if isempty(spatInfo{g}{r})
            continue
        end
        hold on
        %         ratCntr = ratCntr + 1;
        xVal = find(r==sortOrd) + 6*(g-1);
        h = dotplot(xVal, {spatInfo{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([xVal-0.5 xVal-0.5], [0 10], 'Color', 'k')
        line([xVal+0.5 xVal+0.5], [0 10], 'Color', 'k')
        tickLabs{xVal} = group(g).rat(r).name;

        quants = quantile(spatInfo{g}{r}, 3);
        tmpCI = bootci(nboot, {@mean, spatInfo{g}{r}}, 'type', 'per');
        meanData = mean(spatInfo{g}{r});
        %         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
        pt = patch([xVal-0.25 xVal-0.25 xVal+0.25 xVal+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
        line([xVal-0.25 xVal+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group


xlim([0 13])
xticks(1:12)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0 6])
yticks(0:1:6)
ylabel('Spatial information (bits/spike)')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end

%% SPARSITY

figtitle = 'SpatialSpar_byRat_CI';
figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = {};
for g = 1:2
    for r = 1:length(sparsity{g})
        if isempty(sparsity{g}{r})
            continue
        end
        hold on
        ratCntr = ratCntr + 1;
        h = dotplot(ratCntr, {sparsity{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([ratCntr-0.5 ratCntr-0.5], [0 10], 'Color', 'k')
        line([ratCntr+0.5 ratCntr+0.5], [0 10], 'Color', 'k')
        tickLabs = cat(1, tickLabs, group(g).rat(r).name);

        quants = quantile(sparsity{g}{r}, 3);
        tmpCI = bootci(nboot, {@mean, sparsity{g}{r}}, 'type', 'per');
        meanData = mean(sparsity{g}{r});
        %         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
        pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
        line([ratCntr-0.25 ratCntr+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 13])
xticks(1:12)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0 1])
yticks(0:0.25:1)
ylabel('Spatial sparsity')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end

%% SPAT SPAR SORTED

figtitle = 'SpatialSpar_byRat_CI_sorted';
figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = cell(12,1);
for g = 1:2
    tmpMeans = cellfun(@mean, sparsity{g});
    [~, sortOrd] = sort(tmpMeans);
    %     sortOrd = sortOrd - sum(isnan(tmpMeans));
    for r = 1:length(sparsity{g})
        if isempty(sparsity{g}{r})
            continue
        end
        hold on
        %         ratCntr = ratCntr + 1;
        xVal = find(r==sortOrd) + 6*(g-1);
        h = dotplot(xVal, {sparsity{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([xVal-0.5 xVal-0.5], [0 1], 'Color', 'k')
        line([xVal+0.5 xVal+0.5], [0 1], 'Color', 'k')
        tickLabs{xVal} = group(g).rat(r).name;

        quants = quantile(sparsity{g}{r}, 3);
        tmpCI = bootci(nboot, {@mean, sparsity{g}{r}}, 'type', 'per');
        meanData = mean(sparsity{g}{r});
        %         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
        pt = patch([xVal-0.25 xVal-0.25 xVal+0.25 xVal+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
        line([xVal-0.25 xVal+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group


xlim([0 13])
xticks(1:12)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0 1])
yticks(0:0.25:1)
ylabel('Spatial sparsity')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end

%% PLACE FIELD SIZE

figtitle = 'PlaceFieldSize_byRat_CI';
figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = {};
for g = 1:2
    for r = 1:length(fieldSize{g})
        if isempty(fieldSize{g}{r})
            continue
        end
        hold on
        ratCntr = ratCntr + 1;
        h = dotplot(ratCntr, {fieldSize{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([ratCntr-0.5 ratCntr-0.5], [0 150], 'Color', 'k')
        line([ratCntr+0.5 ratCntr+0.5], [0 150], 'Color', 'k')
        tickLabs = cat(1, tickLabs, group(g).rat(r).name);

        quants = quantile(fieldSize{g}{r}, 3);

        tmpCI = bootci(nboot, {@mean, fieldSize{g}{r}}, 'type', 'per');
        meanData = mean(fieldSize{g}{r});
        %         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
        pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
        line([ratCntr-0.25 ratCntr+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 13])
xticks(1:12)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0 125])
yticks(0:25:125)
ylabel('Place field size (cm)')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end

%% PLACE FIELD SIZE SORTED

figtitle = 'PlaceFieldSize_byRat_CI_sorted';
figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = cell(12,1);
for g = 1:2
    tmpMeans = cellfun(@mean, fieldSize{g});
    [~, sortOrd] = sort(tmpMeans);
    %     sortOrd = sortOrd - sum(isnan(tmpMeans));
    for r = 1:length(fieldSize{g})
        if isempty(fieldSize{g}{r})
            continue
        end
        hold on
        %         ratCntr = ratCntr + 1;
        xVal = find(r==sortOrd) + 6*(g-1);
        h = dotplot(xVal, {fieldSize{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([xVal-0.5 xVal-0.5], [0 150], 'Color', 'k')
        line([xVal+0.5 xVal+0.5], [0 150], 'Color', 'k')
        tickLabs{xVal} = group(g).rat(r).name;

        quants = quantile(fieldSize{g}{r}, 3);
        tmpCI = bootci(nboot, {@mean, fieldSize{g}{r}}, 'type', 'per');
        meanData = mean(fieldSize{g}{r});
        %         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
        pt = patch([xVal-0.25 xVal-0.25 xVal+0.25 xVal+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
        line([xVal-0.25 xVal+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group


xlim([0 13])
xticks(1:12)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0 125])
yticks(0:25:150)
ylabel('Place field size (cm)')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end

%% MEAN FIRING RATE

figtitle = 'MeanFiringRate_byRat_CI';
figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = {};
for g = 1:2
    for r = 1:length(meanFirRate{g})
        if isempty(meanFirRate{g}{r})
            continue
        end
        hold on
        ratCntr = ratCntr + 1;
        h = dotplot(ratCntr, {meanFirRate{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([ratCntr-0.5 ratCntr-0.5], [0 100], 'Color', 'k')
        line([ratCntr+0.5 ratCntr+0.5], [0 100], 'Color', 'k')
        tickLabs = cat(1, tickLabs, group(g).rat(r).name);

        quants = quantile(meanFirRate{g}{r}, 3);
        tmpCI = bootci(nboot, {@mean, meanFirRate{g}{r}}, 'type', 'per');
        meanData = mean(meanFirRate{g}{r});
        %         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
        pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
        line([ratCntr-0.25 ratCntr+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 13])
xticks(1:12)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0 5])
yticks(0:1:5)
ylabel('Mean firing rate (Hz)')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option

%% MEAN FIRING RATE SORTED

figtitle = 'MeanFiringRate_byRat_CI_sorted';
figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = cell(12,1);
for g = 1:2
    tmpMeans = cellfun(@mean, meanFirRate{g});
    [~, sortOrd] = sort(tmpMeans);
    %     sortOrd = sortOrd - sum(isnan(tmpMeans));
    for r = 1:length(meanFirRate{g})
        if isempty(meanFirRate{g}{r})
            continue
        end
        hold on
        %         ratCntr = ratCntr + 1;
        xVal = find(r==sortOrd) + 6*(g-1);
        h = dotplot(xVal, {meanFirRate{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([xVal-0.5 xVal-0.5], [0 100], 'Color', 'k')
        line([xVal+0.5 xVal+0.5], [0 100], 'Color', 'k')
        tickLabs{xVal} = group(g).rat(r).name;

        quants = quantile(meanFirRate{g}{r}, 3);
        tmpCI = bootci(nboot, {@mean, meanFirRate{g}{r}}, 'type', 'per');
        meanData = mean(meanFirRate{g}{r});
        %         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
        pt = patch([xVal-0.25 xVal-0.25 xVal+0.25 xVal+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
        line([xVal-0.25 xVal+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group


xlim([0 13])
xticks(1:12)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0 5])
yticks(0:1:5)
ylabel('Mean firing rate (Hz)')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end

%% PEAK FIRING RATE

figtitle = 'PeakPlaceFieldFiringRate_byRat_CI';
figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = {};
for g = 1:2
    for r = 1:length(peakFieldFirRate{g})
        if isempty(peakFieldFirRate{g}{r})
            continue
        end
        hold on
        ratCntr = ratCntr + 1;
        h = dotplot(ratCntr, {peakFieldFirRate{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([ratCntr-0.5 ratCntr-0.5], [0 100], 'Color', 'k')
        line([ratCntr+0.5 ratCntr+0.5], [0 100], 'Color', 'k')
        tickLabs = cat(1, tickLabs, group(g).rat(r).name);

        quants = quantile(peakFieldFirRate{g}{r}, 3);
        tmpCI = bootci(nboot, {@mean, peakFieldFirRate{g}{r}}, 'type', 'per');
        meanData = mean(peakFieldFirRate{g}{r});
        %         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
        pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
        line([ratCntr-0.25 ratCntr+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 13])
xticks(1:12)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0 30])
yticks(0:5:30)
ylabel('Peak in-field firing rate (Hz)')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option


%% PEAK FIRING RATE SORTED

figtitle = 'PeakPlaceFieldFiringRate_byRat_CI_sorted';
figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = cell(12,1);
for g = 1:2
    tmpMeans = cellfun(@mean, peakFieldFirRate{g});
    [~, sortOrd] = sort(tmpMeans);
    %     sortOrd = sortOrd - sum(isnan(tmpMeans));
    for r = 1:length(peakFieldFirRate{g})
        if isempty(peakFieldFirRate{g}{r})
            continue
        end
        hold on
        %         ratCntr = ratCntr + 1;
        xVal = find(r==sortOrd) + 6*(g-1);
        h = dotplot(xVal, {peakFieldFirRate{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([xVal-0.5 xVal-0.5], [0 100], 'Color', 'k')
        line([xVal+0.5 xVal+0.5], [0 100], 'Color', 'k')
        tickLabs{xVal} = group(g).rat(r).name;

        quants = quantile(peakFieldFirRate{g}{r}, 3);
        tmpCI = bootci(nboot, {@mean, peakFieldFirRate{g}{r}}, 'type', 'per');
        meanData = mean(peakFieldFirRate{g}{r});
        %         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
        pt = patch([xVal-0.25 xVal-0.25 xVal+0.25 xVal+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
        line([xVal-0.25 xVal+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 13])
xticks(1:12)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0 30])
yticks(0:5:30)
ylabel('Peak in-field firing rate (Hz)')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end


%% AVG IN-FIELD FIRING RATE

figtitle = 'AvgPlaceFieldFiringRate_byRat_CI';
figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = {};
for g = 1:2
    for r = 1:length(meanFieldFirRate{g})
        if isempty(meanFieldFirRate{g}{r})
            continue
        end
        hold on
        ratCntr = ratCntr + 1;
        h = dotplot(ratCntr, {meanFieldFirRate{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([ratCntr-0.5 ratCntr-0.5], [0 100], 'Color', 'k')
        line([ratCntr+0.5 ratCntr+0.5], [0 100], 'Color', 'k')
        tickLabs = cat(1, tickLabs, group(g).rat(r).name);

        quants = quantile(meanFieldFirRate{g}{r}, 3);
        tmpCI = bootci(nboot, {@mean, meanFieldFirRate{g}{r}}, 'type', 'per');
        meanData = mean(meanFieldFirRate{g}{r});
        %         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
        pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
        line([ratCntr-0.25 ratCntr+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 13])
xticks(1:12)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0 20])
yticks(0:5:30)
ylabel('Average in-field firing rate (Hz)')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option

%% AVG IN-FIELD SORTED

figtitle = 'AvgPlaceFieldFiringRate_byRat_CI_sorted';
figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = cell(12,1);
for g = 1:2
    tmpMeans = cellfun(@mean, meanFieldFirRate{g});
    [~, sortOrd] = sort(tmpMeans);
    %     sortOrd = sortOrd - sum(isnan(tmpMeans));
    for r = 1:length(meanFieldFirRate{g})
        if isempty(meanFieldFirRate{g}{r})
            continue
        end
        hold on
        %         ratCntr = ratCntr + 1;
        xVal = find(r==sortOrd) + 6*(g-1);
        h = dotplot(xVal, {meanFieldFirRate{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([xVal-0.5 xVal-0.5], [0 100], 'Color', 'k')
        line([xVal+0.5 xVal+0.5], [0 100], 'Color', 'k')
        tickLabs{xVal} = group(g).rat(r).name;

        quants = quantile(meanFieldFirRate{g}{r}, 3);
        tmpCI = bootci(nboot, {@mean, meanFieldFirRate{g}{r}}, 'type', 'per');
        meanData = mean(meanFieldFirRate{g}{r});
        %         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
        pt = patch([xVal-0.25 xVal-0.25 xVal+0.25 xVal+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
        line([xVal-0.25 xVal+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 13])
xticks(1:12)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0 20])
yticks(0:5:30)
ylabel('Average in-field firing rate (Hz)')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end

end %function