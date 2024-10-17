function fmr1CircTrack_x_phase_precession(group)
% function fmr1CircTrack_x_phase_precession(group)
%
% PURPOSE:
%   Plot phase precession of neurons on the circle track from WT and KO
%   groups.
%
% INPUT:
%   group struct
%
% OUTPUT:
%   Figures:
%       F(n-n) (optional): Plot showing the phase precession from the whole
%           day for each unit (only spikes included in analysis).
%       F1: Heat plots of spike counts for theta phase x normalized
%           position in place field (as in Bieri et al. 2014). Note: Spike
%           counts will be different for each group since number of cells
%           and cell firing rates are not necessarily the same.
%       F2: Slope of circular-linear regression line.
%       F3: Phase offset of circular-linear regression line.
%       F4: r^2 value of circular-linear regression line.
%
% OPTIONS:
%   See code for options.
%
% MMD
% 7/2021
% Colgin Lab

%% OPTIONS

saveOrNot = 0; %1 to save figs, 0 to not

makeUnitPlots = 0; %1 to make the plots, 0 to not

downSampSpkCnts = 0; %to down sample spike counts on the first plot to be equal between groups

degBinSz = 20; %degrees per bin for theta phase
posBinSz = 0.01; %normalized distance into place field


%% INITIALIZE

velFilt = 1;
durCrit = 1;
runThresh = 5; %cm/s

minSpks = 50; %minimum number of spikes across the pf for each begin for cell to be included
spatCorMin = 0; %minimum spat coherence between first and last run sessions

numDegBins = 720/degBinSz; %two theta cycles
numPosBins = 1/posBinSz;

spatBinSz = 4;
degCmConv = (pi*100)/360; %track has 1 m diameter

degxPos = zeros(numDegBins,numPosBins,2); %initialize, by group
spkInfo = cell(2,1); %for down sampling spikes

ppSlopes = cell(2,1); %phase precession slopes by group
ppSlopesCm = cell(2,1);
ppPhaseOffset = cell(2,1);
ppPhaseRange = cell(2,1); %phase range = phase difference between the
% spike with the highest shifted phase and the spike with the lowest
% shifted phase (shifted phase is minus phase offset)
ppPhaseRangeAlt = cell(2,1); %alt phase range - pp slope * pf size
ppR2 = cell(2,1); %r2 values from circular regression

ppSlopesRat = cell(2,1); %phase precession slopes by group
ppSlopesRat(:) = {cell(6,1)};
ppSlopesCmRat = cell(2,1);
ppSlopesCmRat(:) = {cell(6,1)};
ppPhaseOffsetRat = cell(2,1);
ppPhaseOffsetRat(:) = {cell(6,1)};
ppPhaseRangeRat = cell(2,1); %phase range = phase difference between the
ppPhaseRangeRat(:) = {cell(6,1)};
% spike with the highest shifted phase and the spike with the lowest
% shifted phase (shifted phase is minus phase offset)
ppPhaseRangeAltRat = cell(2,1); %alt phase range - pp slope * pf size
ppPhaseRangeAltRat(:) = {cell(6,1)};
ppR2Rat = cell(2,1); %r2 values from circular regression
ppR2Rat(:) = {cell(6,1)};
ppRhoRat = cell(2,1); %rho values from circular regression
ppRhoRat(:) = {cell(6,1)};

cols = {'Blue', 'Red'};
alpha = 0.5;
groupNames = {'WT', 'FXS'};

saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\THETA_GENERAL\phasePrecession';
curDir = pwd;

statAll = [];
ratCntr = 0;

%% GET DATA

for g = 1:2
    fprintf('%s\n', group(g).name)
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d\n', r, length(group(g).rat))
        ratCntr = ratCntr + 1;
         uCntr = 0;
        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day))
            for u = 1:length(group(g).rat(r).day(d).xBegUnitInfo)
                fprintf('\t\t\tUnit %d/%d\n', u, length(group(g).rat(r).day(d).xBegUnitInfo))
                 uCntr = uCntr + 1;
                tetNum = group(g).rat(r).day(d).xBegUnitInfo(u).ID(1);
                clustNum = group(g).rat(r).day(d).xBegUnitInfo(u).ID(2);
                pf = group(g).rat(r).day(d).xBegUnitInfo(u).pf;
                if isempty(pf)
                    continue %to next cell
                end %no pf

                for p = 1:length(pf) %place field ind
                    spkCheck = 1; %initialize as a pass

                    if spkCheck == 0
                        continue %to next pf
                    end %spk check

                    if pf(p).radPos(end) > pf(p).radPos(1)
                        crossZero = 0;
                    else
                        crossZero = 1;
                    end %determine if pf crosses 0

                    pfLen = rad2deg(circ_dist(deg2rad(pf(p).radPos(end)), deg2rad(pf(p).radPos(1))));
                    pfLenCm = pfLen * degCmConv;

                    uPhisAll = [];
                    uNormPosAll = [];
                    tetNum = group(g).rat(r).day(d).xBegUnitInfo(u).ID(1);

                    for b = 1:length(group(g).rat(r).day(d).begin)
                        begPasses = []; %initialize for this begin
                        begSpks = [];

                        radPos = group(g).rat(r).day(d).begin(b).radPos;
                        coords = group(g).rat(r).day(d).begin(b).coords;
                        runSpeed = smooth_runspeed(get_runspeed(coords));
                        spkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;

                        cd(group(g).rat(r).day(d).begin(b).dir)
                        lfpStruct = read_in_lfp(['CSC' num2str(tetNum) '.ncs']);
                        %                         load(['CSC' num2str(tetNum) '_narrowThetaLfp.mat'], 'filtLfp')
                        %                         lfpStruct.narrowThetaLfp = filtLfp;

                        pfPassBnry = zeros(1,size(radPos,1));
                        if crossZero == 0
                            pfPassBnry(radPos(:,2)>= pf(p).radPos(1) & radPos(:,2)<= pf(p).radPos(end)) = 1;
                        else
                            pfPassBnry(radPos(:,2)>= pf(p).radPos(1) & radPos(:,2)<=360) = 1; %from start-360 and 0-end
                            pfPassBnry(radPos(:,2)>=0 & radPos(:,2)<= pf(p).radPos(end)) = 1;
                        end %cross zero

                        pfPassChunks = bwconncomp(pfPassBnry);
                        for c = 1:length(pfPassChunks.PixelIdxList)
                            tmpInds = pfPassChunks.PixelIdxList{c};

                            passDist = abs(rad2deg(circ_dist(deg2rad(radPos(tmpInds(1),2)), deg2rad(radPos(tmpInds(end),2)))));

                            % Make sure rat traverses almost the whole field
                            if passDist >= pfLen - 5 % a little room for bin rounding error
                                begPasses = [begPasses; radPos(tmpInds(1),1) radPos(tmpInds(end),1)];
                            end %through whole field
                        end %chunks

                        for ps = 1:size(begPasses,1)

                            psSpkTms = spkTms(spkTms >= begPasses(ps,1) & spkTms <= begPasses(ps,2));
                            if isempty(psSpkTms)
                                continue
                            end
                            psRadPos = radPos(radPos(:,1)>= begPasses(ps,1) & radPos(:,1)<= begPasses(ps,2),:);
                            psCoords = coords(coords(:,1)>= begPasses(ps,1) & coords(:,1) <= begPasses(ps,2),:);
                            psSpeed = runSpeed(coords(:,1)>= begPasses(ps,1) & coords(:,1) <= begPasses(ps,2),:);

                            if min(psSpeed(:,2)) < runThresh
                                continue
                            end %check run speed through the pass

                            [~,~,~,spkCnts] = get_ratemap_circtrack(psSpkTms, psCoords, psRadPos, spatBinSz, velFilt, durCrit);
                            if sum(spkCnts(pf(p).inds)) < 3 %at least 3 spikes
                                continue
                            end %at least 3 spikes
                            if length(find(spkCnts(pf(p).inds))) < 2 %in at least two different position bins
                                continue
                            end %at least 2 pos bins

                            begSpks = [begSpks psSpkTms'];
                        end %pass

                        if isempty(begSpks)
                            continue
                        end %no spikes

                        pp = phase_precession_circtrack_pf(begSpks, radPos, coords, pf(p), lfpStruct);

                        if ~isempty(pp.spkTms)
                            uPhisAll = [uPhisAll pp.spkPhis];
                            uNormPosAll = [uNormPosAll pp.normSpkPos];
                        end %if
                    end %begin

                    if length(uNormPosAll) < minSpks
                        continue %to next pf
                    end %not enough spikes

                    tmpDegxPos = zeros(numDegBins/2,numPosBins); %initialize %first 360 degrees
                    degBins = 0:degBinSz:360;

                    allPhis = pp.spkPhis; %get all of the spike theta phases

                    for db = 1:numDegBins/2
                        degRange = [degBins(db) degBins(db+1)];
                        pullPos = uNormPosAll(uPhisAll > degRange(1) & uPhisAll <= degRange(2));

                        tmpPosBins = histcounts(pullPos, 0:posBinSz:1);
                        tmpDegxPos(db,:) = tmpPosBins;
                    end %degrees bins

                    tmpTwoCycle = [tmpDegxPos; tmpDegxPos];
                    degxPos(:,:,g) = degxPos(:,:,g) + tmpTwoCycle;

                    beta = CircularRegression(uNormPosAll, deg2rad(uPhisAll));
                    slope = rad2deg(beta(1));
                    phaseOff = rad2deg(wrapTo2Pi(beta(2)));
                    cmSlope = slope / pfLenCm;

                    altPhaseRange = cmSlope * pfLenCm;

                    phaseRange = max(uPhisAll) - min(uPhisAll);
                    
                    rho = circ_corrcl(deg2rad(uPhisAll), uNormPosAll);
                    tmpCal = beta(1)*uNormPosAll + beta(2);
                    R2 = circLinRegress_r2(deg2rad(uPhisAll), tmpCal);
%                     if isnan(R2)
%                         R2 = [];
%                     end

                    if makeUnitPlots == 1
                        cd(saveDir)
                        cd('unitPlots')
                        cd(group(g).name)
                        try cd(group(g).rat(r).name)
                        catch
                            mkdir(group(g).rat(r).name)
                            cd(group(g).rat(r).name)
                        end
                        figtitle = [group(g).rat(r).day(d).name '_TT' num2str(tetNum) '_' num2str(clustNum)];
                        figure('Name', figtitle)

                        hold on
                        plot(uNormPosAll, uPhisAll, 'k.')
                        plot(uNormPosAll, (uPhisAll+360), 'k.')
                        ylim([0 720])

                        calphase = slope*[0 1] + phaseOff;
                        plot([0 1], calphase, 'r')
                        plot([0 1], calphase+360, 'r')
                        if max(calphase+360) < 720
                            plot([0 1], calphase+720, 'r')
                        end %in range

                        ylim([0 720])
                        ylabel('Theta phase (deg)')

                        xlim([0 1])
                        xticks([0:0.5:1])
                        xlabel({'Normalized position in place field', ['slope = ' num2str(round(slope)) ' deg/pf | ' num2str(round(cmSlope,1)) ' deg/cm']})

                        title(['TT' num2str(tetNum) '\_' num2str(clustNum) ])

                        if saveOrNot == 1
                            saveas(gcf, figtitle, 'epsc');
                            saveas(gcf, figtitle, 'fig');
                            saveas(gcf, figtitle, 'png');
                        end %save
                    end %make unit plots

                    ppSlopes{g} = [ppSlopes{g} slope];
                    ppSlopesRat{g}{r} = [ppSlopesRat{g}{r} slope];
                    ppSlopesCm{g} = [ppSlopesCm{g} cmSlope];
                    ppSlopesCmRat{g}{r} = [ppSlopesCmRat{g}{r} cmSlope];
                    ppPhaseOffset{g} = [ppPhaseOffset{g} phaseOff];
                    ppPhaseOffsetRat{g}{r} = [ppPhaseOffsetRat{g}{r} phaseOff];
                    ppPhaseRange{g} = [ppPhaseRange{g} phaseRange];
                    ppPhaseRangeRat{g}{r} = [ppPhaseRangeRat{g}{r} phaseRange];
                    ppPhaseRangeAlt{g} = [ppPhaseRangeAlt{g} altPhaseRange];
                    ppPhaseRangeAltRat{g}{r} = [ppPhaseRangeAltRat{g}{r} altPhaseRange];
                    ppR2{g} = [ppR2{g} R2];
                    ppR2Rat{g}{r} = [ppR2Rat{g}{r} R2];
%                        ppR2{g} = [ppR2{g} R2];
                    ppRhoRat{g}{r} = [ppRhoRat{g}{r} rho];

                    tmpStat = [g ratCntr d uCntr slope cmSlope phaseOff phaseRange altPhaseRange R2 rho];
                    try
                    statAll = [statAll; tmpStat];
                    catch; keyboard; end
                end %place field ind
            end %unit
            close all
        end %day
    end %rat
end %group

cd(saveDir)
keyboard

%% FIG 1 - PHASE X POS HEATMAP

figtitle = 'PhasePrecession_ThetaPhasexPostion';
figure('Name', figtitle, 'Position', [529 461 940 420])

% maxSpkCnt = max(degxPos(:));
% normDegxPos = degxPos ./ maxSpkCnt;

for g = 1:2
    normSpkCnts = degxPos(:,:,g) ./ max(max(degxPos(:,:,g)));

    subplot(1,2,g)
    imagesc(linspace(0, 1, size(degxPos,2)), linspace(0, 720, size(degxPos,1)), normSpkCnts)
    colormap(jet)
    axis xy

    xlabel('Position in place field')

    if g == 1
        ylabel('Theta phase (deg)')
    end

    cbr = colorbar;
    ylabel(cbr, 'Normalized spike Counts')
    %     ax = gca;
    caxis([0 1])
    title([groupNames{g} ' (n = ' num2str(length(ppSlopes{g})) ' cells)'])
end %group

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIG 2 - PHASE X POS HEATMAP RECONSTRUCTED

figtitle = 'PhasePrecession_ThetaPhasexPostion_reconstructed';
figure('Name', figtitle, 'Position', [558 332 748 420])

% maxSpkCnt = max(degxPos(:));
% normDegxPos = degxPos ./ maxSpkCnt;

radBinCtrs = deg2rad(degBinSz/2:degBinSz:360);
for g = 1:2
    subplot(1,2,g)

    com = zeros(1,size(degxPos,2)/2);
    for p = 1:length(com)
        tmpCom = circ_mean(radBinCtrs', degxPos(1:size(degxPos,1)/2,p,g));
        com(p) = wrapTo360(rad2deg(tmpCom));
    end %pos bins

    plot(linspace(0, 1, length(com)), com, 'k.')
    hold on;
    plot(linspace(0, 1, length(com)), com+360, 'k.')

    beta = CircularRegression(linspace(0, 1, length(com)), deg2rad(com));
    calphase = wrapTo360(rad2deg(beta(1)*[0 1] + beta(2)));

    plot([0 1], calphase, 'b')
    plot([0 1], calphase+360, 'b')

    xlabel('Position in place field')

    if g == 1
        ylabel('Theta phase (deg)')
    end

    ylim([0 720])

    title({[groupNames{g} ' (n = ' num2str(length(ppSlopes{g})) ' cells)'],...
        ['slope = ' num2str(round(rad2deg(beta(1)))) ' deg/pf, phase offset = '...
        num2str(round(rad2deg(beta(2)))) ' deg']})
end %group

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end


%% FIG 3 - CM SLOPE

figtitle = 'PhasePrecession_Slope_byRat_CI';
figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = {};
for g = 1:2
    for r = 1:length(ppSlopesCmRat{g})
        if isempty(ppSlopesCmRat{g}{r})
            continue
        end
        hold on
        ratCntr = ratCntr + 1;
        h = dotplot(ratCntr, {ppSlopesCmRat{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([ratCntr-0.5 ratCntr-0.5], [-20 20], 'Color', 'k')
        line([ratCntr+0.5 ratCntr+0.5], [-20 20], 'Color', 'k')
        tickLabs = cat(1, tickLabs, group(g).rat(r).name);

          quants = quantile(ppSlopesCmRat{g}{r}, 3);
          tmpCI = bootci(nboot, {@mean, ppSlopesCmRat{g}{r}}, 'type', 'per');
    meanData = mean(ppSlopesCmRat{g}{r});
%         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
        pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
    line([ratCntr-0.25 ratCntr+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 13])
xticks(1:12)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([-20 20])
ylabel('Slope (deg/cm)')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
        saveas(gcf, figtitle, 'epsc');
end

%% FIGS

figtitle = 'PhasePrecession_Slope_byRat_CI_sorted';
figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = cell(8,1);
for g = 1:2
    tmpMeans = cellfun(@mean, ppSlopesCmRat{g});
    [~, sortOrd] = sort(tmpMeans);
    for r = 1:length(ppSlopesCmRat{g})
        if isempty(ppSlopesCmRat{g}{r})
            continue
        end
        hold on
%         ratCntr = ratCntr + 1;
     xVal = find(r==sortOrd) + 6*(g-1);
        h = dotplot(xVal, {ppSlopesCmRat{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([xVal-0.5 xVal-0.5], [-20 20], 'Color', 'k')
        line([xVal+0.5 xVal+0.5], [-20 20], 'Color', 'k')
        tickLabs{xVal} = group(g).rat(r).name;

          quants = quantile(ppSlopesCmRat{g}{r}, 3);
            tmpCI = bootci(nboot, {@mean, ppSlopesCmRat{g}{r}}, 'type', 'per');
    meanData = mean(ppSlopesCmRat{g}{r});
%         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
pt = patch([xVal-0.25 xVal-0.25 xVal+0.25 xVal+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
    line([xVal-0.25 xVal+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 13])
xticks(1:12)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([-20 20])
ylabel('Slope (deg/cm)')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
        saveas(gcf, figtitle, 'epsc');
end

%% FIGS - SLOPE BY RAT COMBINED

figtitle = 'PhasePrecession_Slope_byRat_combined';
figure('Name', figtitle)

ratCntr = 0;
tmpJit = 0.3;
markers = {'o', 'x', '^', 'square', '*', '+', 'x', 'diamond', 'o', '^', '*', '+'};
fCols = cat(1, [1 1 1], rgb('Blue'), rgb('Blue'), rgb('Blue'), [1 1 1], rgb('Blue'),  [1 1 1], rgb('Red'), rgb('Red'), [1 1 1], rgb('Red'), rgb('Red'), rgb('Red'));
lh = zeros(1,7);
leg = cell(1,7);
for g = 1:2
    for r = 1:length(ppSlopesCmRat{g})
        if isempty(ppSlopesCmRat{g}{r})
            continue
        end
        hold on
        ratCntr = ratCntr + 1;

        xData = rescale(randn(length(ppSlopesCmRat{g}{r}),1), g-tmpJit, g+tmpJit);
        lh(ratCntr) = plot(xData, ppSlopesCmRat{g}{r}, 'Color', rgb(cols{g}), 'Marker', markers{ratCntr}, 'LineStyle', 'None', 'MarkerFaceColor', fCols(ratCntr,:));
        leg{ratCntr} = group(g).rat(r).name;
    end %rat

%     line([g-0.25 g+0.25], [tmpMean(g) tmpMean(g)], 'Color', 'Black', 'LineWidth', 3) %mean lines should be longer than the jitter so you can see them

end %group

xticks(1:2)
xticklabels(groupNames)
ylabel('Slope (deg/cm)')
legend(lh, leg, 'Location', 'northeastoutside', 'AutoUpdate', 'off')

zero_line;

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIG - R2 BY RAT

figtitle = 'PhasePrecession_R2_byRat_CI';
figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = {};
for g = 1:2
    for r = 1:length(ppR2Rat{g})
        if isempty(ppR2Rat{g}{r})
            continue
        end
        hold on
        ratCntr = ratCntr + 1;
        h = dotplot(ratCntr, {ppR2Rat{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([ratCntr-0.5 ratCntr-0.5], [0 1], 'Color', 'k')
        line([ratCntr+0.5 ratCntr+0.5], [0 1], 'Color', 'k')
        tickLabs = cat(1, tickLabs, group(g).rat(r).name);

          quants = quantile(ppR2Rat{g}{r}, 3);
          tmpCI = bootci(nboot, {@mean, ppR2Rat{g}{r}}, 'type', 'per');
    meanData = mean(ppR2Rat{g}{r});
%         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
        pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
    line([ratCntr-0.25 ratCntr+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 13])
xticks(1:12)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0 0.5])
yticks(0:0.1:0.5)
ylabel('R^2')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
        saveas(gcf, figtitle, 'epsc');
end

%% FIGS

figtitle = 'PhasePrecession_R2_byRat_CI_sorted';
figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = cell(8,1);
for g = 1:2
    tmpMeans = cellfun(@nanmean, ppR2Rat{g});
    [~, sortOrd] = sort(tmpMeans);
    for r = 1:length(ppR2Rat{g})
        if isempty(ppR2Rat{g}{r})
            continue
        end
        hold on
%         ratCntr = ratCntr + 1;
     xVal = find(r==sortOrd) + 6*(g-1);
        h = dotplot(xVal, {ppR2Rat{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([xVal-0.5 xVal-0.5], [0 1], 'Color', 'k')
        line([xVal+0.5 xVal+0.5], [0 1], 'Color', 'k')
        tickLabs{xVal} = group(g).rat(r).name;

          quants = quantile(ppR2Rat{g}{r}, 3);
            tmpCI = bootci(nboot, {@mean, ppR2Rat{g}{r}(~isnan(ppR2Rat{g}{r}))}, 'type', 'per');
    meanData = nanmean(ppR2Rat{g}{r});
%         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
pt = patch([xVal-0.25 xVal-0.25 xVal+0.25 xVal+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
    line([xVal-0.25 xVal+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 13])
xticks(1:12)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0 0.5])
yticks(0:0.1:0.5)
ylabel('R^2')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
        saveas(gcf, figtitle, 'epsc');
end

%% FIGS - R2 BY RAT COMBINED

figtitle = 'PhasePrecession_R2_byRat_combined';
figure('Name', figtitle)

ratCntr = 0;
tmpJit = 0.3;
markers = {'o', 'x', '^', 'square', '*', '+', 'x', 'diamond', 'o', '^', '*', '+'};
fCols = cat(1, [1 1 1], rgb('Blue'), rgb('Blue'), rgb('Blue'), [1 1 1], rgb('Blue'),  [1 1 1], rgb('Red'), rgb('Red'), [1 1 1], rgb('Red'), rgb('Red'), rgb('Red'));
lh = zeros(1,7);
leg = cell(1,7);
for g = 1:2
    for r = 1:length(ppR2Rat{g})
        if isempty(ppR2Rat{g}{r})
            continue
        end
        hold on
        ratCntr = ratCntr + 1;

        xData = rescale(randn(length(ppR2Rat{g}{r}),1), g-tmpJit, g+tmpJit);
        lh(ratCntr) = plot(xData, ppR2Rat{g}{r}, 'Color', rgb(cols{g}), 'Marker', markers{ratCntr}, 'LineStyle', 'None', 'MarkerFaceColor', fCols(ratCntr,:));
        leg{ratCntr} = group(g).rat(r).name;
    end %rat

%     line([g-0.25 g+0.25], [tmpMean(g) tmpMean(g)], 'Color', 'Black', 'LineWidth', 3) %mean lines should be longer than the jitter so you can see them

end %group

xticks(1:2)
xticklabels(groupNames)
ylabel('R^2')
legend(lh, leg, 'Location', 'northeastoutside', 'AutoUpdate', 'off')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIGS

figtitle = 'PhasePrecession_rho_byRat_CI_sorted';
figure('Name', figtitle, 'Position',  [680 558 448 420])

ratCntr = 0;
tmpJit = 0.25;
nboot = 1000;
tickLabs = cell(8,1);
for g = 1:2
    tmpMeans = cellfun(@nanmean, ppRhoRat{g});
    [~, sortOrd] = sort(tmpMeans);
    for r = 1:length(ppRhoRat{g})
        if isempty(ppRhoRat{g}{r})
            continue
        end
        hold on
%         ratCntr = ratCntr + 1;
     xVal = find(r==sortOrd) + 6*(g-1);
        h = dotplot(xVal, {ppRhoRat{g}{r}}, tmpJit, rgb(cols{g}), rgb('White'));

        line([xVal-0.5 xVal-0.5], [0 1], 'Color', 'k')
        line([xVal+0.5 xVal+0.5], [0 1], 'Color', 'k')
        tickLabs{xVal} = group(g).rat(r).name;

          quants = quantile(ppRhoRat{g}{r}, 3);
            tmpCI = bootci(nboot, {@mean, ppRhoRat{g}{r}(~isnan(ppRhoRat{g}{r}))}, 'type', 'per');
    meanData = nanmean(ppRhoRat{g}{r});
%         pt = patch([ratCntr-0.25 ratCntr-0.25 ratCntr+0.25 ratCntr+0.25], [quants(1) quants(3) quants(3) quants(1)], 'w', 'FaceAlpha', 0.3);
pt = patch([xVal-0.25 xVal-0.25 xVal+0.25 xVal+0.25], [tmpCI(1) tmpCI(2) tmpCI(2) tmpCI(1)], 'w', 'FaceAlpha', 0.3);
    line([xVal-0.25 xVal+0.25], [meanData meanData], 'Color', 'k')
    end %rat
end %group

xlim([0 13])
xticks(1:12)
xticklabels(tickLabs)
xlabel('Rat #')
ylim([0 0.75])
yticks(0:0.1:0.75)
ylabel('Correlation')

if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
        saveas(gcf, figtitle, 'epsc');
end

end %function