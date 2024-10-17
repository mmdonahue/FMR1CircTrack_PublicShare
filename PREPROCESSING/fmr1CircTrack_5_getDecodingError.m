function group = fmr1CircTrack_5_getDecodingError(group)
% function group = fmr1CircTrack_5_getDecodingError(group)
%
% PURPOSE:
%   Get decoding error to determine if recorded place cell population is sufficient to perform
%   analyses that require Bayesian decoding.
%
% MMD
% Colgin Lab

% From Hwaun & Colgin, 2019:
% For each rat, 50% of the total probability density was required to be
% less than 20° away from the actual position to be included in further
% analyses. The second measure took the differences between decoded positions
% and actual positions as errors, and cumulative distributions of errors were
% determined (Figure 6b and 6c, bottom panels). Only rats with an error
% distribution reaching 50% at error values less than 20° were included for
% further analysis.

%% INITIALIZE

saveOrNot = 0; %decoding figs

errorAxis = 0:.05:pi;

radThresh = 0.35;
radThreshInd = match(radThresh, errorAxis);
percThresh = 0.5;

Fs = 1000; %1 ms bin size
bayesWin = 0.5; %in sec - same as CZ and EH
bayesStep = bayesWin/5;

runThresh = 5; %cm/s

degBinCtrs = group(2).rat(1).day(1).binCtrs; %doesn't change across days/rats
radBinCtrs = deg2rad(degBinCtrs);

groupNames = {'WT', 'FXS'};
cols = {'Blue', 'Red'};
colMap = define_cust_color_map('white', 'black', 30);
saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\DECODING_GENERAL\decodingByRat'; %for figs
dataDir = 'E:\FMR1_CIRCTRACK\RAW_DATA\';
cd(saveDir);

%% DO DECODING

for g = 1:2
    fprintf('Group %d\n', g);

    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);

        figtitle = [group(g).rat(r).name '_decoding'];
        figure('Name', figtitle, 'Position', [509 121 1092 860])
        dayCntr = 0; %days within rat that meet threshold
        ratConfMats = [];

        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day));
            dayErr = [];

            rateMaps = zeros(length(group(g).rat(r).day(d).xBegUnitInfo), length(radBinCtrs));
            badU = []; %initialize
            for u = 1:size(rateMaps,1)
                if group(g).rat(r).day(d).xBegUnitInfo(u).meetMin
                    rateMaps(u,:) = group(g).rat(r).day(d).xBegUnitInfo(u).smRateMap;
                else
                    badU = [badU; u];
                end %meet min
            end %unit
            rateMaps(badU,:) = [];
            rateMaps(rateMaps==0) = 0.0001; %get rid of zeros because our Bayesian decoder can't handle 'em.

            for b = 1:length(group(g).rat(r).day(d).begin)
                spkTmsByCell = cell(length(group(g).rat(r).day(d).xBegUnitInfo),1);
                for u = 1:length(group(g).rat(r).day(d).xBegUnitInfo)
                    spkTmsByCell{u} = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                end %unit
                spkTmsByCell(badU) = [];

                spkRstr = make_spike_raster(spkTmsByCell, [group(g).rat(r).day(d).begin(b).radPos(1,1) group(g).rat(r).day(d).begin(b).radPos(end,1)], Fs);
                pxn = BayesianDecoder(spkRstr, rateMaps, bayesWin, bayesStep, Fs); %Ernie's decoder

                [nWin, winStartInds] = find_num_windows(size(spkRstr,2), bayesWin*Fs, bayesStep*Fs);
                if winStartInds(end)+bayesWin*Fs < size(spkRstr,2)
                    nWin = nWin + 1;
                    winStartInds(end+1) = winStartInds(end)+bayesStep*Fs; %#ok
                end
                winStartTms = group(g).rat(r).day(d).begin(b).radPos(1,1) + winStartInds/Fs;
                winEndTms = winStartTms + bayesWin;
                pxn(:,nWin+1:size(pxn,2)) = [];
                [~, decPosBins] = max(pxn);
                decPos = radBinCtrs(decPosBins);

                % Get the rat's actual position
                actPosn = nan(1,length(winStartTms));
                instRs = get_runspeed(group(g).rat(r).day(d).begin(b).coords);
                smRs = smooth_runspeed(instRs);
                radPos = group(g).rat(r).day(d).begin(b).radPos;
                for i = 1:length(winStartTms)
                    winSpd = mean(smRs(smRs(:,1)>=winStartTms(i) & smRs(:,1)<winEndTms(i),2));
                    if winSpd > runThresh
                        actPosn(i) = wrapTo360(rad2deg(circ_mean(deg2rad(radPos(smRs(:,1)>=winStartTms(i) & smRs(:,1)<winEndTms(i),2))))); %get mean position across this window
                    end %speed above threshold
                end %each window
                actPosn = deg2rad(actPosn);

                decErr = abs(circ_dist(actPosn, decPos)); %for this begin
                dayErr = [dayErr decErr(~isnan(decErr))]; %store for day

            end %begin

            c = histc(dayErr, errorAxis);
            cumErr = cumsum(c)/sum(c);

            save([dataDir group(g).name '\' group(g).rat(r).name '\' group(g).rat(r).day(d).name '\dayDecodeErr.mat'], 'cumErr', 'dayErr')

            if cumErr(radThreshInd) >= percThresh
                group(g).rat(r).day(d).decThresh = 1; %meets threshold
                dayCntr = dayCntr + 1;

                subplot(4,5,dayCntr+2)
                plot(rad2deg(errorAxis), cumErr, 'Color', rgb(cols{g}))
                title(group(g).rat(r).day(d).name)
                xlabel('Error (degrees)')
                xticks([0 180])
                xlim([0 180])
                ylabel('Cumulative Proportion')
                ylim([0 1])
                yticks([0 1])
                axis square

                subplot(4,5,[1 2 6 7])
                hold on
                plot(rad2deg(errorAxis), cumErr, 'Color', rgb(cols{g}), 'LineStyle', '--')
                title(group(g).rat(r).day(d).name)
                xticks([0 180])
                xlim([0 180])
                ylim([0 1])
                yticks([0 1])
                axis square

                confMat = zeros(length(radBinCtrs),length(radBinCtrs)); %Confusion Matrix sums
                posCount = zeros(length(radBinCtrs),1); %For the denominator to calc confMatrix probability

                % Extract the most probable position from pxn
                decodedPosBins = nan(1,size(pxn,2));
                decodedPosns = nan(1,size(pxn,2));
                for i = 1:size(pxn,2)
                    if sum(sum(diff(pxn(:,i)))) ~= 0 %if probability wasn't equal for all bins (no spikes)
                        tmpPpm = pxn(:,i);
                        [maxVal,decodedPosBins(i)] = max(tmpPpm); %max inds
                        if length(find(tmpPpm == maxVal)) > 1
                            keyboard
                        end %check
                    end %if probability wasn't equal for all bins (no spikes)
                end %time bins
                decodedPosns(~isnan(decodedPosBins)) = radBinCtrs(decodedPosBins(~isnan(decodedPosBins)));

                %Convert actual positions to the closest index
                posInds = match(actPosn,radBinCtrs);
                posInds(isnan(actPosn)) = NaN;

                for p = 1:length(radBinCtrs)
                    % Add up the probability the rat is at each location for the time windows in which the rat’s true position is being indexed
                    confMat(:,p) = confMat(:,p) + nansum(pxn(:,posInds==p),2); % Can't be zero due to precision
                    %The number of time windows for which the rat's true position was equal to the indexed position
                    posCount(p) = posCount(p) + nansum(posInds==p);
                end

                % Calculate the confusion matrix for this day
                confMat = confMat./repmat(posCount', length(radBinCtrs),1);
                confMat(isinf(confMat)) = 0;

                subplot(4,5,dayCntr+12)
                imagesc(degBinCtrs, degBinCtrs, confMat)
                colormap(colMap);
                axis xy
                axis square

                title(group(g).rat(r).day(d).name)
                xticks([0 360])
                xlim([0 360])
                yticks([0 360])
                ylim([0 360])

                xlabel('Actual position (deg)');
                ylabel('Decoded position (deg)');
                caxis([0 0.3])

                % Concatenate across days
                ratConfMats = cat(3, ratConfMats, confMat);
            else
                group(g).rat(r).day(d).decThresh = 0; %does not meet theshold
            end
        end %day

        subplot(4,5,[1 2 6 7])
        title(group(g).rat(r).name)
        xlabel('Error (degrees)')
        ylabel('Cumulative Proportion')

        subplot(4,5,[11 12 16 17])

        imagesc(degBinCtrs, degBinCtrs, mean(ratConfMats,3))
        colormap(colMap);
        axis xy
        axis square

        title(group(g).rat(r).name)

        xticks([0 360])
        xlim([0 360])
        yticks([0 360])
        ylim([0 360])

        xlabel('Actual position (deg)');
        ylabel('Decoded position (deg)');
        caxis([0 0.3])

        cbr = colorbar('SouthOutside');
        cbr.Position = [0.1 0.05 0.35 0.01];
        ylabel(cbr, 'Probability')
        cbr.Ticks = [0:0.1:0.4];

        if saveOrNot == 1
            saveas(gcf, figtitle, 'epsc');
            saveas(gcf, figtitle, 'fig');
            saveas(gcf, figtitle, 'png');
        end %save

    end %rat
end %group


end %function