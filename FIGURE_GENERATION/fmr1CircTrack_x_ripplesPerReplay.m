function fmr1CircTrack_x_ripplesPerReplay(group)
% function fmr1CircTrack_x_ripplesPerReplay(group)
%
% PURPOSE:
%   Plot the number of ripples that co-occur with replay events.
%
% INPUT:
%   group = data struct, through function 6.
%
% OUTPUT:
%   Figures.
%
% MMD
% 02/2024
% Colgin Lab

%% OPTIONS

saveOrNot = 1; %to save figs

%% INITIALIZE

saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\REPLAY\ripCoOccur';
curDir = pwd;

ripPerRep = cell(2,1);
ripxEvDur = cell(2,1);
ripxEvDurAlt = cell(2,1);

Fmax = 250; %250Hz
Fmin = 150; %150Hz
Fs = 2000;

[B,A]=butter(3,[Fmin/(Fs*0.5) Fmax/(Fs*0.5)]);

gWinStd = 8/1000; %Davidson, Kloosterman, & Wilson 2009
gWinStd = gWinStd *Fs; %convert based on sampling freq
gKrnl = gausskernel(gWinStd, gWinStd);

pkCut = 3;

r2Thresh = 0.5;
propCloseThresh = 0;
jumpThresh = 0.25*2*pi; %one quarter of the track

cols = {'Blue', 'Red'};
alpha = 1;

ratCntr = 0;

jitter = 0.15;
cols = {'Blue', 'Red'};
groupNames = {'WT', 'FXS'};

statAll = [];

%% GET DATA

for g = 1:2
    fprintf('%s\n', group(g).name)
    for r = 1:length(group(g).rat)
        fprintf('\t%s\n', group(g).rat(r).name)
        ratCntr = ratCntr + 1;
        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day))
            if group(g).rat(r).day(d).decThresh == 0
                continue
            end %not enough to decode

            tetNums = group(g).rat(r).day(d).tetNums;

            for s = 2:5 %start from 2, after experience
                fprintf('\t\t\tSleep %d\n', s)
                if isempty(group(g).rat(r).day(d).sleep(s).unit) %if there is anything in the sleep
                    continue %to next event
                end %if there's any data

                events = group(g).rat(r).day(d).sleep(s).popEv;

                fprintf('\t\t\t\tDetecting ripples\n')
                
                 allIDs = vertcat(group(g).rat(r).day(d).xBegUnitInfo(:).ID);
            tetNum = mode(allIDs(:,1)); %use tetrode w the msot cells
             
                    lfpStruct = read_in_lfp([group(g).rat(r).day(d).sleep(s).dir '\CSC' num2str(tetNum) '.ncs']);
                    tmpFilt = filtfilt(B, A, lfpStruct.data);
                    hilbFilt = abs(hilbert(tmpFilt));

                smData = conv(hilbFilt, gKrnl, 'same');
                zData = zscore(smData);
                zci = find(diff(sign(zData))); %find where zscored firing rate crosses 0
                [~, pkLocs] = findpeaks(zData, 'MinPeakHeight', pkCut); %find where it peaks
                pkTms = lfpStruct.ts(pkLocs);

                fprintf('\t\t\t\tMatching ripples to replay events\n')
                for i = 1:length(events)
                  if isnan(events(i).r2) || events(i).r2 < r2Thresh || events(i).propClose < propCloseThresh || events(i).maxJump > jumpThresh %it doesn't meet threshold
                        continue %to next event
                    end %r2 thresh

                    pkInds = find(pkTms >= events(i).tms(1) & pkTms <= events(i).tms(2));
                    evRips = pkTms(pkInds);
                    ripPerRep{g} = [ripPerRep{g}; length(evRips)];

                    ripDurs = [];
                    for p = 1:length(evRips)

                        startEdgeInd = find(zci < pkLocs(pkInds(p)), 1, 'Last');
                        startEdge = zci(startEdgeInd);

                        endEdgeInd = find(zci > pkLocs(pkInds(p)), 1, 'First');
                        endEdge = zci(endEdgeInd);

                        if endEdge < startEdge
                            keyboard
                        end

                        ripDurs = [ripDurs; (endEdge-startEdge) ./ lfpStruct.Fs];
                    end %rip

                    ripxEvDur{g} = [ripxEvDur{g}; length(evRips) diff(events(i).tms) mean(ripDurs)];

tmpStat = [g r length(evRips) diff(events(i).tms)];
statAll = [statAll; tmpStat]; %group rat #rips repDur ripDur ripDurAlt
                    

                end %i - events
            end %sleep
        end %day
    end %rat
end %group

cd(saveDir)

dummyGroup = statAll(:,1)==2; %dummy coding the group
tmpMeanCent = statAll(:,4) - mean(statAll(:,4)); %mean center the replay event duration
tmpInt = dummyGroup .* tmpMeanCent;
statAll = cat(2, statAll, dummyGroup, tmpMeanCent, tmpInt);
keyboard

%% CORRELATION - NUM RIP

figtitle = 'NumRip_Duration';
% figure('Name', figtitle, 'Position', [680 558 460 420])
figure('Name', figtitle)

for g = 1:2
    hold on
subplot(1,2,g)
    xData = ripxEvDur{g}(:,2);
    yData = ripxEvDur{g}(:,1);

    tmpJitter = -jitter + (jitter*2)*rand(length(yData),1);
    tmpJitter = abs(tmpJitter);
yData = yData + tmpJitter;
    scatter(xData, yData, 20, 'MarkerFaceColor', rgb(cols{g}), 'MarkerFaceAlpha', 0.1, 'MarkerEdgeColor', rgb(cols{g}), 'MarkerEdgeAlpha', 0.5)


hold on
     fit = polyfit(ripxEvDur{g}(:,2), ripxEvDur{g}(:,1), 1);
     [rho, pVal] = corr(ripxEvDur{g}(:,2), ripxEvDur{g}(:,1));
    xFit = [min(ripxEvDur{g}(:,2)) max(ripxEvDur{g}(:,2))];
    yFit = xFit * fit(1) + fit(2);
    plot(xFit, yFit, 'Color', rgb(cols{g}), 'LineWidth', 1, 'LineStyle', '--')

    title(['R = ' num2str(rho) ', p = ' num2str(pVal)])
 axis square
    ylabel('Number of ripples')
ax = gca;
ylim([-0.5 ax.YLim(2)])
xlabel('Replay event duration (s)')

end %g

same_axes;


if saveOrNot == 1
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    set(gcf,'renderer','Painters')
    saveas(gcf, figtitle, 'epsc');
end



end %function