function fmr1CircTrack_x_plotReplayEvents(group)
% function fmr1CircTrack_x_plotReplayEvents(group)
%
% PURPOSE: 
%   Plots the replay events with the highest r^2 values for each genotype.
% 
% MMD
% 05/2024 - but based on code originally written by MMD 06/2020
% Colgin Lab


%% INITIALIZE

pxnForPlot = cell(5,2); %for storing the ppms for plotting
r2ForPlot = zeros(5,2); %for storing the matching r2 values
slopeForPlot = zeros(5,2); %for storing the matching r2 values
ratNameForPlot = cell(5,2);

r2Thresh = 0.5;
propCloseThresh = 0;
jumpThresh = 0.25*2*pi;

radCmConv = (pi*100) / (2*pi); %convert all track measurements from deg to cm before plotting

% Decoding parameters
bayesWin = 20/1000; %20 ms time window, as in Hwaun & Colgin 2019
bayesStep = 10/1000; %10 ms time step, as in Hwaun & Colgin 2019

cols = {'Blue', 'Red'};
groupNames = {'WT', 'FXS'};

saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\REPLAY\replayProperties';

%% GET DATA

for g = 1:2
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
            if group(g).rat(r).day(d).decThresh == 0
                continue
            end
            for s = 2:length(group(g).rat(r).day(d).sleep) %after experience, so just sleep 2 on
                events = group(g).rat(r).day(d).sleep(s).popEv; %shorten variable name
                for i = 1:length(events)

                    r2 = events(i).r2;
                    if isnan(events(i).r2) || events(i).r2 < r2Thresh || events(i).propClose < propCloseThresh || events(i).maxJump > jumpThresh %it doesn't meet threshold
                        continue %to next event
                    end %r2 thresh

                    if r2 > min(r2ForPlot(:,g))
                        [~, storeInd] = min(r2ForPlot(:,g));
                        pxnForPlot{storeInd,g} = events(i).pxn;
                        r2ForPlot(storeInd,g) = r2;
                        slopeForPlot(storeInd,g) = events(i).slope*radCmConv;
                        ratNameForPlot{storeInd,g}  = group(g).rat(r).name;
                    end %use this event

                end %i - events
            end %sleep
        end %day

    end %rat
end %groupo

cd(saveDir)
keyboard

figtitle = 'ReplayFidelity_examples_highest';
%     figtitle = [figtitle '_' methodNames{detectMethod}];
figure('Name', figtitle, 'Position', [248 448 1349 532])
spMap = [1:5; 6:10];

for g = 1:2

    tmpPxnForPlot = pxnForPlot(:,g);
    tmpR2ForPlot = r2ForPlot(:,g);
    tmpSlopeForPlot = slopeForPlot(:,g);
    rNamesForPlot = ratNameForPlot(:,g);

    [~, sortOrd] = sort(tmpR2ForPlot, 'descend');

    for i = 1:5
        subplot(2,5,spMap(g,find(sortOrd == i)))

        maxTm = size(tmpPxnForPlot{i},2) * bayesStep + bayesWin/2;

        imagesc(0:1, 0:360, tmpPxnForPlot{i})
        axis xy
        axis square
        colormap hot
        caxis([0 0.2])

        xticks([0 1])
        xticklabels({'', num2str(maxTm)})
        if g == 2
            xlabel('Time (s)')
        end

        yticks([0 180 360])
        if find(sortOrd == i) == 1
            ylabel({groupNames{g}, 'Angular position (deg)'})
        end

        title({rNamesForPlot{i}; ['r^2 = ' num2str(round(tmpR2ForPlot(i),2)) ' | slope = ' num2str(round(tmpSlopeForPlot(i))) ' cm/s']})

    end %which map

end %group
cbr = colorbar;
set(cbr, 'Position', [.92 .60 .02 .32])
ylabel(cbr, 'Probability')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

% alt colors
figtitle = 'ReplayFidelity_examples_highest_altColors';
%     figtitle = [figtitle '_' methodNames{detectMethod}];
figure('Name', figtitle, 'Position', [248 448 1349 532])
spMap = [1:5; 6:10];

for g = 1:2

    custColMap = define_cust_color_map('White', cols{g}, 200, 0);

    tmpPxnForPlot = pxnForPlot(:,g);
    tmpR2ForPlot = r2ForPlot(:,g);
    tmpSlopeForPlot = slopeForPlot(:,g);
     rNamesForPlot = ratNameForPlot(:,g);

    [~, sortOrd] = sort(tmpR2ForPlot, 'descend');

    for i = 1:5
        subplot(2,5,spMap(g,find(sortOrd == i)))

        maxTm = size(tmpPxnForPlot{i},2) * bayesStep + bayesWin/2;

        imagesc(0:1, 0:360, tmpPxnForPlot{i})
        axis xy
        axis square
        %             colormap hot
        colormap(gca, custColMap)
        caxis([0 0.2])

        xticks([0 1])
        xticklabels({'', num2str(maxTm)})
        if g == 2
            xlabel('Time (s)')
        end

        yticks([0 180 360])
        if find(sortOrd == i) == 1
            ylabel({groupNames{g}, 'Angular position (deg)'})
        end

      title({rNamesForPlot{i}; ['r^2 = ' num2str(round(tmpR2ForPlot(i),2)) ' | slope = ' num2str(round(tmpSlopeForPlot(i))) ' cm/s']})

    end %which map

    cbr = colorbar;
    if g == 1
        set(cbr, 'Position', [.92 .60 .02 .32])
    else
        set(cbr, 'Position', [.92 .125 .02 .32])
    end %colorbar position
    ylabel(cbr, 'Probability')

end %group


if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end



end %function