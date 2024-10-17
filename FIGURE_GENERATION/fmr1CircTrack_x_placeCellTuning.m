function fmr1CircTrack_x_placeCellTuning(group)
% function fmr1CircTrack_x_placeCellTuning(group)
% 
% PURPOSE:
%   Plot tuning curves for CA1 place cells on the circle track.
% 
% MMD
% 08/2024
% Colgin Lab


groupNames = {'WT', 'FXS'};
spMap = [1:4; 5:8; 9:12; 13:16; 17:20];


saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\PLACE_CELLS\tuningCurves';

g = 1; r = 5; d = 2; tt = 12; % WT example
% g = 2; r = 6; d = 1; tt = 5; % FXS example

figtitle = [groupNames{g} '_' group(g).rat(r).name];
figure('Name', figtitle, 'Position', [336 108 1103 847])

findUs = vertcat(group(g).rat(r).day(d).xBegUnitInfo(:).ID);
findUs = findUs(:,1);
useUs = find(findUs == tt);

for uInd = 1:5
    for b = 1:4
        subplot(5, 4, spMap(uInd,b))

        u = useUs(uInd);
        plot(2:4:360, group(g).rat(r).day(d).begin(b).unit(u).smRateMap, 'Color', 'k')
        axis square
        xlim([0 360])
        xticks([0 360])

        if b == 1
            ylabel({['TT' num2str(group(g).rat(r).day(d).begin(b).unit(u).ID(1)) '\_' num2str(group(g).rat(r).day(d).begin(b).unit(u).ID(2))]...
                'Firing rate (Hz)'})
        end
        if uInd == 5
            xlabel('Position (deg)')
        end

    end %begin

end %unitInd

for uInd = 1:5
    same_axes(spMap(uInd,:));
end %uInd

cd(saveDir)

saveas(gcf, figtitle, 'fig');
saveas(gcf, figtitle, 'png');
saveas(gcf, figtitle, 'epsc');


end %function