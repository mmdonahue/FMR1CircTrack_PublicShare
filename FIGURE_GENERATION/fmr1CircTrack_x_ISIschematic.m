function fmr1CircTrack_x_ISIschematic(group)
% function fmr1CircTrack_x_ISIschematic(group)
%
% PURPOSE:
%   Create a schematic to clarify the difference between the first spike ISI and population ISI
%   measures.
%
% MMD
% 06/2024
% Colgin Lab

preEvTm = 0.02;

g = 1;
r = 6;
d = 1;
s = 4;
i = 60;

figtitle = 'ISI_schematic';
figure('Name', figtitle, 'Position', [680 106 560 872])

load([group(g).rat(r).day(d).sleep(s).dir '\PopEventOut_CA1.mat'], 'spkRstr', 'timeMap', 'rateMaps', 'uIDs')
inds = group(g).rat(r).day(d).sleep(s).popEv(i).inds;
actCell = group(g).rat(r).day(d).sleep(s).popEv(i).actCell;
tetNum = mode(actCell(:,1));

lfpStruct = read_in_lfp([group(g).rat(r).day(d).sleep(s).dir '\CSC' num2str(tetNum) '.ncs']);
%                    filtLFP{tInd} = filtfilt(B, A, lfpStruct.data);

subplot(3,1,1)

postStartTm = ceil(diff(group(g).rat(r).day(d).sleep(s).popEv(i).tms) * 100) ./ 100;

evLFP = lfpStruct.data(match(group(g).rat(r).day(d).sleep(s).popEv(i).tms, lfpStruct.ts)-...
    (preEvTm*Fs):match(group(g).rat(r).day(d).sleep(s).popEv(i).tms, lfpStruct.ts)+(postStartTm*Fs));
xVals = linspace(-preEvTm, postStartTm, length(evLFP));

plot(xVals, evLFP, 'Color', rgb('Black'))
xlim([-preEvTm postStartTm])
axis square
axis off

ax = gca;
line([-preEvTm (0.1-preEvTm)], [ax.YLim(1) ax.YLim(1)], 'Color', 'Black')
text(-preEvTm, ax.YLim(1)-50, '100 ms')

line([-preEvTm -preEvTm], [ax.YLim(1) ax.YLim(1)+200], 'Color', 'Black')
h = text(-(preEvTm + preEvTm/3), ax.YLim(1) , '200 uV');
set(h, 'Rotation', 90);

subplot(3,1,2)
[~, cellInds] = ismember(actCell, uIDs, 'rows');
actCellRstr = full(spkRstr(cellInds,inds(1)-(preEvTm/mean(diff(timeMap))):...
    inds(1)+(postStartTm/mean(diff(timeMap)))));

evTmMat = timeMap(inds(1)-(preEvTm/mean(diff(timeMap))):inds(1)+(postStartTm/mean(diff(timeMap))));
evTmMat = evTmMat - group(g).rat(r).day(d).sleep(s).popEv(i).tms(1);

[~, fsInd] = max(actCellRstr, [], 2);
[~, sortOrd] = sort(fsInd);

hold on
cMap = jet(length(actCell));


for u = 1:length(actCell)
    yVal = find(u==sortOrd);
    spkInds = find(actCellRstr(u,:));
    xVals = evTmMat(spkInds);
    yVals = repmat(yVal,length(spkInds),1);


    scatter(xVals, yVals, '|', 'MarkerEdgeColor', cMap(yVal,:))

end %unit

xlim([-preEvTm postStartTm])
ylim([-0.5 length(actCell)+0.5])
axis square
axis off

subplot(3,1,3)

hold on
allSpkTms = [];
fsSkTms = [];
for u = 1:length(actCell)
    yVal = find(u==sortOrd);
    spkInds = find(actCellRstr(u,:));
    xVals = evTmMat(spkInds);
    yVals = repmat(1,length(spkInds),1);
    allSpkTms = [allSpkTms xVals];

    scatter(xVals, yVals, '|', 'MarkerEdgeColor', cMap(yVal,:))

    scatter(xVals(1), 0.5,  '|', 'MarkerEdgeColor', cMap(yVal,:))
    fsSkTms = [fsSkTms xVals(1)];
end %unit

xlim([-preEvTm postStartTm])
ylim([0.5 1.25])
% ylim([-0.5 length(actCell)+0.5])
axis square
axis off

h = text((preEvTm + preEvTm/2), 1.1, ['Population ISI = ' num2str(mean(diff(allSpkTms)))]);
h = text((preEvTm + preEvTm/2), 0.6, ['First spike ISI = ' num2str(mean(diff(fsSkTms)))]);

cd('C:\Users\mdonahue\Documents\Lab\Generic figs')

saveas(gcf, figtitle, 'epsc');
saveas(gcf, figtitle, 'fig');
saveas(gcf, figtitle, 'png');


end %function