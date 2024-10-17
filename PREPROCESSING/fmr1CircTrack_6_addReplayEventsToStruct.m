function group = fmr1CircTrack_6_addReplayEventsToStruct(group)
% function group = fmr1CircTrack_6_addReplayEventsToStruct(group)
%
% PURPOSE:
%   Detect replay events via two measures:
%       Detect population events from the binned firing rate.
%       Detect SWRs from the LFP.
%
% INPUT:
%   group data struct
%
% OUTPUT:
%   group data struct, with added for all beings and sleeps (with data):
%       popEvents: tms, pxn, r2, slope, bins2use, xAx, com, calphase,
%           maxJump, propClose
%
% MMD
% 7/2021
% Colgin Lab

%% INITIALIZE

runThresh = 5; %cm/s
Fs = 1000; %1 ms bin size
durBounds = [50/1000 2000/1000]; %50-2000 ms
pkCut = 3;
cellCrit = 5;

% BAYES PARAMETERS ARE CHANGED IN HELPER FUNCTION
% bayesWin = 20/1000; %20 ms time window, as in Hwaun & Colgin 2019
% bayesStep = 10/1000; %10 ms time step, as in Hwaun & Colgin 2019

stdCut = 5; %standard deviations above the mean for rip power cut
durThresh = 0.05; %for ripples


%% DO THE DATA

for g = 1:2
    fprintf('%s\n', group(g).name)
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d\n', r, length(group(g).rat))
        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day))

            if group(g).rat(r).day(d).decThresh == 0 %if meets decoding threshold
                for b = 1:length(group(g).rat(r).day(d).begin)
                    group(g).rat(r).day(d).begin(b).popEv = [];
                    group(g).rat(r).day(d).begin(b).rip = [];
                end %begin
                for s = 1:5
                    group(g).rat(r).day(d).sleep(s).popEv = [];
                    group(g).rat(r).day(d).sleep(s).rip = [];
                end %sleep
                continue %to next day
            end %if enough cells

            rateMaps = zeros(length(group(g).rat(r).day(d).xBegUnitInfo), length(group(g).rat(r).day(d).xBegUnitInfo(1).smRateMap));

            badU = [];
            uIDs = zeros(length(group(g).rat(r).day(d).xBegUnitInfo),2);
            for u = 1:length(group(g).rat(r).day(d).xBegUnitInfo)
                if group(g).rat(r).day(d).xBegUnitInfo(u).meetMin
                    rateMaps(u,:) = group(g).rat(r).day(d).xBegUnitInfo(u).smRateMap; %Smoothed ratemap
                    uIDs(u,:) = group(g).rat(r).day(d).xBegUnitInfo(u).ID;
                else
                    badU = [badU u]; %#ok
                end %if meets max
            end %units
            rateMaps(badU,:) = [];
            rateMaps(rateMaps==0) = 0.0001; %get rid of zeros because our Bayesian decoder can't handle 'em.
            uIDs(badU,:) = [];

            for s = 1:5
                if isempty(group(g).rat(r).day(d).sleep(s).unit)
                    continue
                end
                group(g).rat(r).day(d).sleep(s).popEv = []; %initialize
                coords = group(g).rat(r).day(d).sleep(s).coords;

                spkTmsByCell = cell(1,length(group(g).rat(r).day(d).xBegUnitInfo));
                for u = 1:length(group(g).rat(r).day(d).xBegUnitInfo)
                    if group(g).rat(r).day(d).xBegUnitInfo(u).meetMin
                        tmpVelSpks = group(g).rat(r).day(d).sleep(s).unit(u).spkTms; %sleep
                        spkTmsByCell{u} = tmpVelSpks;
                    end %meets min
                end %unit
                spkTmsByCell(badU) = [];

                [spkRstr, timeMap] = make_spike_raster(spkTmsByCell, [coords(1,1) coords(end,1)], Fs); %Fs = 1000 for 1 ms bin size
                [popEventInds, popEventTms] = find_population_events(spkRstr, Fs, timeMap, durBounds, pkCut, cellCrit);
                save([group(g).rat(r).day(d).sleep(s).dir '\PopEventOut_CA1'], 'spkRstr', 'timeMap', 'popEventInds', 'popEventTms', 'uIDs', 'rateMaps')

                events = get_events_info(popEventInds, popEventTms, spkRstr, rateMaps, uIDs);
                group(g).rat(r).day(d).sleep(s).popEv = events;

                if ~isfile([group(g).rat(r).day(d).sleep(s).dir '\RippleOut_CA1.mat'])
                    allLfp = cell(length(group(g).rat(r).day(d).tetNums),1);
                    for tt = 1:length(group(g).rat(r).day(d).tetNums)
                        lfpStruct = read_in_lfp([group(g).rat(r).day(d).sleep(s).dir '\CSC' num2str(group(g).rat(r).day(d).tetNums(tt)) '.ncs']);
                        allLfp{tt} = lfpStruct.data;
                    end %tetNums

                    [ripOnInds, ripOffInds] = DetectRipples_EH(allLfp, lfpStruct.Fs, stdCut, durThresh);
                    ripInds = cat(2, ripOnInds, ripOffInds);
                    ripTms = lfpStruct.ts(ripInds);

                    save([group(g).rat(r).day(d).sleep(s).dir '\RippleOut_CA1.mat'], 'ripInds', 'ripTms')
                end
            end %sleep

            for b = 1:length(group(g).rat(r).day(d).begin)
                group(g).rat(r).day(d).begin(b).popEv = []; %initialize
                %                 timesMat(bInds(b),:) = [group(g).rat(r).day(d).begin(b).radPos(1,1) group(g).rat(r).day(d).begin(b).radPos(end,1)];
                coords = group(g).rat(r).day(d).begin(b).coords;
                runSpeed = smooth_runspeed(get_runspeed(coords));

                spkTmsByCell = cell(1,length(group(g).rat(r).day(d).xBegUnitInfo));
                for u = 1:length(group(g).rat(r).day(d).xBegUnitInfo)
                    if group(g).rat(r).day(d).xBegUnitInfo(u).meetMin
                        spkInds = match(group(g).rat(r).day(d).begin(b).unit(u).spkTms, runSpeed(:,1));
                        spkSpds = runSpeed(spkInds,2);

                        tmpVelSpks = group(g).rat(r).day(d).begin(b).unit(u).spkTms(spkSpds < runThresh); %get velocity filtered spikes
                        spkTmsByCell{u} = tmpVelSpks;
                    end %meets min
                end %unit
                spkTmsByCell(badU) = [];

                [spkRstr, timeMap] = make_spike_raster(spkTmsByCell, [coords(1,1) coords(end,1)], Fs); %Fs = 1000 for 1 ms bin size
                if full(sum(spkRstr(:))) ~= 0
                    [popEventInds, popEventTms] = find_population_events(spkRstr, Fs, timeMap, durBounds, pkCut, cellCrit);
                    events = get_events_info(popEventInds, popEventTms, spkRstr, rateMaps, uIDs);
                    group(g).rat(r).day(d).begin(b).popEv = events;
                else
                    popEventInds = [];
                    popEventTms = [];
                end
                save([group(g).rat(r).day(d).begin(b).dir '\PopEventOut_CA1'], 'spkRstr', 'timeMap', 'popEventInds', 'popEventTms', 'uIDs')

                %because it's the begins, need to go through and get actual position
                for i = 1:length(group(g).rat(r).day(d).begin(b).popEv)
                    startInd = match(group(g).rat(r).day(d).begin(b).popEv(i).tms(1), group(g).rat(r).day(d).begin(b).radPos(:,1)); %match to times for pos  (different Fs)
                    endInd = match(group(g).rat(r).day(d).begin(b).popEv(i).tms(2), group(g).rat(r).day(d).begin(b).radPos(:,1));
                    group(g).rat(r).day(d).begin(b).popEv(i).actPos = rad2deg(wrapTo2Pi(circ_mean(deg2rad(group(g).rat(r).day(d).begin(b).radPos(startInd:endInd,2)))));
                end %i - events in this begin

            end %begin

        end %day
    end %rat
end %group

end %function

%% HELPER FUNCTINS

function events = get_events_info(popEventInds, popEventTms, spkRstr, rateMaps, uIDs)

bayesWin = 20/1000; %20 ms time window, as in Hwaun & Colgin 2019
bayesStep = 10/1000; %10 ms time step, as in Hwaun & Colgin 2019

Fs = 1000;

degBinCtrs = 2:4:360; %doesn't change across days/rats
radBinCtrs = deg2rad(degBinCtrs);

events = []; %initialize

for i = 1:size(popEventInds,1)

    onsetTm = popEventTms(i,1);
    onset = popEventInds(i,1);

    offsetTm = popEventTms(i,2);
    offset = popEventInds(i,2);

    events(i).inds = [onset offset];
    events(i).tms = [onsetTm offsetTm];

    evSpkRstr = spkRstr(:,onset:offset); %pull out spike raster across this event
    events(i).nCell = length(find((sum(evSpkRstr,2)>0))); %number of active cells in this event
    events(i).actCell = uIDs((sum(full(evSpkRstr),2)>0),:); %which cells were active

    pxn = BayesianDecoder(evSpkRstr, rateMaps, bayesWin, bayesStep, Fs);
    [nWin, winStartInds] = find_num_windows(size(evSpkRstr,2), bayesWin*Fs, bayesStep*Fs);
    if winStartInds(end)+bayesWin*Fs <= size(evSpkRstr,2)
        nWin = nWin + 1;
        winStartInds(end+1) = winStartInds(end)+bayesStep*Fs; %#ok
    end
    pxn(:,nWin+1:size(pxn,2)) = [];

    events(i).pxn = pxn;
    events(i).empBinProp =  sum(isnan(pxn(1,:)))/length(pxn(1,:));
    if events(i).empBinProp >= 0.2
        continue %it will get deleted later so skip the next stuff
    end %too many empty bins

    timeAx = 0:bayesStep:bayesStep*size(pxn,2)-bayesStep; %for calculating slope
    bins2use = find(~isnan(sum(pxn,1)));
    events(i).bins2use = bins2use;
    %         [events(i).r2, events(i).calphase, events(i).xAx, ~, events(i).slope] = Cir_reg(pxn, radBinCtrs', timeAx, bins2use);
    [events(i).xAx, events(i).com, events(i).calphase, events(i).slope, events(i).r2, p, events(i).maxJump, events(i).propClose] = circLinRegress_replay(pxn, radBinCtrs, timeAx, bins2use);

    if isnan(pxn(1,1)) || isnan(pxn(1,end))
        keyboard
    end %shouldn't happen

end %events

if ~isempty(events)
    pullCheck = vertcat(events(:).empBinProp);
    badEvs = pullCheck >= 0.2;
    events = events(~badEvs); %delete events that didn't meet criteria
end %there are events to check

end %helper function - event info