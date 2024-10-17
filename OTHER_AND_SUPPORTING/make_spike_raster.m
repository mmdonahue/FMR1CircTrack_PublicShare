function [spkRstr, timeMat] = make_spike_raster(spkTmsByCell, tms, Fs)
% function [spkRstr, timeMat] = make_spike_raster(spkTmsByCell, tms, Fs)
% 
% PURPOSE:
%   Generate a spike raster.
% 
% INPUT:
%   spkTmsByCell = cell array containing spike times for each cell
%   tms = start and end times for raster (will ignore spikes that occur
%       outside those times)
%   Fs = sampling frequency of raster
% 
% OUPUT:
%   spkRstr = sparse matrix that contains the time of each spike from
%       each cell in cell ID X sampling time point for each cell content
%   timeMat = times that match each time bin of the spike raster
% 
% MMD (based on EH code from Zheng et al. 2021)
% Colgin Lab
% 06/2022 - EDITED 01/2024


lenMat = round(diff(tms) * Fs); %apprx size of spike raster (in time dim)
timeMat = tms(1):1/Fs:tms(2);

%check that all spikes fall within the tms
for u = 1:length(spkTmsByCell)
    spkTmsByCell{u} = spkTmsByCell{u}(spkTmsByCell{u}>=tms(1) & spkTmsByCell{u}<=tms(2));
end %unit


numSpks = length(vertcat(spkTmsByCell{:})); %initialize for sparse
cellInd = zeros(numSpks,1);
spkTms = zeros(numSpks,1);

cellInd2 = [];
spkTms2 = [];
spkPerBin = [];

spkCntr = 0;
for u = 1:length(spkTmsByCell)
    spkRange = spkCntr+1:spkCntr+length(spkTmsByCell{u}); %which spikes
    cellInd(spkRange) = u;

    relStartSpks = spkTmsByCell{u} - tms(1); %times relative to start (in s)
    binSpks = round(relStartSpks * Fs); %bin according to Fs
    binSpks(binSpks==0) = 1;
    spkTms(spkRange) = binSpks;

    uBins = unique(binSpks); %unique bins
    cellInd2 = [cellInd2; repmat(u, length(uBins),1)];
    spkTms2 = [spkTms2; uBins];
    for b = 1:length(uBins)
            spkPerBin = [spkPerBin; sum(binSpks==uBins(b))];
    end %bins

    spkCntr = spkCntr + length(spkTmsByCell{u});
end %unit




% spkRstr = sparse(cellInd, spkTms, true(size(cellInd)), length(spkTmsByCell), lenMat);
spkRstr = sparse(cellInd2, spkTms2, spkPerBin, length(spkTmsByCell), lenMat);
% sparse(rowVal, columnVal, vectorVal, rowMat, colMat) 


end %function