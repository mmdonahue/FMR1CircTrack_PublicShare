function spatialSparsity = get_spatial_sparsity(rateMap, timePerBin)
% function spatialSparsity = get_spatial_sparsity(rateMap, timePerBin)
%
% PURPOSE:
%   Computes spatial sparsity - "how diffuse the unit firing is in the
%   spatial domain ... if a unit fires equally all over the apparatus, the info 
%   per spike is 0 and sparsity is 1; if a unit fires evenly over one half of 
%   the apparatus but nevers fires in the other half, the info per spike is 1 bit 
%   and sparsity is 0.5." Reference: Jung, Wiener, & 
%   McNaughton (J Neurosci 1994). 
%
% INPUT:
%   rateMap = n x m matrix of firing rate by position bin
%   timePerBin = n x m matrix of time spent in each position bin
%
% OUTPUT:
%   spatialSparsity = sparsity of unit's place field
%
% Code is adapted from code written by AM, CZ, & MMD  for the Colgin Lab.
%
% PG Demetrovich
% Colgin Lab
% 10/2021

if size(rateMap) ~= size(timePerBin)
   warning('Both inputs must be the same size!') 
   keyboard
end

rateMap = reshape(rateMap,1,[]); %reshape to linearize
timePerBin = reshape(timePerBin,1,[]); 

nanInds = find(isnan(timePerBin)); %find any bins that the rat did not visit
rateMap(nanInds) = [];
timePerBin(nanInds) = [];

nanInds = find(isnan(rateMap)); %find any bins that the rat did not visit
rateMap(nanInds) = [];
timePerBin(nanInds) = [];

numBins = length(rateMap);
% meanFr = mean(rateMap); %same as lambda in formula
totTime = sum(timePerBin(:));

% spatInfo = 0;
pBinLambI   = [];
pBinLambISq = [];
for i = 1:numBins
        
    pBin = timePerBin(i)/totTime; %prob of being in i bin
    lambI = rateMap(i); %mean firing rate of unit in i bin
    pBinLambI(i)   = pBin*lambI;
    pBinLambISq(i) = pBin*(lambI.^2); 
    
%     binSpatInfo = pBin * (lambI/meanFr) * log2(lambI/meanFr);

%     if ~isnan(binSpatInfo)
%          spatInfo = spatInfo + binSpatInfo; %in AM and CZ code, ignore where bin Fr = 0, time per bin = 0
%     end %real number

end %bins

sparsityNum   = ( sum(pBinLambI) ).^2; % numerator in sparsity formula
sparsityDenom = sum(pBinLambISq);      % denominator in sparsity formula

spatialSparsity = sparsityNum / sparsityDenom;
end %function