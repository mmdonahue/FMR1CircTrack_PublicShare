function [p_resid, p_r2, r2, slope, xSpan, calphase, propPostClose, com] = get_seq_pVal_com(pxn,posBins,tBins,bins2use, varargin)
% function [p_resid, p_r2, r2, slope, xSpan, calphase, propPostClose, com] = get_seq_pVal(Pxn,xbins,tbins,bins2use,bound)
%
% PURPOSE:
%   Get the p-value of a sequence event as in Zheng et al. 2016, but adapted for the circle track.
%
% INPUT:
%   pxn = decoded probability distribution across the replay event (output of BayesianDecoder
%       function)
%   posBins = radial positions corresponding to each bin in pxn
%   tBins = time in seconds correspoding to each bin in pxn
%   bins2use = which bins to use (aka which bins have deocded info in them, aka which bins are not
%       NaN from the decoder)
%
% OUTPUT:
%   p_resid = p-value when comparing the shuffled residual sum of squares values
%   p_r2 = p-value when comparing to the shuffled r2 values
%
% Written by CZ 2015, modified for circle track by MMD 08/2024

% CZ: shuffle the position bins instead of time bins, which was uesd in
% Davison et al Neuron 2009 paper

%changed to use actual time and space dimensions rather than
%just bin numbers
%also added in com measeure and its difference from the predicted line as
%an additional measure of accuracy of prediction

bound = 8;
if nargin > 4
    bound = varargin{1};
end

num2shuf = 1000;

%% CALCULATE

com = NaN(1,size(pxn,2)); % initialize center of mass
maxPos = NaN(1,size(pxn,2)); % initialize pos of max decoding
for t = 1:length(bins2use) %across time bins to use
    bin = bins2use(t); %which time bin it is
    
    com(bin) = circ_mean(posBins, pxn(:,bin)); %weighted circular mean based on prob distribution

    [~, maxInd] = max(pxn(:,bin));
    maxPos(bin) = posBins(maxInd);
end %bins - i

com = wrapTo2Pi(com); %now we have the center of mass of each time bin

comUse = com(bins2use);
maxUse = maxPos(bins2use);
xAx = tBins(bins2use);

decUse = comUse;
para  = circ_lin_regress(xAx(~isnan(comUse)), decUse(~isnan(comUse)), bound); %CZ code
calphase = wrapTo2Pi(2*pi*para(1,1) * xAx+para(1,2));
slope = para(1,1)*2*pi;

xSpan = circ_dist(calphase(end), calphase(1));
if xSpan < 0 && slope > 0
    xSpan = 2*pi + xSpan;
elseif xSpan > 0 && slope < 0
    xSpan = 2*pi - xSpan;
end %dist is less than 0 but forward slope

if round(abs(para(1)),2) == bound
%     warning('Slope exceeds bounds')
    p_resid = NaN;
    p_r2 = NaN;
    r2 = NaN;
    slope = NaN;
    xSpan = NaN;
    calphase = NaN;
    propPostClose = NaN;
    return
end %maxes out

sumClose = 0; %initialize
sumAll = 0;
closeDef = deg2rad(20);
for t = 1:length(decUse)
  
    if isnan(com(bins2use(t)))
        continue
    end
      calVal = calphase(t);
    [~,binMin] = min(abs(circ_dist(posBins, (calVal - closeDef))));
    [~,binMax] = min(abs(circ_dist(posBins, (calVal + closeDef))));
    
    if binMax > binMin
        pullBins = pxn(binMin:binMax,bins2use(t));
    else
        pullBins = [pxn(binMin:end,bins2use(t)); pxn(1:binMax,bins2use(t))];
    end %whether cross 0
    
    sumClose = sumClose + sum(pullBins);
    sumAll = sumAll + sum(pxn(:,bins2use(t)));
end %time bins

propPostClose = sumClose / sumAll;

[r2, SSres] = circLinRegress_r2(decUse(~isnan(decUse)), calphase(~isnan(decUse)));


% shuffle the positions in each time bin
rsqdShuf = nan(num2shuf,1);
SSresidShuf = nan(num2shuf,1);

for j = 1:num2shuf %get other regression stats for shuffled versions
    ind = randsample(size(pxn,1),size(pxn,2));
    Pxn_Shuf = [];
shuffCom = NaN(1,size(pxn,2)); % initialize center of mass
for t = 1:length(bins2use) %across time bins to use
    bin = bins2use(t); %which time bin it is
      pdtemp = circshift(pxn(:,bin),ind(t))';
    
    shuffCom(bin) = circ_mean(posBins, pdtemp'); %weighted circular mean based on prob distribution

end %bins - t

comShuffUse = shuffCom(bins2use);
xAx = tBins(bins2use);

decUseShuff = comShuffUse;
paraShuff  = circ_lin_regress(xAx(~isnan(comUse)), decUse(~isnan(comUse)), bound); %CZ code
calphaseShuff = 2*pi*paraShuff(1,1) * xAx+paraShuff(1,2);

[rsqdShuf(j), SSresidShuf(j)] = circLinRegress_r2(decUseShuff(~isnan(decUseShuff)), calphaseShuff(~isnan(decUseShuff)));
end %j - shuffles

p_resid = length(find(SSresidShuf<SSres))/length(SSresidShuf);
p_r2 = length(find(rsqdShuf>r2))/length(rsqdShuf);
if isnan(r2)
%     keyboard
    p_r2 = 1;
end



end %function
