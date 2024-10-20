function pf = get_circtrack_pfs(rateMap, rmBinSz, minFr, minPfLen)
% function pf = get_circtrack_pfs(rateMap, rmBinSz, minPkFr, minPfLen)
%
% PURPOSE:
%   To find place-fields on the circular track. This differs from the linear track
%   algorithm because place-fields can wrap around from -pi to pi.
%   Place-fields here are defined as peaks above 2 SDs, with edges at 0.5 SDs, covering
%   at least 18 degrees by default(equivalent to 25 cms, pulled from Bieri et al., 2014).
%    NOTE: The Bieri et al. (2014) method didn't work here at all, producing fields that
%          were absurdly large, among other issues. Hence the switch to a z-score based threshold.
%
% INPUT:
%   rateMap = 1xn ratemap, NOT SMOOTHED. Output of 'get_ratemap_circtrack'
%   rmBinSz = bin size employed in 'get_ratemap_circtrack'
%     minFr = minimum average firing rate. If absent or empty, function defaults to 2.5 Hz (from Bieri et al., 2014)
%  minPfLen = minimum place-field length, in degrees. If absent or empty, function defaults to 18 degrees as noted above.
%
% OUTPUT:
%       pf = structure array of length #placeFields with subfields:
%            -  pf.pkFr = peak firing rate within the place-field
%            -  pf.pkPos = place-field peak firing rate position in degrees
%            -  pf.radPos = place-field position in degrees
%            -  pf.inds = place-field position in rate-map indices
%
% JB Trimper
% 10/2019
% Colgin Lab

pkZThresh = 2; %STDs  --  z-scorefor peak of place-field
edgeZThresh = 0.5; %STDs  --  z-score for edges of place-field

gWinStd = 8; %degrees, as in Zheng et al. 2021
gWinSz = gWinStd * 2;

if nargin < 3 || ~exist('minPkFr', 'var')
    minFr = 0.25; %Hz
end

if nargin < 4 || ~exist('minPfLen', 'var')
    minPfLen = 18; %degrees
    %               Bieri et al. (2014) used 3 contiguous 5 cm bins
    %               on the linear track (15 cm). The circumference of
    %               our 100 cm diameter circular track is 314.16 cm.
    %               Therefore, the circular track would be divided into
    %               ~21 bins using 15 cm bins. A 360 degree track
    %               divided by ~21 is ~17 degrees, and we're going
    %               to round that up to 18 degrees.
    %               In sum, a 15cm bin = ~18 degree bin on our track.
    %
end


% Centers of rate-map angular bins
rmBinAngs = edges_to_x_vals(linspace(0,360,(360/rmBinSz)+1));


pf = [];
pCntr = 1;
if mean(rateMap) > minFr %If minimum average firing rate criteria met
    
    zRm = zscore(rateMap); %zscore the ratemap
    pkInds = find_lfp_peaks(zRm'); %find high firing peaks
    pkInds(zRm(pkInds)<pkZThresh) = []; %keep only the peaks above threshold
    
    zSmRm = zscore(smooth_circtrack_ratemap(rateMap, rmBinSz, gWinSz, gWinStd)); %z-score the smoothed rate-map
    %                             We want peaks from the actual map so no spatial distortion,
    %                             but edges from the smoothed so small-scale fluctuations don't impact the field size
    
    if ~isempty(pkInds) %if there are any
        
        % Concatenate the final 3rd of the rate-map to the beginning and the 1st third of the rate-map
        %  to the end so we don't miss place-fields at the 'ends' of the circular track
        rmThird = ceil(length(rateMap)/3);
        extRmBinAngs = [rmBinAngs(end-rmThird+1:end) rmBinAngs rmBinAngs(1:rmThird)];
        zExtRm = [zSmRm(end-rmThird+1:end) zSmRm zSmRm(1:rmThird)];
        
        extPkInds = pkInds + rmThird; %align pk inds with extended ratemap
        
        for p = 1:length(pkInds) %for each peak
            
            %get the putative field edges (where the ratemap dropped below threshold)
            startField = find(zExtRm(1:extPkInds(p))<=edgeZThresh, 1, 'Last');
%             startField = startField + 1; %go up 1 to actually be above threshold
            endField = find(zExtRm(extPkInds(p):end)<=edgeZThresh, 1, 'First');
%             endField = endField + extPkInds(p) - 2; %go back 2 to align with threshold and re-align with whole map
endField = endField + extPkInds(p) - 1; %go back 2 to align with threshold and re-align with whole map
            
            %Get place field length in degrees (they're already in degrees, but circ_dist operates on radians, hence the conversions)
            pfLen = abs(rad2deg(circ_dist(deg2rad(extRmBinAngs(startField)), deg2rad(extRmBinAngs(endField)))));
            
            % If the place-field is big enough...
            if pfLen >= minPfLen
                
                pfInds = startField:endField; %indices of the place-field
                
                % Since the map was made artifically longer, if indices were above the rate-map size, wrap them back around
                %  i.e., account for the circular track
                overInds = find(pfInds>length(rateMap)+ rmThird);
                if ~isempty(overInds)
                    pfInds(overInds) = pfInds(overInds) - length(rateMap);
                end
                
                % Same as above, but if the indices were in the first part of the appended rate-map
                %  Confusing wording. Sorry.
                underInds = find(pfInds<rmThird+1);
                if ~isempty(underInds)
                    distFromRmStart = rmThird - pfInds(underInds);
                    pfInds(underInds) = length(rateMap) + rmThird - distFromRmStart;
                end
                
                % Subtract the size of that appended third to get indices scaled back to rate-map's original size
                pfInds = pfInds - rmThird;
                
                pkFr = rateMap(pkInds(p)); %Peak firing within the place-field
                
                % Save all of that for output
                pf(pCntr).pkFr = pkFr; %#ok
                pf(pCntr).pkPos = rmBinAngs(pkInds(p)); %#ok
                pf(pCntr).radPos = rmBinAngs(pfInds); %#ok
                pf(pCntr).inds = pfInds; %#ok
                
                % And advance the counter
                pCntr = pCntr + 1;
                
            end
            
        end %peak
        
        % If there is more than one place-field, make sure they don't overlap
        if pCntr > 2
            badPf = []; %for storing indices of overlapping fields
            for p1 = 1:length(pf)
                pf1Pos = pf(p1).radPos; %indices for the first pf
                for p2 = p1+1:length(pf)
                    pf2PkPos = pf(p2).pkPos; %firing peak of second pf
                    if ismember(pf2PkPos, pf1Pos)
                        
                        %The two place-fields will have different peaks, and we only want to keep that one that is greater
                        if pf(p1).pkFr >= pf(p2).pkFr
                            maxInd = p1;
                        else
                            maxInd = p2;
                        end
                        
                        pf(p1).pkFr = pf(maxInd).pkFr; %#ok
                        pf(p1).pkPos = pf(maxInd).pkPos; %#ok
                        
                        badPf = [badPf p2]; %#ok - if peak of second pf was within field of 1st pf, keep index for removal
                    end
                end
            end
            if ~isempty(badPf)
                pf(badPf) = []; %remove bad ones
            end
        end
        
        
        
    end %if peaks>pkZThresh found
end %if min avg firing rate criteria met

end %fnctn