function velSpkTms = velocity_filter_spikes(spkTms, smRunSpeed, minSpeed, maxSpeed)
% function velSpkTms = velocity_filter_spikes(spkTms, smRunSpeed, minSpeed, maxSpeed)
%
% PURPOSE:
%   Velocity filter the spike times to be within min speed-max speed range.
%
% INPUTS:
%   spkTms = spike times in seconds
%   smRunSpeed = run speed, smoothed, with (:,1) = time in seconds and
%       (:,2) = speed in cm/s
%   minSpeed = self-explanatory
%   maxSpeed = self-explanatory
%
% OUTPUT:
%   velSpkTms = spike times in seconds, velocity filtered within range
%
% MMD
% 09/2023
% Colgin Lab



spkInds = match(spkTms, smRunSpeed(:,1));
spkSpds = smRunSpeed(spkInds,2);

velSpkTms = spkTms(spkSpds < maxSpeed & spkSpds > minSpeed); %get velocity filtered spikes


end %function