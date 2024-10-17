function group = fmr1CircTrack_0_hardCodeTestData(group)
% function group = fmr1CircTrack_0_hardCodeTestData(group)
%
% Function adds the names, dates, reward locs, and theta tet to use to the
% group structure. It's just separated to this function vs it's parent because
% it takes less space.
%
% JBT / MMD
% 03/21
% Colgin Lab




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%                 KO RATS                 %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% RAT 316
group(2).rat(1).name = 'rat316';
group(2).rat(1).day(1).name = '2020-11-09'; %Day 2
group(2).rat(1).day(1).rewLocs = [0 180]; % E/W
group(2).rat(1).day(1).thetaTet = 4; %7 units, highest amp/cleanest theta

%% RAT 330

group(2).rat(2).name = 'rat330';
group(2).rat(2).day(1).name = '2021-03-04'; %Day 1
group(2).rat(2).day(1).rewLocs = 45; % NE (only 1 rew loc)
group(2).rat(2).day(1).thetaTet = 3; %8 units, tied for most units, cleanest theta

group(2).rat(2).day(2).name = '2021-03-08'; %Day 2 FOR MEG for Emma this is different
group(2).rat(2).day(2).rewLocs = [135 315]; % NW/SE
group(2).rat(2).day(2).thetaTet = 3; %selected by Emma

group(2).rat(2).day(3).name = '2021-03-15'; %Day 3 FOR MEG for Emma this is different
group(2).rat(2).day(3).rewLocs = [135 315]; % NW/SE
group(2).rat(2).day(3).thetaTet = 4; %selected by Emma

%% RAT 394

group(2).rat(3).name = 'rat394';
group(2).rat(3).day(1).name = '2022-10-15'; %Day 1
group(2).rat(3).day(1).rewLocs = 315; %green leg
group(2).rat(3).day(1).thetaTet = 8; %most units

%% RAT 395

group(2).rat(4).name = 'rat395';
group(2).rat(4).day(1).name = '2023-03-14'; %Day 1
group(2).rat(4).day(1).rewLocs = 180; %green leg
group(2).rat(4).day(1).thetaTet = 2; %only option

group(2).rat(4).name = 'rat395';
group(2).rat(4).day(2).name = '2023-03-17'; %Day 2
group(2).rat(4).day(2).rewLocs = 225; %yellow leg
group(2).rat(4).day(2).thetaTet = 2; %only option

group(2).rat(4).name = 'rat395';
group(2).rat(4).day(3).name = '2023-03-21'; %Day 3
group(2).rat(4).day(3).rewLocs = 135; %purple leg
group(2).rat(4).day(3).thetaTet = 17; %cleaner and more units

group(2).rat(4).name = 'rat395';
group(2).rat(4).day(4).name = '2023-03-29'; %Day 4
group(2).rat(4).day(4).rewLocs = 180; %green leg
group(2).rat(4).day(4).thetaTet = 2; %cleaner and more units

%% RAT 445

group(2).rat(5).name = 'rat445';
group(2).rat(5).day(1).name = '2024-04-14';
group(2).rat(5).day(1).rewLocs = 300; %green
group(2).rat(5).day(1).thetaTet = 2; %most units

group(2).rat(5).day(2).name = '2024-04-18'; 
group(2).rat(5).day(2).rewLocs = 120; %yellow
group(2).rat(5).day(2).thetaTet = 6; %second most units, larger amp theta

%% RAT 442

group(2).rat(6).name = 'rat442';
group(2).rat(6).day(1).name = '2024-05-02';
group(2).rat(6).day(1).rewLocs = [0 180]; %NW/SE or yellow/orange
group(2).rat(6).day(1).thetaTet = 8; %most units

group(2).rat(6).day(2).name = '2024-05-06';
group(2).rat(6).day(2).rewLocs = [45 225]; %NE/SW or purple/green
group(2).rat(6).day(2).thetaTet = 4; %most units

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%                  WT RATS                %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% RAT 326

group(1).rat(1).name = 'rat326';
group(1).rat(1).day(1).name = '2021-06-14A'; %Day 1
group(1).rat(1).day(1).rewLocs = [0 180]; % E/W
group(1).rat(1).day(1).thetaTet = [8]; %5 units, cleaner theta

group(1).rat(1).day(2).name = '2021-06-14B'; %Day 2
group(1).rat(1).day(2).rewLocs = [45 225]; %NE/SW AKA red and blue legs
group(1).rat(1).day(2).thetaTet = [7]; %4 units, but less noise

group(1).rat(1).day(3).name = '2021-06-15'; %Day 3
group(1).rat(1).day(3).rewLocs = [135 315]; %NW/SE AKA green and yellow legs
group(1).rat(1).day(3).thetaTet = [5]; %EXPLAIN HERE


%% RAT 334

group(1).rat(2).name = 'rat334';
group(1).rat(2).day(1).name = '2021-11-05'; %Day 1
group(1).rat(2).day(1).rewLocs = 135; % yellow leg
group(1).rat(2).day(1).thetaTet = 4; %most units

%% RAT 335 

group(1).rat(3).name = 'rat335';
group(1).rat(3).day(1).name = '2022-03-03'; %Day 1
group(1).rat(3).day(1).rewLocs = [155 310]; % 
group(1).rat(3).day(1).thetaTet = 2; %more units

group(1).rat(3).day(2).name = '2022-03-07'; %Day 2
group(1).rat(3).day(2).rewLocs = [135 275]; % 
group(1).rat(3).day(2).thetaTet = 4; %only 


%% RAT 392

group(1).rat(4).name = 'rat392';
group(1).rat(4).day(1).name = '2023-03-23'; %Day 1
group(1).rat(4).day(1).rewLocs = 315;  %blue
group(1).rat(4).day(1).thetaTet = 3; %2nd most units and pretty clean theta

group(1).rat(4).day(2).name = '2023-03-24'; %Day 1
group(1).rat(4).day(2).rewLocs = 0; %yellow
group(1).rat(4).day(2).thetaTet = 10; %2nd most units and pretty clean theta

%% RAT 416 

group(1).rat(5).name = 'rat416';
group(1).rat(5).day(1).name = '2023-10-27'; %Day 2
group(1).rat(5).day(1).rewLocs = 120; %yellow
group(1).rat(5).day(1).thetaTet = 6; %most units

group(1).rat(5).day(2).name = '2023-10-29'; %Day 2
group(1).rat(5).day(2).rewLocs = 225; %pink
group(1).rat(5).day(2).thetaTet = 7; %most units

%% RAT 418

group(1).rat(6).name = 'rat418';
group(1).rat(6).day(1).name = '2023-11-04'; %Day 1
group(1).rat(6).day(1).rewLocs = [0 180]; %NW/SE
group(1).rat(6).day(1).thetaTet = 18; %2nd most units and cleanest theta

group(1).rat(6).day(2).name = '2023-11-07'; %Day 1
group(1).rat(6).day(2).rewLocs = [45 225]; %NE/SW
group(1).rat(6).day(2).thetaTet = 18; %most units

end %function