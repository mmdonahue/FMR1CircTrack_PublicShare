function group = fmr1CircTrack_1_buildDataStruct
% function group = fmr1CircTrack_1_buildDataStruct(dataDir)
%
% PURPOSE:
%  Function to build the main data structure for circle track FMR1 project.
%
%
% OUTPUT:
%  group = structure with data for each group (FXS vs WT)
%
% JBT / MMD
% 11/29/2020
% Colgin Lab


dataDir = 'E:\FMR1_CIRCTRACK\RAW_DATA';

group(1).name = 'WT';
group(2).name = 'KO';

group = fmr1CircTrack_0_hardCodeTestData(group); %Add names, dates, theta tet, reward locs

curDir = pwd;

ttList = 'CA1.txt';

for g = 1:2
    fprintf('%s\n', group(g).name)
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);
        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day));

            
            cd([dataDir '\' group(g).name '\' group(g).rat(r).name '\' group(g).rat(r).day(d).name]);
%             fid = fopen('CA1.txt');
            fid = fopen(ttList);
            
            tmp = textscan(fid, '%s', 'delimiter', '\n');
            unitList = tmp{1};
            fclose(fid);
            
            for b = 1:4
                fprintf('\t\t\tBegin %d\n', b);
                
                group(g).rat(r).day(d).begin(b).dir = [dataDir '\' group(g).name '\' group(g).rat(r).name '\' group(g).rat(r).day(d).name '\begin' num2str(b)];
                try
                cd(group(g).rat(r).day(d).begin(b).dir);
                catch
                    group(g).rat(r).day(d).begin(b) = []; %no begin 4 for one rat
                    continue
                end
                [post,posy,posx] = LoadCircPos('VT1.nvt');
                
                radPos = circpos(posx,posy); %radial position
                radPos = [post' radPos'];
                
                coords = zeros(length(post), 3);
                coords(:,1) = post;
                coords(:,2) = posx;
                coords(:,3) = posy;
                
                group(g).rat(r).day(d).begin(b).radPos = radPos;
                group(g).rat(r).day(d).begin(b).coords = coords;
                
                tetNums = [];

                for u = 1:length(unitList)
                         hyphInd = strfind(unitList{u}, '_');
                    tetNum = str2double(unitList{u}(3:hyphInd-1));
                    extInd = strfind(unitList{u}, '.t');
                    clustNum = str2double(unitList{u}(hyphInd+1:extInd-1));
                    try
                        if isempty(~strfind(unitList{u}, '64'))
                            spkTms = Readtfile(unitList{u});
                        else
                            spkTms = Readtfile(unitList{u}, 'uint64');
                        end
                    catch
                        if ~isfile(unitList{u}) && isfile(['TT' num2str(tetNum) '.clusters'])
                            spkTms = [];
                        else
                            keyboard
                        end
                    end

                    spkTms = spkTms ./ 10^4; %convert to seconds
                    
                    group(g).rat(r).day(d).begin(b).unit(u).ID = [tetNum clustNum];
                    group(g).rat(r).day(d).begin(b).unit(u).spkTms = spkTms;
                    
                    tetNums = [tetNums tetNum]; %#ok

                end %unit

                fprintf('\t\t\t\t%d Units Attached\n', length(unitList));
                
                tetNums = unique(tetNums);
                group(g).rat(r).day(d).tetNums = tetNums;
                
                for tt = 1:length(tetNums)
                    cscFn = ['CSC' num2str(tetNums(tt)) '.ncs'];
                    if ~isfile([cscFn(1:end-4) '_narrowThetaLfp.mat'])
                        fprintf('\t\t\t\tFiltering LFP for Tetrode #%d\n', tetNums(tt));
                        lfpStruct = read_in_lfp(cscFn);
                        filtLfp = filter_lfp(lfpStruct, 6, 10);
                        save([cscFn(1:end-4) '_narrowThetaLfp'], 'filtLfp');
                    end %if the filtered lfp isn't saved yet
                    if ~isfile([cscFn(1:end-4) '_deltaLfp.mat'])
      fprintf('\t\t\t\tFiltering LFP for Tetrode #%d\n', tetNums(tt));
                        lfpStruct = read_in_lfp(cscFn);
                        filtLfp = filter_lfp(lfpStruct, 2, 4);
                        save([cscFn(1:end-4) '_narrowThetaLfp'], 'filtLfp');
                    end
                end %tetrodes that had cells
                
                
            end %begin
            
        end %day
    end
end

cd(curDir);

end %fnctn
