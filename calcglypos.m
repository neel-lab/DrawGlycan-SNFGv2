function [plotinfoout,specialoptions] = calcglypos(sgpseq,options,specialoptions)
% CALCGLYPOS: calculate monosaccharide positions of glycan
%
% Syntax:
% [plotinfoout,specialoptions] = calcglypos(sgpseq,options,specialoptions)
%
% Input:
% sgpseq: string, the glycans to be analyzed.
% options: structure, contains global options for drawing glycans. See
% DRAWGLYCAN for details.
% specialoptions: structure, pre-compiled local options. Each field is a
% type of local option, which contains a n x 2 cell array, 1st column is
% the index number of monosac. that's been modified, 2nd column is the
% detail of the modification.
%
% Output:
% plotinfoout: structure containing the following fields:
% "mspos": n x 2 numerical matrices, monosaccharide positions in the glycan.
% n equals to the number of monosaccharides in the glycan.
% "bondmap": n x n numerical matrices, the linkage information of monosac. in glycan.
% "alllinkout": n x 1 cell array of strings, the glycosidic bond information of glycans.
% "allms": n x 1 cell array of strings, the name of the monosaccharides in glycan.
% "directionseq": n x 1 numerical matrices, the orientation of monosaccharides.
% FIELDS BELOW APPEARS ONLY WHEN INPUT "SPECIALOPTIONS" IS NOT EMPTY
% "bonddescription": modification of individual monosaccharides, such as
% -U, -D, and bond modifications, such as -ZIG, -BOLD.
% other fields: all fields of input "specialoptions" will be copied here.
%
% Note:
% Solo perpendicular monosaccharides will be  plotted up,
%e.g. the top vertex of a solitary Fuc will point up
% whatever the orientation of the glycan is.
%
% Example:
% N/A. Run examples in DRAWGLYCAN and set breakpoints.
%
% Children function:
% DRAWRAWTREE, BRANCHEQUALIZER, SMARTPMDISPENSER,
% BONEADDPM, ADDSTUBPM

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2019, Research Foundation for State University of New York. All rights reserved
%


sortbranches = false;
if ischar(options.sortbranches)
    if strcmpi(options.sortbranches,'yes') || strcmpi(options.sortbranches,'true')
        sortbranches = true;
    end
elseif isnumeric(options.sortbranches)
    if options.sortbranches > 0
        sortbranches = true;
    end
elseif islogical(options.sortbranches)
    sortbranches = options.sortbranches;
end
if sortbranches
    [sgpseq,msindfix] = branchreform(sgpseq);
    if ~isempty(specialoptions)
        fldnms = fieldnames(specialoptions);
        for i = 1:length(fldnms)
            for j = 1:size(specialoptions.(fldnms{i}),1)
                if ~any(strfind(fldnms{i},'pep'))
                    specialoptions.(fldnms{i}){j,2} = specialoptions.(fldnms{i}){j,2} + msindfix(specialoptions.(fldnms{i}){j,2});
                end
            end
        end
    end
end
alllinkout = regexp(sgpseq,DrawGlycanPara.regexp_monosaclinkage,'match');
for i = 1:length(alllinkout)
    alllinkout{i} = alllinkout{i}(2:end-1);
end
%% Start calculation

%% sep. monosac. and bond
thisgly = sgpseq;
thisgly = strrep(thisgly,'[','');
thisgly = strrep(thisgly,']','');
indtemp = strfind(thisgly,'{');
levelindex = zeros(2,length(thisgly));
indtemp2 = strfind(thisgly,'}');
if ~ismember(1,indtemp)
    error('Did you forget something?');
end
levelindex(1,indtemp) = 1;
levelindex(1,indtemp2) = -1;
levelindex(2,1) = 1;
for j = 2:size(levelindex,2)
    levelindex(2,j) = levelindex(2,j-1) + levelindex(1,j);
end
wholestr = zeros(1,length(thisgly));
wholestr(indtemp) = 1;
wholestr(indtemp2) = 1;
readind = 1;
allms = cell(sum(wholestr)/2,1);
writeind = 1;
tempstr = '';
while readind <= length(wholestr)
    if wholestr(readind) ~= 1  % gather character to form the monosac. string
        tempstr = [tempstr,thisgly(readind)];
    else
        if ~isempty(tempstr)
            bondinfo = strfind(tempstr,')');
            if ~isempty(bondinfo)  % using parenthesis
                allms{writeind} = tempstr(max(bondinfo) + 1:end);
            else  % No parenthesis means there's only monosac.
                allms{writeind} = tempstr;
            end
            tempstr = '';
            writeind = writeind + 1;
        end
    end
    readind = readind + 1;
end
%% calculate distance of each monosac
letterindex = zeros(1,length(thisgly));
letterindex(regexp(thisgly,'[^{}]')) = 1;
distance = letterindex.*levelindex(2,:);
distance = distance(indtemp+1);  % all monosac's

%% build adjacency matrix "bondmap"
bondmap = zeros(length(distance));
readind = 1;
while readind < length(distance)
    if (distance(readind + 1) - distance(readind)) == 1
        bondmap(readind,readind + 1) = 1;  % consecutive numbers indicate bond
        readind = readind + 1;
    elseif (distance(readind + 1) - distance(readind)) < 1  % if chain is broken, go back to find its fork point
        thisind = distance(readind + 1);  % where it's broken
        itsforkpt = find(distance(1:readind) == thisind - 1,1,'last');  % where is the fork point
        bondmap(itsforkpt,readind + 1) = 1;  % mark this bond
        readind = readind + 1;  % keep going on
    end
end

ispm = find(ismember(lower(allms),lower(options.perpendicularmonosac)));
%% new version: monosac. is PM only when it's at the terminal of a branch
isterminal = zeros(size(ispm));
for i = 1:length(ispm)
    isterminal(i) = ~any(bondmap(ispm(i),:));
end
ispm = ispm(logical(isterminal));
%% new version end

ispm = [ispm,zeros(size(ispm))];  % 1: perpendicular to the right (default)  2: to the left. Assuming parent structure points upwards
%% customized perpendicular monosac.
ppmopt = {'P','PL','PR'};
ppmlr = [0 1 2];
% specialpm = [];
if ~isempty(specialoptions)
    spoptfld = fieldnames(specialoptions);
    if any(ismember(ppmopt,upper(spoptfld)))
        for i = 1:length(ppmopt)
            if ismember(ppmopt{i},upper(spoptfld))
                specialpm = specialoptions.(ppmopt{i});
                specialpm = cell2mat(specialpm(:,2));
                ispm = [[specialpm(:),ones(length(specialpm),1)*ppmlr(i)];ispm];
            end
        end
    end
end

if ~isempty(ispm)
    [~,ind,~] = unique(ispm(:,1));
    ispm = ispm(ind,:);
end
% ispm = [ispm;specialpm(:)];
directionseq = zeros(size(distance));


%% Decide working mode
if isempty(ispm)
    drawglycanmode = 1;  % bone only
else
    if ismember(1,ispm(:,1))
        drawglycanmode = 2;  % PMS only
    else
        drawglycanmode = 3;  % bone+PMS
    end
end

%% Handling 3 different situations
switch drawglycanmode
    case 1
        mspos = drawrawtree(bondmap,distance,[]);
        mspos = branchequalizer(mspos,bondmap);
        
    case 2
        mspos = drawrawtree(bondmap,distance,[]);
        mspos = branchequalizer(mspos,bondmap);
        
    case 3
        bonebondmap = bondmap;  % bone part
        bonedistance = distance;
        pmchildlist = {};  % Each element in 'pmschildlist' is a sub-tree connected to the main structure,
        % the first monosac in each element is 'thispms'
        while ~isempty(ispm)  % retrieve the monosac index of PMS containing sub-trees
            tobetracked = ispm(find(ispm(:,1),1,'first'),1);
            temppmchild = [];
            while ~isempty(tobetracked)
                if any(bondmap(tobetracked(1),:))
                    tobetracked = [tobetracked find(bondmap(tobetracked(1),:))];
                end
                temppmchild = [temppmchild;tobetracked(1)];
                tobetracked(1) = [];
            end
            [~,pmchildind] = ismember(ispm(:,1),temppmchild);
            pmchildind = pmchildind(pmchildind > 0);
            tempispm = ispm(pmchildind,:);
            temppmchild = [temppmchild,ones(size(temppmchild))*tempispm(1,2)];
            pmchildlist = [pmchildlist;temppmchild];
            [~,ispmind] = setdiff(ispm(:,1),temppmchild);
            ispm = ispm(ispmind,:);
        end
        parentms = zeros(size(pmchildlist));
        allpmchild = [];  % this is used later to draw the bone structure of glycan
        for j = 1:size(pmchildlist,1)
            parentms(j,1) = find(bondmap(:,pmchildlist{j}(1)));
            allpmchild = [allpmchild;pmchildlist{j}];
        end
        [parentmsnum,~,parentmsind] = unique(parentms);
        pmgroup = cell(max(parentmsind),2);  % The number of elements in 'pmsgroup' equals the number of
        for j = 1:size(pmgroup,1)
            pmgroup{j,1} = parentmsnum(j);  % 1st column contains the parent monosac
            pmgroup{j,2} = pmchildlist(parentmsind == j);  % 2nd column contains children PMS
        end
        bonebondmap(:,allpmchild(:,1)) = 0;
        bonebondmap(allpmchild(:,1),:) = 0;
        bonedistance(allpmchild(:,1)) = 0;
        tempbonebondmap = bonebondmap;
        tempbonebondmap(allpmchild(:,1),:) = [];
        tempbonebondmap(:,allpmchild(:,1)) = [];
        mspos = drawrawtree(bondmap,distance,allpmchild(:,1));  % Here mspos contains all pos info of bone monosac
        bonemspos = mspos;
        bonemspos(allpmchild(:,1),:) = [];
        bonemspos = branchequalizer(bonemspos,tempbonebondmap);
        mspos(setdiff(1:length(distance),allpmchild(:,1)),:) = bonemspos;
        %% start dealing with PM
        [pmposstorage,pmindstorage,directionseq] = smartPMdispenser(pmgroup,bondmap,distance,directionseq);
        allfork = find(sum(bonebondmap,2) > 1);
        if ~isempty(allfork)
            % BONEADDPM adds perp. monosac. to branches of the glycan tree structure
            [mspos,safezone,PMroot] = boneaddPM(allfork,mspos,bonebondmap,bondmap,pmgroup,pmindstorage,pmposstorage,allpmchild);
            stub = min(find(sum(bonebondmap,2) > 1,1,'first')) - 1;  % find the operating region for the next step, beyond this point everything stays at where they are now.
            % ADDSTUBPM adds perp. monosac. to the "trunk" of the tree structure
            mspos = addstubPM(mspos,bondmap,stub,bonebondmap,PMroot,pmposstorage,pmindstorage,safezone,directionseq);
        else
            safezone = find(bonedistance);
            stub = safezone(end);
            PMroot = cell2mat(pmgroup(:,1));
            [mspos,directionseq] = addstubPM(mspos,bondmap,stub,bonebondmap,PMroot,pmposstorage,pmindstorage,safezone,directionseq);
        end
end
mspos(:,2) = -mspos(:,2);  % glycan has been drawn upside down, flip it up.
mspos = mspos - repmat(mspos(1,:),size(mspos,1),1);  % normalize all coordinates


%% Option: orientation
if strcmpi(options.orientation,'left')
    mspos = -mspos;
    mspos = mspos - repmat(mspos(1,:),size(mspos,1),1);
elseif strcmpi(options.orientation,'up')
    tmp = mspos(:,1);
    mspos(:,1) = -mspos(:,2);
    mspos(:,2) = tmp;
    mspos = mspos - repmat(mspos(1,:),size(mspos,1),1);
elseif strcmpi(options.orientation,'down')
    tmp = mspos(:,1);
    mspos(:,1) = mspos(:,2);
    mspos(:,2) = -tmp;
    mspos = mspos - repmat(mspos(1,:),size(mspos,1),1);
end
plotinfoout.mspos = mspos;
plotinfoout.bondmap = bondmap;
plotinfoout.alllinkout = alllinkout(:);
plotinfoout.allms = allms;
plotinfoout.directionseq = directionseq(:);
%% pass info from specialoptions to plotinfoout
if ~isempty(specialoptions)
    spfieldnames = fieldnames(specialoptions);
    bonddesname = spfieldnames(ismember(spfieldnames,[DrawGlycanPara.intglybondinfo,DrawGlycanPara.intglymodinfo]));
    bonddescription = {};
    for i = 1:length(bonddesname)
        tempbonddes = specialoptions.(bonddesname{i});
        if ~isempty(tempbonddes)
            bonddescription = [bonddescription;tempbonddes(:,2),repmat(bonddesname(i),size(tempbonddes,1),1),tempbonddes(:,1)];
        end
    end
    plotinfoout.bonddescription = bonddescription;
    spfieldnames = setdiff(spfieldnames,bonddesname);
    for i = 1:length(spfieldnames)
        if strcmpi(spfieldnames{i},'PEPC') || strcmpi(spfieldnames{i},'PEPN')
            plotinfoout.(spfieldnames{i}(4)) = specialoptions.(spfieldnames{i});
        else
            plotinfoout.(spfieldnames{i}) = specialoptions.(spfieldnames{i});
        end
    end
end
end

