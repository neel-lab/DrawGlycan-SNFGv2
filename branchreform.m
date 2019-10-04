function [reformed,msposfix] = branchreform(thisgly)
% BRANCHREFORM: Rearrange branches in glycan structure according to linkage
% information.
%
% Syntax:
% [reformed,msposfix] = branchreform(thisgly)
%
% Input:
% thisgly: string, SGP2.0 sequence of glycan.
%
% Output:
% reformed: string, condensed IUPAC sequence of glycan after rearrangement.
% msposfix: 1 x n numerical array, the position shift of each
% monosaccharide. For example, a glycan "{n{s}{h{s}}" was rearranged into
% "{n{h{s}}{s}}, the first "n" stays, so msposfix(1) = 0, first "{s}" was
% moved to the end, position moved from 2 to 4, so msposfix(2) = 2.
% "{h{s}}" changed from 3rd and 4th to 2nd and 3rd, so msposfix(3) and
% msposfix(4) both = -1.
% If no rearrangement happens, this variable is [0,0,0,...].
%
% Note:
% Glycan structure is sorted by ascending order, first anomeric carbon
% (alpha, beta,...) then hydroxyl group(1, 2, 3,...).
%
% Example:
% reformed,msposfix]=branchreform('{GalNAc(??-?){Neu5Ac(a2-6)}{Gal(b1-3){Neu5Ac(a2-6)}}}')
% 
% reformed =
% 
%     '{GalNAc(??-?){Gal(b1-3){Neu5Ac(a2-6)}}{Neu5Ac(a2-6)}}'
% 
% 
% msposfix =
% 
%      0     2    -1    -1
%
% Children function:
% N/A
%

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2019, Research Foundation for State University of New York. All rights reserved
%


thisgly = strrep(thisgly,'[','');
thisgly = strrep(thisgly,']','');
thisglynocurl = thisgly;
[optvalstart,optvalend] = regexp(thisgly,'["''].*?["'']','start','end');
for j = 1:length(optvalstart)
    thisglynocurl(optvalstart(j):optvalend(j)) = ['"',repmat('?',1,optvalend(j)-optvalstart(j)-1),'"'];
end
indtemp = strfind(thisglynocurl,'{');
levelindex = zeros(2,length(thisgly));
indtemp2 = strfind(thisglynocurl,'}');
if ~ismember(1,indtemp)
    error('Does glycan starts with a "{"?');
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
allbond = cell(sum(wholestr)/2,1);
writeind = 1;
tempstr = '';
readind = 1;
while readind <= length(wholestr)
    if wholestr(readind) ~= 1  % gather character to form the monosac. string
        tempstr = [tempstr,thisgly(readind)];
    else
        if ~isempty(tempstr)
            bondinfo = regexp(tempstr,'[a-z?][\d?]+-[\d?]+','match');
            allbond{writeind} = bondinfo{1};
            tempstr = '';
            writeind = writeind + 1;
        end
    end
    readind = readind + 1;
end

%% calculate distance of each monosac
letterindex = zeros(1,length(thisgly));
letterindex(regexp(thisglynocurl,'[^{}]')) = 1;
distance = letterindex.*levelindex(2,:);
distance = distance(indtemp+1);  % all monosac's
if length(indtemp) > 1
    for j = 2:length(indtemp)
        letterindex(indtemp(j):end) = letterindex(indtemp(j):end)/(j-1)*j;
    end
end
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
isfork = find(sum(bondmap,2) > 1);
msseq = 1:size(bondmap,1);
if any(isfork)
    for j = length(isfork):-1:1
        forkchildren = find(bondmap(isfork(j),:));
        branches = cell(size(forkchildren));
        bonds = cell(length(forkchildren),1);
        branchpos = zeros(length(forkchildren),2);
        for k = 1:length(forkchildren)
            branches{k} = glytreetracker(bondmap,forkchildren(k),[],'down');
            branchpos(k,1) = min(find(letterindex == forkchildren(k))) - 1;
            branchpos(k,2) = find(levelindex(2,branchpos(k,1):end) == levelindex(2,branchpos(k,1))-1,1) + branchpos(k,1)-1;
            tempbond = allbond{forkchildren(k)};
            carbon = regexp(tempbond,'[a-z?]','match');
            nrnum = regexp(tempbond,'[0-9?]+','match');
            bonds{k,1} = carbon{1};
            bonds{k,2} = nrnum{end};
        end
        tempind = 1:size(bonds,1);
        branchseq = zeros(size(tempind));
        ind1 = ones(size(bonds,1),1);
        writeind = 1;
        for k = 1:max(ind1)
            [~,~,ind2] = unique(bonds(ind1 == k,2));
            indperletter = tempind(ind1 == k);
            for l = 1:max(ind2)
                indpernumber = indperletter(ind2 == l);
                branchseq(writeind:writeind + length(indpernumber) - 1) = indpernumber;
                writeind = writeind + length(indpernumber);
            end
        end
        newbranches = branches(branchseq);
        thisglyhead = thisgly(1:min(min(branchpos))-1);
        thisglytail = thisgly(max(max(branchpos))+1:end);
        thisglymiddle = thisgly(min(min(branchpos)):max(max(branchpos)));
        letterindhead = letterindex(1:min(min(branchpos))-1);
        letterindtail = letterindex(max(max(branchpos))+1:end);
        letterindmiddle = letterindex(min(min(branchpos)):max(max(branchpos)));
        lvlindhead = levelindex(:,1:min(min(branchpos))-1);
        lvlindtail = levelindex(:,max(max(branchpos))+1:end);
        lvlindmiddle = levelindex(:,min(min(branchpos)):max(max(branchpos)));
        msseqhead = msseq(1:min(cellfun(@min,newbranches))-1);
        msseqtail = msseq(max(cellfun(@max,newbranches))+1:end);
        msseqmiddle = [];
        for k = 1:length(newbranches)
            msseqmiddle = [msseqmiddle,msseq(newbranches{k})];
        end
        msseq = [msseqhead,msseqmiddle,msseqtail];
        minmsind = cellfun(@min,branches);
        maxmsind = cellfun(@max,branches);
        bondmaphead = bondmap(1:min(minmsind)-1,:);
        bondmaptail = bondmap(max(maxmsind)+1:end,:);
        branchpos = branchpos-min(min(branchpos))+1;
        newglymiddle = [];
        newletterindmiddle = [];
        newlvlindmiddle = [];
        newbondmapmiddle = [];
        newbranchlength = cellfun(@length,branches(branchseq));
        newforkchildren = forkchildren;
        for k = 2:length(forkchildren)
            newforkchildren(k) = newforkchildren(k-1) + newbranchlength(k-1);
        end
        bondmapforkline = zeros(1,size(bondmap,2));
        bondmapforkline(newforkchildren) = 1;
        bondmaphead(end,:) = bondmapforkline;
        bondmapfix = newforkchildren - forkchildren(branchseq);
        
        for k = 1:length(branchseq)
            newglymiddle = [newglymiddle,thisglymiddle(branchpos(branchseq(k),1):branchpos(branchseq(k),2))];
            newletterindmiddle = [newletterindmiddle,letterindmiddle(branchpos(branchseq(k),1):branchpos(branchseq(k),2))];
            newlvlindmiddle = [newlvlindmiddle,lvlindmiddle(:,branchpos(branchseq(k),1):branchpos(branchseq(k),2))];
            tempbondmapmiddle = bondmap(branches{branchseq(k)},:);
            if bondmapfix(k) < 0
                tempbondmapmiddle = [tempbondmapmiddle(:,-bondmapfix(k)+1:end),zeros(size(tempbondmapmiddle,1),-bondmapfix(k))];
            elseif bondmapfix(k) > 0
                tempbondmapmiddle = [zeros(size(tempbondmapmiddle,1),bondmapfix(k)),tempbondmapmiddle(:,1:size(tempbondmapmiddle,2)-bondmapfix(k))];
            end
            newbondmapmiddle = [newbondmapmiddle;tempbondmapmiddle];
        end
        thisgly = [thisglyhead,newglymiddle,thisglytail];
        letterindex = [letterindhead,newletterindmiddle,letterindtail];
        levelindex = [lvlindhead,newlvlindmiddle,lvlindtail];
        bondmap = [bondmaphead;newbondmapmiddle;bondmaptail];
    end
end
reformed = thisgly;
[~,ind] = sort(msseq);
msposfix = ind - (1:length(msseq));
end