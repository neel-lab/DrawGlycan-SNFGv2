function [mspos,safezoneind,PMroot] = boneaddPM(allfork,mspos,bonebondmap,bondmap,pmgroup,pmindstorage,pmposstorage,allpmchild)
% BONEADDPM: collate bone structure with perpendicular monosac. subtrees to
% create full glycan
%
% Syntax:
% [mspos,safezoneind,PMSroot] = boneaddPM(allfork,mspos,bonebondmap,bondmap,pmgroup,pmindstorage,pmposstorage,allpmchild)
%
% Input:
% allfork: m x 1 numerical array, the serial number of monosac. in the bone
% structure that has PM attached.
% mspos: n x 2 numerical array, the monosac. position.
% bonebondmap: n x n numerical array, the linakge map of the bone
% structure.
% bondmap: n x n numerical array, the linkage map of the whole structure.
% pmgroup: m x 2 cell array, 1st column contains serial number of monosac.
% in the glycan that has PM attached to it, 2nd column contains k x 1 cell
% array of numbers, k equals to the number of arms directly connected
% to the bone monosac, each element contains the serial number of monosac.
% in the arm.
% pmposstorage: m x 2 cell array of k x 2 numerical array, the position
% info of PM subtrees.
% pmindstorage: m x 2 cell array of k x 1 numerical array, the serial
% number of monosac in PM subtrees.
% allpmchild: m x 1 numerical array, the serial number of all monosac. that
% belong to a PM subtree.
%
% Output:
% mspos: modified 'mspos' after adding PM subtrees.
% safezoneind: n x 1 numerical array, indicator of adjusted monosac.
% PMroot: n x 1 numerical array, the serial number of monosac. in the bone
% structure that has PM attached to it.
%
% Note:
% N/A
%
% Example:
% N/A. Set breakpoints in main program.
%
% Children function:
% GLYTREETRACKER, CALCSAFEREGION, ANTICOLLISION
%

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2019, Research Foundation for State University of New York. All rights reserved
%


%%  Start combining PMS and bone structure, arm by arm
% 1. retrieve fork pos from bonebondmap

forkpos = mspos(sum(bonebondmap,2) > 1,:);
forkseq = [];
[~,~,xind] = unique(forkpos(:,1));
for i = max(xind):-1:1
    [~,yind] = sort(forkpos(xind == i,2));
    localforkseq = find(xind == i);
    forkseq = [forkseq;localforkseq(yind)];
end
PMroot = cell2mat(pmgroup(:,1));
safezoneind = false(size(bondmap,1),1);
for i = 1:length(forkseq)  % one fork at a time
    %% inside a fork
    thisfork = allfork(forkseq(i));
    forkchild = glytreetracker(bondmap,thisfork,[],'down');
    thesebranches = find(bonebondmap(thisfork,:));
    [~,branchvertseq] = sort(mspos(thesebranches,2));
    frklocalsafeind = false(size(bondmap,1),1);  % recording of finished fork and anticollision operations
    brnhlocalsafeind = false(size(bondmap,1),1);  % recording of finished branch, also used by in-branch PM anticollision
    for j = 1:length(thesebranches)
        %% inside a branch
        currentbranch = glytreetracker(bondmap,thesebranches(branchvertseq(j)),[],'down');
        currentbranchuncombined = glytreetracker(bondmap,thesebranches(branchvertseq(j)),allfork,'down');
        % stop at a finished fork, so if this fork contains a PM it won't
        % trigger an adding process.
        if ~safezoneind(thesebranches(branchvertseq(j)))
            pmsrootinext = ismember(PMroot,currentbranchuncombined);  % PMS roots in extension
            if ~any(pmsrootinext)  % a clean branch
                safezone = calcsaferegion(mspos,find(brnhlocalsafeind));  % existing safezone
                currentzone = calcsaferegion(mspos,currentbranch);  % this branch's safezone
                yfix = anticollision(safezone,currentzone,'up');
                mspos(currentbranch,2) = mspos(currentbranch,2) + yfix;
                brnhlocalsafeind(currentbranch) = 1;
            else  % branch containing PM
                rootind = find(ismember(PMroot,currentbranchuncombined));
                for k = 1:length(rootind)
                    tempPMSstruct = pmposstorage(rootind(k),:);
                    if ~isempty(tempPMSstruct{1})  % if there is a upper PMS subtree
                        upperpmpos = tempPMSstruct{1} + repmat(mspos(PMroot(rootind(k)),:),size(tempPMSstruct{1},1),1);
                        safezone = calcsaferegion(mspos,find(brnhlocalsafeind));
                        currentzone = calcsaferegion(upperpmpos,1:size(upperpmpos,1));
                        yfix = anticollision(safezone,currentzone,'down');
                        upperpmpos(:,2) = upperpmpos(:,2) + yfix;
                        brnhlocalsafeind(pmindstorage{rootind(k),1}) = 1;
                        mspos(pmindstorage{rootind(k),1},:) = upperpmpos;
                    end
                    if ~isempty(tempPMSstruct{2})  % if there is a lower PMS subtree
                        lowerpmpos = tempPMSstruct{2} + repmat(mspos(PMroot(rootind(k)),:),size(tempPMSstruct{2},1),1);
                        safezone = calcsaferegion(mspos,find(brnhlocalsafeind));
                        currentzone = calcsaferegion(lowerpmpos,1:size(lowerpmpos,1));
                        yfix = anticollision(safezone,currentzone,'up');
                        lowerpmpos(:,2) = lowerpmpos(:,2) + yfix;
                        brnhlocalsafeind(pmindstorage{rootind(k),2}) = 1;
                        mspos(pmindstorage{rootind(k),2},:) = lowerpmpos;
                    end
                end
            end
        else
            brnhlocalsafeind(currentbranch) = 1;
        end
        %% Now a branch is finished, consider its relationship with other branches in this fork
        if any(frklocalsafeind)
            dist = mspos(thesebranches(branchvertseq(j-1)),2) - mspos(thesebranches(branchvertseq(j)),2) + 1;
            safezone = calcsaferegion(mspos,find(frklocalsafeind));
            currentzone = calcsaferegion(mspos,currentbranch);
            currentzone = currentzone + [ones(size(currentzone,1),1)*dist,zeros(size(currentzone,1),2)];
            yfix = anticollision(safezone,currentzone,'up');
            mspos(currentbranch,2) = mspos(currentbranch,2)+ dist + yfix;
        end
        frklocalsafeind(currentbranch) = 1;  % update the safe zone index of the fork
        brnhlocalsafeind = false(size(bondmap,1),1);
    end
    safezoneind(frklocalsafeind) = 1;
    safezoneind(thisfork) = 1;
    % now working on fork point
    if mod(length(thesebranches),2)
        mspos(thisfork,2) = median(mspos(thesebranches,2));
    else
        mspos(thisfork,2) = mean(mspos(thesebranches,2));
    end
    rootpmind = thisfork;
    if ismember(thisfork,PMroot)  % if fork has PMS attached
        tempPMSstruct = pmposstorage(ismember(PMroot,thisfork),:);
        if ~isempty(tempPMSstruct{1})  % if there is a upper PMS subtree
            upperpmpos = tempPMSstruct{1}+ repmat(mspos(thisfork,:),size(tempPMSstruct{1},1),1);
            safezone = calcsaferegion(mspos,find(frklocalsafeind));
            currentzone = calcsaferegion(upperpmpos,1:size(upperpmpos,1));
            yfix = anticollision(safezone,currentzone,'down');
            upperpmpos(:,2) = upperpmpos(:,2) + yfix;
            mspos(pmindstorage{ismember(PMroot,thisfork),1},:) = upperpmpos;
            rootpmind = [rootpmind;pmindstorage{ismember(PMroot,thisfork),1}];
        end
        if ~isempty(tempPMSstruct{2})  % if there is a lower PMS subtree
            lowerpmpos = tempPMSstruct{2}+ repmat(mspos(thisfork,:),size(tempPMSstruct{2},1),1);
            safezone = calcsaferegion(mspos,find(frklocalsafeind));
            currentzone = calcsaferegion(lowerpmpos,1:size(lowerpmpos,1));
            yfix = anticollision(safezone,currentzone,'up');
            lowerpmpos(:,2) = lowerpmpos(:,2) + yfix;
            mspos(pmindstorage{ismember(PMroot,thisfork),2},:) = lowerpmpos;
            rootpmind = [rootpmind;pmindstorage{ismember(PMroot,thisfork),2}];
        end
    end
    % Fix branch (horizontal)
    currentzone = calcsaferegion(mspos,find(frklocalsafeind));
    safezone = calcsaferegion(mspos,rootpmind);
    xfix = anticollision(safezone,currentzone,'right');
    mspos(frklocalsafeind,1) = mspos(frklocalsafeind,1) + xfix;
    %% align with parent monosac
    parentms = find(bondmap(:,thisfork));
    if ~ismember(parentms,allfork)
        yfix = mspos(parentms,2) - mspos(thisfork,2);
        mspos(forkchild,2) = mspos(forkchild,2) + yfix;
    end
end
%% after adjusting everything after fork point, start adjusting those before
stub = setdiff(glytreetracker(bondmap,thisfork,[],'up'),thisfork);
for i = 1:length(stub)
    mspos(stub(i),2) = mspos(thisfork,2);
end
mspos = mspos - repmat(mspos(1,:),size(mspos,1),1);  % normalize all coordinates
end