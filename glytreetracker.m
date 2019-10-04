function children = glytreetracker(bondmap,start,stop,direction)
% GLYTREETRACKER: tracks up/down the glycan tree structure and return the
% serial number of the parent/children of specified monosac.
%
% Syntax:
% children = glytreetracker(bondmap,start,stop,'down')
% parent = glytreetracker(bondmap,start,stop,'up')
%
% Input:
% bondmap: n x n numerical array, the linkage map of the glycan.
% start: intiger, the serial number where the tracking starts.
% stop: m x 1 numerical array, the serial number(s) of the monosac.(s)
% where the tracking process should stop at.
% direction: 'up' means track the parent of the specified monosac. until
% the initial monosac. or the ones specified in input 'stop'. 'down' mans
% track the children of the specified monosac. until the last monosac. in
% the branch and stop at any monosac. included in 'stop'.
%
% Output:
% children/parent: the serial number of the children/parent monosac. of the
% specified one.
%
% Note:
% Monosacs. that exist in 'stop' will not be included in output.
%
% Example:
% N/A. Set breakpoints in parent program.
%
% Children function:
% N/A


%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2019, Research Foundation for State University of New York. All rights reserved
%


children = zeros(size(bondmap,1),1);
tobetracked = zeros(size(bondmap,1),1);
tobetracked(start) = 1;
if strcmpi(direction,'down')
    if ~isempty(stop)  % in order to improve speed,
        % the situation where 'stop' needs to be considered is written separately.
        while any(tobetracked)
            where = find(tobetracked,1,'first');
            if any(bondmap(where,:))
                tobetracked((bondmap(where,:) ~= 0)) = 1;
            end
            tobetracked(where) = 0;
            children(where) = 1;
            if any(ismember(stop,find(tobetracked)))
                tobetracked(stop) = 0;
            end
        end
    else
        while any(tobetracked)
            where = find(tobetracked,1,'first');
            if any(bondmap(where,:))
                tobetracked((bondmap(where,:) ~= 0)) = 1;
            end
            tobetracked(where) = 0;
            children(where) = 1;
        end
    end
elseif strcmpi(direction,'up')
    if ~isempty(stop)
        while any(tobetracked)
            where = find(tobetracked,1,'first');
            if any(bondmap(:,where))
                tobetracked((bondmap(:,where) ~= 0)) = 1;
            end
            tobetracked(where) = 0;
            children(logical(tobetracked)) = 1;
            if any(ismember(stop,find(tobetracked)))
                tobetracked(stop) = 0;
            end
        end
    else
        while any(tobetracked)
            where = find(tobetracked,1,'first');
            if any(bondmap(:,where))
                tobetracked((bondmap(:,where) ~= 0)) = 1;
            end
            tobetracked(where) = 0;
            children(logical(tobetracked)) = 1;
        end
    end
end
children = find(children);
end