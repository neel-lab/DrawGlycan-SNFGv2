function fixval = anticollision(exsafezone,currentzone,direction)
% ANTICOLLISION: generate fix values by performing collision test
%
% Syntax:
% fixval = anticollision(exsafezone,currentzone,direction)
%
% Input:
% exsafezone: n x 3 numerical array, the outline of the existing area. 1st
% column contains y values, 2nd column contains the x values of the left
% boundary, 3rd column contains the x values of the right boundary. Y
% values increases from up to down.
% currentzone: n x 3 numerical array, the outline of the area to be
% adjusted. Structure is same as 'exsafezone'.
% direction: 'up', 'down', 'left' or 'right', the direction of the dodging
% behavior.
%
% Output:
% fixval: number, the adjustment value.
%
% Note:
% The output fix value has direction '+/-', it just needs to be added to
% 'currentzone'.
% The description of the two outlines must be
% 'dense' enough, i.e. elements in the 1st column must have difference less
% than or equal to 1, otherwise program cannot detect collision.
%
% Example:
% fixval = anticollision([1 0 2;1.5 0,2],[0,1,3;1,1,3],'up')
%
% fixval =
%
%                        2.5
%
% Children function:
% COLLISIONEXAM
%

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2019, Research Foundation for State University of New York. All rights reserved
%


[iscleared,vertoverlay,horioverlay] = collisionexam(exsafezone,currentzone);
fixedcurrentzone = currentzone;
if ~iscleared
    if strcmpi(direction,'up')
        while ~iscleared
            fix = vertoverlay;
            fixedcurrentzone(:,1) = fixedcurrentzone(:,1) + fix;
            [iscleared,vertoverlay,~] = collisionexam(exsafezone,fixedcurrentzone);
        end
        fixval = fixedcurrentzone(1,1) - currentzone(1,1);
    elseif strcmpi(direction,'down')
        while ~iscleared
            fix = -vertoverlay;
            fixedcurrentzone(:,1) = fixedcurrentzone(:,1) + fix;
            [iscleared,vertoverlay,~] = collisionexam(exsafezone,fixedcurrentzone);
        end
        fixval = fixedcurrentzone(1,1) - currentzone(1,1);
    elseif strcmpi(direction,'right')
        while ~iscleared
            fix = horioverlay;
            fixedcurrentzone(:,2) = fixedcurrentzone(:,2) + fix;
            fixedcurrentzone(:,3) = fixedcurrentzone(:,3) + fix;
            [iscleared,~,horioverlay] = collisionexam(exsafezone,fixedcurrentzone);
        end
        fixval = fixedcurrentzone(1,2) - currentzone(1,2);
    elseif strcmpi(direction,'left')
        while ~iscleared
            fix = -horioverlay;
            fixedcurrentzone(:,2) = fixedcurrentzone(:,2) + fix;
            fixedcurrentzone(:,3) = fixedcurrentzone(:,3) + fix;
            [iscleared,~,horioverlay] = collisionexam(exsafezone,fixedcurrentzone);
        end
        fixval = fixedcurrentzone(1,2) - currentzone(1,2);
    end
else
    fixval = 0;
end
end


function [iscleared,vertoverlay,horioverlay] = collisionexam(exsafezone,currentzone)
iscleared = 1;
vertcol = ones(size(currentzone,1),size(exsafezone,1));
vertoverlay = zeros(size(currentzone,1),size(exsafezone,1));
horicol = vertcol;
horioverlay = vertoverlay;

for i = 1:size(currentzone,1)
    currentlayer = currentzone(i,:);
    for j = 1:size(exsafezone,1)
        if currentlayer(1) >= exsafezone(j,1) + 1 || currentlayer(1) <= exsafezone(j,1) - 1
            vertcol(i,j) = 0;
        else
            vertoverlay(i,j) = 1 + exsafezone(j,1) - currentlayer(1);
        end
        if currentlayer(3) <= exsafezone(j,2) - 1 || currentlayer(2) >= exsafezone(j,3) + 1
            horicol(i,j) = 0;
        else
            horioverlay(i,j) = 1 + min(exsafezone(j,3) - currentlayer(2),currentlayer(3)-exsafezone(j,2));
        end
    end
end
actualcollision = vertcol.*horicol;
if any(actualcollision(:))
    iscleared = 0;
end
vertoverlay = max(vertoverlay(logical(actualcollision)));
horioverlay = max(horioverlay(logical(actualcollision)));
end