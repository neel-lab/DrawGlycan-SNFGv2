function paintlink(mspos,bondmap,bondescription,options)
% PAINTLINK: draw the glycosidic bonds.
%
% Syntax:
% paintlink(mspos,bondmap,bondescription,options)
%
% Input:
% mspos: n x 2 numerical array, the position of monosaccharides.
% bondmap: n x n numerical array, how each monosac. link to each other.
% bonddescription: n x 2 cell array, instructions for drawing specific
% bonds, 1st column tells which link, 2nd column tells how.
% options: global options for drawing.
%
% Output:
% N/A
%
% Note:
% Available line styles are: 'BOLD','ZIG','DASH','DASHB','DOUBLE','
% TRIPLE','WAVY','WEDGE'
%
% Example:
% N/A. Set breakpoints in main program.
%
% Children function:
% N/A

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2019, Research Foundation for State University of New York. All rights reserved
%


bondtyp = {'BOLD','ZIG','DASH','DASHB','DOUBLE','TRIPLE','WAVY','WEDGE'};
theta = pi/4;
ziglen = options.monosacsize/4;
thisms = 1;
bondinfopos = [];
bondinfotyp = {};
fig = options.figurehandle;
if ~isempty(bondescription)
    bondinfopos = cell2mat(bondescription(:,1));
    bondinfotyp = bondescription(:,2);
    bondinfopos = bondinfopos(ismember(bondinfotyp,bondtyp));
    bondinfotyp = bondinfotyp(ismember(bondinfotyp,bondtyp));
end
while thisms <= size(mspos,1)
    children = bondmap(thisms,:);
    if any(children)
        children = find(children);
        for i = 1:length(children)
            yn = ismember(bondinfopos,children(i));
            if any(yn)
                temp_loc = find(yn);
                for j = 1:length(temp_loc)
                    temp_opttyp = bondinfotyp{temp_loc(j)};
                    alpha = atan((mspos(children(i),2) - mspos(thisms,2))/(mspos(children(i),1) - mspos(thisms,1)));
                    if mspos(children(i),1) - mspos(thisms,1) < 0
                        alpha = alpha + pi;
                    elseif mspos(children(i),2) - mspos(thisms,2) == 0
                        if mspos(children(i),1) - mspos(thisms,1) < 0
                            alpha = alpha + pi;
                        end
                    end
                    eqn = @(x,b) tan(alpha - pi/2) * x + b;
                    switch temp_opttyp
                        case 'BOLD'
                            plot(fig,[mspos(thisms,1),mspos(children(i),1)],[mspos(thisms,2),mspos(children(i),2)],'k','linewidth',options.bondwidth*2)
                        case 'ZIG'
                            midpoint = [mspos(children(i),1) + mspos(thisms,1),mspos(children(i),2) + mspos(thisms,2)]/2;
                            waypointfixa = abs(ziglen*cos(theta-alpha));
                            waypointfixb = abs(ziglen*sin(theta-alpha));
                            if alpha >= 0 && alpha <= pi/4  % Ia
                                waypointfix = [-waypointfixa,waypointfixb];
                            elseif alpha > pi/4 && alpha <= pi*3/4  % Ib - IIa
                                waypointfix = [-waypointfixa,-waypointfixb];
                            elseif alpha > pi*3/4 && alpha <= pi*5/4  % IIb - IIIa
                                waypointfix = [waypointfixa,-waypointfixb];
                            elseif (alpha >= pi*5/4 && alpha <= pi*3/2) || (alpha >= -pi/2 && alpha < -pi/4)  % IIIb - IVa
                                waypointfix = [waypointfixa,waypointfixb];
                            elseif alpha >= -pi/4 && alpha <0  % IVb
                                waypointfix = [-waypointfixa,waypointfixb];
                            end
                            waypoint1 = midpoint + waypointfix;
                            waypoint2 = midpoint - waypointfix;
                            plot(fig,[mspos(thisms,1),waypoint1(1),waypoint2(1),mspos(children(i),1)],[mspos(thisms,2),waypoint1(2),waypoint2(2),mspos(children(i),2)],'k','linewidth',options.bondwidth)
                        case 'DASH'
                            plot(fig,[mspos(thisms,1),mspos(children(i),1)],[mspos(thisms,2),mspos(children(i),2)],'--k','linewidth',options.bondwidth)
                        case 'DASHB' % 10 lines
                            n=10;
                            waypoints = zeros(1,4);
                            centerpoint = [linspace(mspos(thisms,1),mspos(children(i),1),n);linspace(mspos(thisms,2),mspos(children(i),2),n)]';
                            if mod(alpha,pi) == 0
                                for k = 1:n
                                    waypoints = [centerpoint(k,1),centerpoint(k,1),centerpoint(k,2) + ziglen*2/3*k/n,centerpoint(k,2) - ziglen*2/3*k/n];
                                    plot(fig,waypoints(1:2),waypoints(3:4),'k','linewidth',options.bondwidth)
                                end
                            elseif mod(alpha,pi/2) == 0
                                for k = 1:n
                                    waypoints = [centerpoint(k,1) - ziglen*2/3*k/n,centerpoint(k,1) + ziglen*2/3*k/n,centerpoint(k,2),centerpoint(k,2)];
                                    plot(fig,waypoints(1:2),waypoints(3:4),'k','linewidth',options.bondwidth)
                                end
                            else
                                for k = 1:n
                                    b = centerpoint(k,2) - tan(alpha - pi/2) * centerpoint(k,1);
                                    waypoints(1) = centerpoint(k,1) - ziglen*2/3*k/n * cos(alpha - pi/2);
                                    waypoints(2) = centerpoint(k,1) + ziglen*2/3*k/n * cos(alpha - pi/2);
                                    waypoints(3) = eqn(waypoints(1),b);
                                    waypoints(4) = eqn(waypoints(2),b);
                                    plot(fig,waypoints(1:2),waypoints(3:4),'k','linewidth',options.bondwidth)
                                end
                            end
                        case 'DOUBLE' % 10 lines
                            %                             n=10;
                            %                             waypoints = zeros(1,4);
                            %                             centerpoint = [linspace(mspos(thisms,1),mspos(children(i),1),n);linspace(mspos(thisms,2),mspos(children(i),2),n)]';
                            startpoint = mspos(thisms,:);
                            endpoint = mspos(children(i),:);
                            if mod(alpha,pi) == 0
                                plot(fig,[startpoint(1),endpoint(1)],[startpoint(2),endpoint(2)] - ziglen/2,'k','linewidth',options.bondwidth)
                                plot(fig,[startpoint(1),endpoint(1)],[startpoint(2),endpoint(2)] + ziglen/2,'k','linewidth',options.bondwidth)
                                %                                 for k = 1:n
                                %                                     waypoints = [centerpoint(k,1),centerpoint(k,1),centerpoint(k,2) + ziglen/2,centerpoint(k,2) - ziglen/2];
                                %                                     plot(fig,waypoints(1:2),waypoints(3:4),'k','linewidth',options.bondwidth)
                                %                                 end
                            elseif mod(alpha,pi/2) == 0
                                plot(fig,[startpoint(1),endpoint(1)] - ziglen/2,[startpoint(2),endpoint(2)],'k','linewidth',options.bondwidth)
                                plot(fig,[startpoint(1),endpoint(1)] + ziglen/2,[startpoint(2),endpoint(2)],'k','linewidth',options.bondwidth)
                                %                                 for k = 1:n
                                %                                     waypoints = [centerpoint(k,1) - ziglen/2,centerpoint(k,1) + ziglen/2,centerpoint(k,2),centerpoint(k,2)];
                                %                                     plot(fig,waypoints(1:2),waypoints(3:4),'k','linewidth',options.bondwidth)
                                %                                 end
                            else
                                plot(fig,[startpoint(1),endpoint(1)] - ziglen/2 * cos(alpha - pi/2),[startpoint(2),endpoint(2)] - ziglen/2 * sin(alpha - pi/2),'k','linewidth',options.bondwidth)
                                plot(fig,[startpoint(1),endpoint(1)] + ziglen/2 * cos(alpha - pi/2),[startpoint(2),endpoint(2)] + ziglen/2 * sin(alpha - pi/2),'k','linewidth',options.bondwidth)
                                %                                 for k = 1:n
                                %                                     b = centerpoint(k,2) - tan(alpha - pi/2) * centerpoint(k,1);
                                %                                     waypoints(1) = centerpoint(k,1) - ziglen/2 * cos(alpha - pi/2);
                                %                                     waypoints(2) = centerpoint(k,1) + ziglen/2 * cos(alpha - pi/2);
                                %                                     waypoints(3) = eqn(waypoints(1),b);
                                %                                     waypoints(4) = eqn(waypoints(2),b);
                                %                                     plot(fig,waypoints(1:2),waypoints(3:4),'k','linewidth',options.bondwidth)
                                %                                 end
                            end
                        case 'TRIPLE'
                            startpoint = mspos(thisms,:);
                            endpoint = mspos(children(i),:);
                            if mod(alpha,pi) == 0
                                plot(fig,[startpoint(1),endpoint(1)],[startpoint(2),endpoint(2)],'k','linewidth',options.bondwidth)
                                plot(fig,[startpoint(1),endpoint(1)],[startpoint(2),endpoint(2)] - ziglen*2/3,'k','linewidth',options.bondwidth)
                                plot(fig,[startpoint(1),endpoint(1)],[startpoint(2),endpoint(2)] + ziglen*2/3,'k','linewidth',options.bondwidth)
                            elseif mod(alpha,pi/2) == 0
                                plot(fig,[startpoint(1),endpoint(1)],[startpoint(2),endpoint(2)],'k','linewidth',options.bondwidth)
                                plot(fig,[startpoint(1),endpoint(1)] - ziglen*2/3,[startpoint(2),endpoint(2)],'k','linewidth',options.bondwidth)
                                plot(fig,[startpoint(1),endpoint(1)] + ziglen*2/3,[startpoint(2),endpoint(2)],'k','linewidth',options.bondwidth)
                            else
                                plot(fig,[startpoint(1),endpoint(1)] - ziglen*2/3 * cos(alpha - pi/2),[startpoint(2),endpoint(2)],'k','linewidth',options.bondwidth)
                                plot(fig,[startpoint(1),endpoint(1)] - ziglen*2/3 * cos(alpha - pi/2),[startpoint(2),endpoint(2)] - ziglen*2/3 * sin(alpha - pi/2),'k','linewidth',options.bondwidth)
                                plot(fig,[startpoint(1),endpoint(1)] + ziglen*2/3 * cos(alpha - pi/2),[startpoint(2),endpoint(2)] + ziglen*2/3 * sin(alpha - pi/2),'k','linewidth',options.bondwidth)
                            end
                        case 'WAVY'
                            n = 11;  % 11 half units
                            m = 5;  % 5 points per half unit
                            %% half-circle
                            centerpoint = [linspace(mspos(thisms,1),mspos(children(i),1),n*2+1);linspace(mspos(thisms,2),mspos(children(i),2),n*2+1)]';
                            radius = sqrt((centerpoint(2,2) - centerpoint(1,2))^2 + (centerpoint(2,1) - centerpoint(1,1))^2);
                            centerpoint = centerpoint(2:2:end,:);
                            for k = 1:n
                                tempwaypoints = zeros(5,2);
                                if mod(k,2)
                                    tempangles = linspace(alpha-pi,alpha,m);
                                else
                                    tempangles = linspace(alpha,alpha+pi,m);
                                end
                                for l = 1:m
                                    tempwaypoints(l,1) = radius * cos(tempangles(l)) + centerpoint(k,1);
                                    tempwaypoints(l,2) = radius * sin(tempangles(l)) + centerpoint(k,2);
                                end
                                plot(fig,tempwaypoints(:,1),tempwaypoints(:,2),'k','linewidth',options.bondwidth)
                            end
                            
                            %% sine wave
                            %                             len = sqrt((mspos(thisms,1)-mspos(children(i),1))^2 + (mspos(thisms,2)-mspos(children(i),2))^2)/cos(alpha);
                            %                             x = linspace(mspos(thisms,1),mspos(thisms,1)+len,n*m);
                            %                             waypoints = [x;arrayfun(@(x) sin(pi/m*x)*ziglen,0:n*m-1)]';
                            %                             waypoints = waypoints * [cos(alpha-pi),sin(alpha-pi);-sin(alpha-pi),cos(alpha-pi)];
                            %                             waypoints(:,1) = waypoints(:,1) + mspos(thisms,1) - waypoints(1,1);
                            %                             waypoints(:,2) = waypoints(:,2) + mspos(thisms,2) - waypoints(1,2);
                            %                             plot(fig,waypoints(:,1),waypoints(:,2),'k','linewidth',options.bondwidth)
                        case 'WEDGE'
                            waypoints = zeros(3,2);
                            waypoints(1,1) = mspos(thisms,1) + options.monosacsize/2 * cos(alpha);
                            waypoints(1,2) = (waypoints(1,1) - mspos(thisms,1)) * tan(alpha) + mspos(thisms,2);
                            tempbotcenter(1,1) = mspos(children(i),1) - options.monosacsize/2 * cos(alpha);
                            if mod(alpha,pi/2) == 0
                                tempbotcenter(1,2) = mspos(children(i),2) - sin(alpha) * options.monosacsize/2;
                            else
                                tempbotcenter(1,2) = (tempbotcenter(1,1) - mspos(thisms,1)) * tan(alpha) + mspos(thisms,2);
                            end
                            b = tempbotcenter(1,2) - tan(alpha + pi/2) * tempbotcenter(1,1);
                            if mod(alpha,pi) == 0
                                waypoints(2,:) = [tempbotcenter(1,1),tempbotcenter(1,2) - ziglen/2];
                                waypoints(3,:) = [tempbotcenter(1,1),tempbotcenter(1,2) + ziglen/2];
                            elseif mod(alpha,pi/2) == 0
                                waypoints(2,:) = [tempbotcenter(1,1) - ziglen/2,tempbotcenter(1,2)];
                                waypoints(3,:) = [tempbotcenter(1,1) + ziglen/2,tempbotcenter(1,2)];
                            else
                                waypoints(2,:) = [tempbotcenter(1,1) - ziglen/2*sin(alpha),eqn(tempbotcenter(1,1) - ziglen/2*sin(alpha),b)];
                                waypoints(3,:) = [tempbotcenter(1,1) + ziglen/2*sin(alpha),eqn(tempbotcenter(1,1) + ziglen/2*sin(alpha),b)];
                            end
                            patch(fig,waypoints(:,1),waypoints(:,2),'k','EdgeColor','none')
                    end
                end
            else
                plot(fig,[mspos(thisms,1),mspos(children(i),1)],[mspos(thisms,2),mspos(children(i),2)],'k','linewidth',options.bondwidth)
            end
        end
    end
    thisms = thisms + 1;
end
end