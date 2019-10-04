function glyvis(mspos,bondmap,linkage,bondescription,mspara,directionseq,identity,options)
% GLYVIS: Glycan visualizer is the main function for displaying glycans
%
% Syntax:
% glyvis(mspos,bondmap,linkage,bondbreak,mspara,directionseq,identity,options)
%
% Input:
% mspos: n x 2 numerical array, the monosaccharide position in glycan
% bondmap: m x m numerical array, the map showing connection between
% monosaccharides.
% linkage: n x 1 cell array of strings, the glycosidic bond information.
% bondbreak: 1 x 3 cell array, 1st element contains marker position, 2nd
% position the type, 3rd the content.
% mspara: structure, the shape, color size, perimeter width and monosac.
% name of the monosac. to be drawn.
% directionseq: 1 x n numerical array, the orientation of monosac. symbols
% identity: the nature of each symbol, can be SHAPE, ANOMER or CHAR
% options: global drawing options
%
% Output:
% N/A
%
% Note:
% N/A
%
% Example:
% N/A. Set breakpoints in main program.
%
% Children function:
% DRAWSHAPES
%

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2019, Research Foundation for State University of New York. All rights reserved
%


identityseq = cell(size(mspos,1),1);
identityseq(:) = {'SHAPE'};
fig = options.figurehandle;

if ~isempty(identity)
    identityinfo = cell2mat(identity);
    identityseq(identityinfo(identityinfo(:,1) == 1,2)) = deal({'CHAR'});
    identityseq(identityinfo(identityinfo(:,1) == 2,2)) = deal({'ANOMER'});
end
paintlink(mspos,bondmap,bondescription,options)
if ~isempty(bondescription)
    bondbreakpos = cell2mat(bondescription(:,1));
    bondbreaktyp = bondescription(:,2);
    bondbreakinfo = bondescription(:,3);
else
    bondbreakpos = [];
    bondbreaktyp = '';
    bondbreakinfo = '';
end
if ischar(bondbreakinfo)
    bondbreakinfo = {bondbreakinfo};
end
if ismember(lower(options.orientation),{'right','down'})
    alignment = {'left','right'};
elseif ismember(lower(options.orientation),{'left','up'})
    alignment = {'right','left'};
end
if strcmpi(options.orientation,'left') || strcmpi(options.orientation,'right')
    msaddtextrot = 0;
elseif strcmpi(options.orientation,'up') || strcmpi(options.orientation,'down')
    msaddtextrot = 90;
end
for i = 1:size(mspos,1)
    linkdist = mspara.size(i)/2*(1+options.linkinfodist);  % distance: outer rim + 20% radius (default)
    drawshapes(mspara.shape{i},mspos(i,:),mspara.size(i),mspara.color{i},...
        mspara.periwidth(i),directionseq(i),mspara.name{i},identityseq{i},options);
    if strcmpi(mspara.shape{i},'flat hexagon') && strcmpi(mspara.color{i},'White') && ...
            strcmpi(identityseq{i},'SHAPE') && ~strcmpi(mspara.name{i},'unknown')
        text(fig,mspos(i,1),mspos(i,2),mspara.name{i}(1),'fontsize',options.fontsize*1.5,...
            'horizontalalignment','center');
    end
end
for i = 1:size(mspos,1)
    linktyp = linkage{i}(1);
    if strcmpi(linktyp,'a')
        linktyp = 945;
        %         linktyp = '\alpha';
    elseif strcmpi(linktyp,'b')
        linktyp = 946;
        %         linktyp = '\beta';
    elseif strcmpi(linktyp,'?')
        linktyp = '';
    end
    linknum = regexp(linkage{i}(2:end),'-[?\d/\\]+','match');
    linknum = linknum{1}(2:end);
    if strcmpi(linknum,'?')
        linknum = '';
    end
    
    if i == 1
        if strcmpi(options.orientation,'left')
            lastmonosac = mspos(1,:) + [1,0];
        elseif strcmpi(options.orientation,'right')
            lastmonosac = mspos(1,:) - [1,0];
        elseif strcmpi(options.orientation,'up')
            lastmonosac = mspos(1,:) - [0,1];
        elseif strcmpi(options.orientation,'down')
            lastmonosac = mspos(1,:) + [0,1];
        end
    else
        lastmonosac = mspos(find(bondmap(:,i),1,'first'),:);
    end
    if mspos(i,1)-lastmonosac(1) ~= 0
        alpha = atan((mspos(i,2)-lastmonosac(2))/(mspos(i,1)-lastmonosac(1)));
    else
        alpha = pi/2;
    end
    xdif = mspos(i,1)-lastmonosac(1);
    ydif = mspos(i,2)-lastmonosac(2);
    msdistance = sqrt(xdif^2+ydif^2);
    if xdif > 0 && ydif <= 0
        alpha = alpha + 2*pi;
    elseif xdif <= 0 && ydif < 0
        alpha = alpha + pi;
    elseif xdif < 0 && ydif >= 0
        alpha = alpha + pi;
    end  % fix alpha value
    typx = mspos(i,1) - linkdist * cos(alpha - options.linkinfotheta/180*pi);
    typy = mspos(i,2) - linkdist * sin(alpha - options.linkinfotheta/180*pi);
    numx = lastmonosac(1) + linkdist * cos(alpha + options.linkinfotheta/180*pi);
    numy = lastmonosac(2) + linkdist * sin(alpha + options.linkinfotheta/180*pi);
    if strcmpi(options.showlink,'yes')
        text(fig,typx,typy,char(double(linktyp)),'fontsize',options.linkfontsize,'horizontalalignment','center');
        text(fig,numx,numy,char(double(linknum)),'fontsize',options.linkfontsize,'horizontalalignment','center');
    end
    whichbondbreak = ismember(bondbreakpos,i);
    if any(whichbondbreak)
        breakpos = bondbreakpos(whichbondbreak);
        breaktype = bondbreaktyp(whichbondbreak);
        breakinfo = bondbreakinfo(whichbondbreak);
        centx = mean([lastmonosac(1),mspos(i,1)]);
        centy = mean([lastmonosac(2),mspos(i,2)]);
        linestart = [centx + options.bondbreaksiglength/1.4 * cos(alpha + pi/2),...
            centy + options.bondbreaksiglength/1.4 * sin(alpha + pi/2)];
        lineend = [centx - options.bondbreaksiglength/1.4 * cos(alpha + pi/2),...
            centy - options.bondbreaksiglength/1.4 * sin(alpha + pi/2)];
        for j = 1:numel(breakpos)
            typelinefix = [options.bondbreaksiglength/ 5 * cos( 0 /180*pi+alpha),...
                options.bondbreaksiglength/ 5 * sin( 0 /180*pi+alpha)];
            typetextfix = [options.bondbreaksiglength/ 5 * cos( 90 /180*pi+alpha),...
                options.bondbreaksiglength/ 5 * sin( 90 /180*pi+alpha)];
            
            if strcmpi(breaktype{j},'NR')
                typelineend = linestart + typelinefix;
                plot(fig,[typelineend(1),linestart(1)],[typelineend(2),linestart(2)],'k');
                text(fig,typelineend(1)+  typetextfix(1),typelineend(2)+ typetextfix(2),...
                    char(double(breakinfo{j})),'horizontalalignment',alignment{1},...
                    'fontsize',options.fontsize);
                plot(fig,[linestart(1),lineend(1)],[linestart(2),lineend(2)],'k');
            elseif strcmpi(breaktype{j},'R')
                typelineend = lineend - typelinefix;
                plot(fig,[typelineend(1),lineend(1)],[typelineend(2),lineend(2)],'k');
                text(fig,typelineend(1)-  typetextfix(1),typelineend(2)- typetextfix(2),...
                    char(double(breakinfo{j})),'horizontalalignment',alignment{2},...
                    'fontsize',options.fontsize);
                plot(fig,[linestart(1),lineend(1)],[linestart(2),lineend(2)],'k');
            elseif strcmpi(breaktype{j},'U')
                if msaddtextrot  % vertical glycan
                    text(fig,mspos(i,1) - options.monosacsize*1.02/2,mspos(i,2),...
                        char(double(breakinfo{j})),'horizontalalignment','center',...
                        'verticalalignment','middle','fontsize',options.fontsize * 1.25);
                else  % horizontal glycan
                    text(fig,mspos(i,1),mspos(i,2) + options.monosacsize*1.02/2,char(double(breakinfo{j})),...
                        'horizontalalignment','center',...
                        'verticalalignment','bottom','fontsize',options.fontsize * 1.25);
                end
            elseif strcmpi(breaktype{j},'D')
                if msaddtextrot  % vertical glycan
                    text(fig,mspos(i,1) + options.monosacsize*1.02/2,mspos(i,2),char(double(breakinfo{j})),...
                        'horizontalalignment','center',...
                        'verticalalignment','middle','fontsize',options.fontsize * 1.25);
                else  % horizontal glycan
                    text(fig,mspos(i,1),mspos(i,2) - options.monosacsize*1.02/2,char(double(breakinfo{j})),...
                        'horizontalalignment','center',...
                        'verticalalignment','top','fontsize',options.fontsize * 1.25);
                end
            elseif strcmpi(breaktype{j},'C')
                if ismember(breakinfo{j},{'f','p','o'})
                    fontangle = 'italic';
                else
                    fontangle = 'normal';
                end
                text(fig,mspos(i,1),mspos(i,2),char(double(breakinfo{j})),'horizontalalignment','center',...
                    'verticalalignment','middle',...
                    'fontsize',options.fontsize * 1.25,'FontAngle',fontangle);
            end
        end
    end
end
end