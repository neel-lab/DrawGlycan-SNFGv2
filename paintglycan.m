function [plotinfoout,tcp] = paintglycan(plotinfoout,options,tcp)
% PAINTGLYCAN: Draws monosaccharides using symbols, glycosidic bonds using lines,
% and also depicts glycan fragmentations
%
% Syntax:
% tcp = paintglycan(plotinfoout,options,tcp)
%
% Input:
% plotinfoout: the info of glycan structure to be drawn.
% options: structure, global options.
% tcp: 2 x 2 numerical array, coordinates of the lower left and upper right
% corner of the figure. Can be left empty.
%
% Output:
% tcp: updated input "tcp" after drawing the glycan structure.
%
% Note:
% N/A.
%
% Example:
% N/A. Set breakpoints in program.
%
% Children function:
% GLYVIS
%

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2019, Research Foundation for State University of New York. All rights reserved
%


msposout = {plotinfoout.mspos};
bondmapout = {plotinfoout.bondmap};
alllinkout = {plotinfoout.alllinkout};
allmonosacout = {plotinfoout.allms};
directionseqout = {plotinfoout.directionseq};
identity = cell(size(plotinfoout));
for i = 1:length(plotinfoout)
    if isfield(plotinfoout(i),'CHAR') && ~isempty(plotinfoout(i).CHAR)
        thesechar = plotinfoout(i).CHAR;
        thesechar(:,1) = deal({1});
        identity{i} = [identity{i};thesechar];
    end
    if isfield(plotinfoout(i),'ANOMER') && ~isempty(plotinfoout(i).ANOMER)
        theseanomer = plotinfoout(i).ANOMER;
        theseanomer(:,2) = deal({1});
        theseanomer(:,1) = deal({2});
        identity{i} = [identity{i};theseanomer];
    end
end
if isfield(plotinfoout,'bonddescription')
    bondescription = {plotinfoout.bonddescription};
else
    bondescription = cell(size(plotinfoout));
end

if isempty(tcp)
    tcp = zeros(2,2);
end


if strcmpi(options.workingmode,'g')
    switch lower(options.orientation)
        case 'up'
            previousmspos = -[min(msposout{1}(:,1)) + options.structspacing,0];
        case 'down'
            previousmspos = -[min(msposout{1}(:,1)) + options.structspacing,0];
        case 'left'
            previousmspos = -[0,min(msposout{1}(:,2)) + options.structspacing];
        case 'right'
            previousmspos = -[0,min(msposout{1}(:,2)) + options.structspacing];
    end
    availableMSlist = properties(MonoColor);
    for i = 1:numel(msposout)
        mscolor = cell(size(allmonosacout{i}));
        msshape = mscolor;
        %% special case: deoxyhexose 6dAlt & 6dTal
        for j = 1:numel(allmonosacout{i})
            if ismember(allmonosacout{i}{j},availableMSlist)
                mscolor{j} = MonoColor.(allmonosacout{i}{j});
                msshape{j} = MonoShape.(allmonosacout{i}{j});
            elseif strcmp(allmonosacout{i}{j},'6dTal')
                mscolor{j} = 'Light Blue';
                msshape{j} = 'Filled Triangle';
            elseif strcmp(allmonosacout{i}{j},'6dAlt')
                mscolor{j} = 'Pink';
                msshape{j} = 'Filled Triangle';
            elseif strcmp(allmonosacout{i}{j},'6dGul')
                mscolor{j} = 'Orange';
                msshape{j} = 'Filled Triangle';
            elseif strcmp(allmonosacout{i}{j},'6dAltNAc')
                mscolor{j} = 'Pink';
                msshape{j} = 'Divided Triangle';
            elseif strcmp(allmonosacout{i}{j},'6dTalNAc')
                mscolor{j} = 'Light Blue';
                msshape{j} = 'Divided Triangle';
            elseif strcmp(allmonosacout{i}{j},'4eLeg')
                mscolor{j} = 'Light Blue';
                msshape{j} = 'Flat Diamond';
            else
                mscolor{j} = 'White';
                msshape{j} = 'Flat Hexagon';
            end
        end
        mspara.shape = msshape;
        mspara.color = mscolor;
        mspara.size = ones(size(mscolor)) * options.monosacsize;
        mspara.periwidth = ones(size(mscolor)) * options.msperiwidth;
        mspara.name = allmonosacout{i};
        switch lower(options.orientation)
            case 'up'
                previousmspos = previousmspos + [min(msposout{i}(:,1)) + options.structspacing,0];
            case 'down'
                previousmspos = previousmspos + [min(msposout{i}(:,1)) + options.structspacing,0];
            case 'left'
                previousmspos = previousmspos + [0,min(msposout{i}(:,2)) + options.structspacing];
            case 'right'
                previousmspos = previousmspos + [0,min(msposout{i}(:,2)) + options.structspacing];
        end
        mspos = msposout{i} + repmat(previousmspos,size(msposout{i},1),1);
        glyvis(mspos,bondmapout{i},alllinkout{i},bondescription{i},mspara,directionseqout{i},identity{i},options)
        switch lower(options.orientation)
            case 'up'
                previousmspos = [max(mspos(:,1)),0];
            case 'down'
                previousmspos = [max(mspos(:,1)),0];
            case 'left'
                previousmspos = [0,max(mspos(:,2))];
            case 'right'
                previousmspos = [0,max(mspos(:,2))];
        end
        tcp = [min([mspos(:,1);tcp(1,1)]),min([mspos(:,2);tcp(1,2)]);...
            max([mspos(:,1);tcp(2,1)]),max([mspos(:,2);tcp(2,2)])];  % update 2 corner pos.
        ax = options.figurehandle;
        xlim = [tcp(1,1),tcp(2,1)];
        set(ax,'Units','points');
        axwidth = get(ax,'Position');
        axwidth = axwidth(3);
        fontsizeconst = options.pointsperunit * DrawGlycanPara.fontsizeconst;
        anomerhorifix = 0;
        if isfield(plotinfoout(i),'ANOMER') && ~isempty(plotinfoout(i).ANOMER)
            switch lower(options.orientation)
                case {'up','down'}
                    % INTENTIONALLY LEFT BLANK
                case 'left'
                    %                     letterwidth = diff(xlim)/axwidth;
                    letterwidth = .05;
                    anomerhorifix = letterwidth * options.fontsize * 1.5*...
                        fontsizeconst * length(plotinfoout(i).ANOMER{1,1});
                    tcp(2,1) = tcp(2,1) + anomerhorifix;
                case 'right'
                    %                     letterwidth = diff(xlim)/axwidth;
                    letterwidth = .05;
                    anomerhorifix = letterwidth * options.fontsize *1.5 *...
                        fontsizeconst * length(plotinfoout(i).ANOMER{1,1});
                    tcp(1,1) = tcp(1,1) - anomerhorifix;
            end
        end
        if isfield(plotinfoout(i),'ADDUCT') && isfield(plotinfoout(i),'bracketspos') &&...
                isfield(plotinfoout(i),'adducttextpos')
            switch lower(options.orientation)
                case 'left'
                    plotinfoout(i).bracketspos(5:8,1) = plotinfoout(i).bracketspos(5:8,1) + anomerhorifix;
                    plotinfoout(i).adducttextpos(1) = plotinfoout(i).adducttextpos(1) + anomerhorifix;
                case 'right'
                    plotinfoout(i).bracketspos(1:4,1) = plotinfoout(i).bracketspos(1:4,1) - anomerhorifix;
                    plotinfoout(i).adducttextpos(1) = plotinfoout(i).adducttextpos(1) - anomerhorifix;
            end
        end
    end
    set(ax,'Units','normalized');
elseif strcmpi(options.workingmode,'gp')
    availableMSlist = properties(MonoColor);
    for i = 1:numel(msposout)
        mscolor = cell(size(allmonosacout{i}));
        msshape = mscolor;
        for j = 1:numel(allmonosacout{i})
            if ismember(allmonosacout{i}{j},availableMSlist)
                mscolor{j} = MonoColor.(allmonosacout{i}{j});
                msshape{j} = MonoShape.(allmonosacout{i}{j});
            elseif strcmp(allmonosacout{i}{j},'6dTal')
                mscolor{j} = 'Light Blue';
                msshape{j} = 'Filled Triangle';
            elseif strcmp(allmonosacout{i}{j},'6dAlt')
                mscolor{j} = 'Pink';
                msshape{j} = 'Filled Triangle';
            else
                mscolor{j} = 'White';
                msshape{j} = 'Flat Hexagon';
            end
        end
        mspara.shape = msshape;
        mspara.color = mscolor;
        mspara.size = ones(size(mscolor)) * options.monosacsize;
        mspara.periwidth = ones(size(mscolor)) * options.msperiwidth;
        mspara.name = allmonosacout{i};
        mspos = msposout{i};
        glyvis(mspos,bondmapout{i},alllinkout{i},bondescription{i},mspara,directionseqout{i},identity{i},options)
        tcp = [min([mspos(:,1);tcp(1,1)]),min([mspos(:,2);tcp(1,2)]);...
            max([mspos(:,1);tcp(2,1)]),max([mspos(:,2);tcp(2,2)])];  % update 2 corner pos.
    end
end
end