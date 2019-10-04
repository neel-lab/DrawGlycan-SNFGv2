function estimatefigsize(plotinfoout,peppos,tcp,options)
% ESTIMATEFIGSIZE: Designate figure size according to the size of glycans
% and peptide backbone.
%
% Syntax:
% estimatefigsize(plotinfoout,peppos,options)
%
% Input:
% plotinfoout: m x n cell array, contains the plot information used in
% painting program, 1st column contains the positioning info for each
% monosaccharide in each glycan.
% peppos: 1 x n numerical array, the position of each amino acid in the
% peptide backbone. The numbers are the x position of amino acid. They will
% be put at [x,-1] with "verticalalignment" set to "top".
% options: struct, options.
%
%
% Output:
% horirec: coordinates for the left and right boundaries of the figure
% maxvert: figure height
% Figure size will be adjusted based on the estimation program
%
% Note:
% Figure size estimation is based on 70 pixels per unit. Program will
% read proposed glycan structure size and peptide backbone
% length, then reset figure size.
%
% Example:
% N/A
%
% Children function:
% N/A
%

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2019, Research Foundation for State University of New York. All rights reserved
%


pointsperunit = options.pointsperunit;
currentfig = options.figurehandle;
switch lower(options.workingmode)
    case 'g'
        mainmspos = plotinfoout.mspos;
        if isfield(plotinfoout,'CBITEM')
            cbitemmspos = plotinfoout.CBITEM.cbplotinfosto;
        else
            cbitemmspos = [];
        end
        cbxlim = [min(mainmspos(:,1)),max(mainmspos(:,1))];  % initial guess
        cbylim = [min(mainmspos(:,2)),max(mainmspos(:,2))];
        switch lower(options.orientation)
            case 'up'
                if ~isempty(cbitemmspos)
                    structwidth = 0;
                    structheight = 0;
                    for i = 1:length(cbitemmspos)
                        structwidth = structwidth + max(cbitemmspos(i).mspos(:,1)) - min(cbitemmspos(i).mspos(:,1)) + 1;
                        structheight = max([structheight,max(cbitemmspos(i).mspos(:,2)-cbitemmspos(i).mspos(1,2))]);
                    end
                    cbxlim = [min(cbitemmspos(1).mspos(:,1)),min(cbitemmspos(1).mspos(:,1)) + structwidth - 1];
                    cbylim(2) = cbylim(2) + structheight + 1;
                    cbxlim = [min([cbxlim(1),cbxlim(1)]),max([cbxlim(2),cbxlim(2)])];
                end
            case 'down'
                if ~isempty(cbitemmspos)
                    structwidth = 0;
                    structheight = 0;
                    for i = 1:length(cbitemmspos)
                        structwidth = structwidth + max(cbitemmspos(i).mspos(:,1)) - min(cbitemmspos(i).mspos(:,1)) + 1;
                        structheight = max([structheight,-min(cbitemmspos(i).mspos(:,2)-cbitemmspos(i).mspos(1,2))]);
                    end
                    cbxlim = [min(cbitemmspos(1).mspos(:,1)),min(cbitemmspos(1).mspos(:,1)) + structwidth - 1];
                    cbylim(1) = cbylim(1) - structheight - 1;
                    cbxlim = [min([cbxlim(1),cbxlim(1)]),max([cbxlim(2),cbxlim(2)])];
                end
            case 'left'
                if ~isempty(cbitemmspos)
                    structwidth = 0;
                    structheight = 0;
                    for i = 1:length(cbitemmspos)
                        structwidth = max([structwidth,-min(cbitemmspos(i).mspos(:,1)-cbitemmspos(i).mspos(1,1))]);
                        structheight = structheight + max(cbitemmspos(i).mspos(:,2)) - min(cbitemmspos(i).mspos(:,2)) + 1;
                    end
                    cbylim = [min(cbitemmspos(1).mspos(:,2)),min(cbitemmspos(1).mspos(:,2)) + structheight - 1];
                    cbxlim(1) = cbxlim(1) - structwidth - 1;
                    cbylim = [min([cbylim(1),cbylim(1)]),max([cbylim(2),cbylim(2)])];
                end
            case 'right'
                if ~isempty(cbitemmspos)
                    structwidth = 0;
                    structheight = 0;
                    for i = 1:length(cbitemmspos)
                        structwidth = max([structwidth,max(cbitemmspos(i).mspos(:,1)-cbitemmspos(i).mspos(1,1))]);
                        structheight = structheight + max(cbitemmspos(i).mspos(:,2)) - min(cbitemmspos(i).mspos(:,2)) + 1;
                    end
                    cbylim = [min(cbitemmspos(1).mspos(:,2)),min(cbitemmspos(1).mspos(:,2)) + structheight - 1];
                    cbxlim(2) = cbxlim(2) + structwidth + 1;
                    cbylim = [min([cbylim(1),cbylim(1)]),max([cbylim(2),cbylim(2)])];
                end
        end
        xlim = [min(cbxlim(1),tcp(1,1)),max(cbxlim(2),tcp(2,1))];
        ylim = [min(cbylim(1),tcp(1,2)),max(cbylim(2),tcp(2,2))];
        set(currentfig,'XLim',xlim + [-1,1],'YLim',ylim + [-1,1]);
        if ~options.isaxes
            set(currentfig,'Units','points')
            reshapedaxwidth = (diff(xlim) + 2)*pointsperunit;
            reshapedaxheight = (diff(ylim) + 2)*pointsperunit;
            set(currentfig,'position',[currentfig.Position(1:2),reshapedaxwidth,reshapedaxheight]);
            if options.resize
                parentfig = get(currentfig,'Parent');
                set(parentfig,'Units','points')
                pause(0.01);
                parentfigpos = get(parentfig,'Position');
                set(parentfig,'Position',[0,0,reshapedaxwidth*1.11,reshapedaxheight*1.11]);
                set(parentfig,'Units','normalized')
                set(currentfig,'Units','normalized','Position',[.05,.05,.9,.9]);
            end
        end
        
    case 'gp'
        LR = [min(peppos),max(peppos)];
        maxvert = 0;
        for i = 1:length(plotinfoout)
            thismspos = plotinfoout{i}.mspos;
            maxvert = max([maxvert,max(thismspos(:,2))]);
            thismsLR = [min(thismspos(:,1)),max(thismspos(:,1))];
            fldnames = lower(fieldnames(plotinfoout{i}));
            if ismember('adducttextpos',fldnames)
                adduct = plotinfoout{i}.ADDUCT;
                adducttext = '';
                for j = 1:size(adduct,1)
                    adducttext = [adducttext,adduct{j,1},' '];
                end
                adductLR = [min(plotinfoout{i}.bracketspos(:,1)),plotinfoout{1}.adducttextpos(1) + (length(adducttext)-1) * DrawGlycanPara.letterwidthpxl/DrawGlycanPara.pixelsperunit];
            else
                adductLR = [];
            end
            if ismember('cbitem',fldnames)
                cbLR = [0,0];
                cbplotinfos = plotinfoout{i}.CBITEM.cbplotinfosto;
                for j = 1:size(cbplotinfos,1)
                    cbitempos = plotinfoout{i}.CBITEM.cbplotinfosto(j).mspos;
                    cbLR = [min([cbLR(1),min(cbitempos(:,1))]),max([cbLR(2),max(cbitempos(:,1))])];
                    maxvert = max([maxvert,max(cbitempos(:,2))]);
                end
            else
                cbLR = [];
            end
            tempallLR = [thismsLR;adductLR;cbLR];
            LR = [min([LR(1),min(tempallLR(:,1))]),max([LR(2),max(tempallLR(:,2))])];
        end
        if maxvert == 0
            maxvert = 1;  % compensation for glycans with max. dist. = 0, one monosac or one monosac with fucose attached.
        end
        endcomp = diff(peppos);  % compensate space for both ends of peptide backbone
        if ~isempty(endcomp)  % need to consider this only when num. of AA > 1
            if LR(1) > min(peppos) - endcomp(1)
                LR(1) = min(peppos) - endcomp(1);  % if peptide backbone is long
            end
            if LR(2) < max(peppos) + endcomp(end)
                LR(2) = max(peppos) + endcomp(end);
            end
        end
        set(currentfig,'XLim',LR+[-1,1])
        set(currentfig,'YLim',[-2.5,maxvert+1])
        if ~options.isaxes
            set(currentfig,'position',[currentfig.Position(1:2)/4,(sum(abs(LR))+2)*pixelsperunit,(maxvert+3)*pixelsperunit])
        end
end
end
