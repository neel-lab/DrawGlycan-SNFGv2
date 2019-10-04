function [cleanerseq,specialoptions] = getglydressing(GLYseq)
% GETGLYDRESSING: analyze glycan sequence (glycan + customization),
% return separated information
%
% Syntax:
% [cleanerseq,specialoptions] = getglydressing(GLYseq)
%
% Input:
% GLYseq: string, IUPAC sequence of glycan, sequence may contain
% customization information: anomer, adduct, curly bracket, etc.
%
% Output:
% cleanerseq: string, IUPAC sequence of glycan, deprived of customization.
% specialoptions: struct, customization information of glycan, values are
% 1 x 2 cell array, 1st element is value of the option, 2nd is the monosac.
% number it applies to.
%
% ------------------------ available fieldnames are ------------------------------
% CURLYBRACKET       structures appear outside curly bracket
% RS                             position where repeated units start
% RE                             position where repeated units end
% U                               text that appear above monosac. symbol
% D                               text that appear below monosac. symbol
% P                                monosac. that is perpendicular to main structure
% CHAR                         item to be shown in text form
% R                                fragment at the reducing end
% NR                              fragment at the non-reducing end
% BOLD                          bond shown as bold line
% ZIG                              bond shown as zigzag line
% DASH                          bond shown as dashed line
% ADDUCT                    glycan adduct
% ANOMER                   anomer
%--------------------------------------------------------------------------------------
%
% Note:
% This function handles glycan only.
% This function is specifically designed to handle curly bracket item. All
% other options are handled through GETSEQOPTIONS.
% Curly bracket items must be placed at the beginning of the main structure.
%
% Example:
% [cleanerseq,specialoptions] = getglydressing('{Neu5Ac,Neu5Ac}{Man[Man]Man[GlcNAc[GlcNAc]Man]Man(-adduct "Na+")GlcNAc()GlcNAc()UDP(-char)}')
%
% cleanerseq =
%
% Man[Man]Man[GlcNAc[GlcNAc]Man]Man()GlcNAc()GlcNAc()UDP()
%
%
% specialoptions =
%
%     CURLYBRACKET: {'Neu5Ac,Neu5Ac'  [1 2 3 4 5 6 7 8 9 10]}
%           ADDUCT: {'Na+'  [7]}
%             CHAR: {[]  [10]}
%
% Children function:
% GETSEQOPTIONS, SPLIT
%

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2019, Research Foundation for State University of New York. All rights reserved
%


%% locate content in quotation mark
[optcontent,optcontentstart,optcontentend] = regexp(GLYseq,DrawGlycanPara.regexp_optionvalue,'match','start','end');
tempGPseq = GLYseq;
if ~isempty(optcontent)  % block customization info ("noise")
    for i = 1:length(optcontentstart)
        tempGPseq(optcontentstart(i):optcontentend(i)) = blanks(optcontentend(i) - optcontentstart(i) + 1);
    end
end
specialoptions = [];
%%  curly bracket handling
indtemp = strfind(tempGPseq,'{');
if any(indtemp)
    levelindex = zeros(2,length(tempGPseq));
    indtemp2 = strfind(tempGPseq,'}');
    levelindex(1,indtemp) = 1;
    levelindex(1,indtemp2) = -1;
    levelindex(2,1) = 1;
    for i = 2:size(levelindex,2)
        levelindex(2,i) = levelindex(2,i-1) + levelindex(1,i);
    end
    numsection = length(indtemp)/2+1;
    sectionpos = zeros(numsection,2);
    sectionend = find(levelindex(2,:) == 0)';
    sectionpos(:,2) = [sectionend(1:numsection-1);sectionend(end)];
    sectionpos(1,1) = 1;
    sectionpos(2:end,1) = sectionpos(1:end-1,2)+1;  % starting, ending position of curly bracket units, last one is main sturcture
    cbitems = cell(size(sectionpos,1)-1,2);
    for i = 1:size(cbitems,1)
        cbitems{i,1} = GLYseq(sectionpos(i,1)+1:sectionpos(i,2)-1);
    end  % original plan was to develop nested curly bracket representation, cancelled due to complexity
    
    mainstructGPseq = tempGPseq(sectionpos(end,1):sectionpos(end,2));
    indtempmain = strfind(mainstructGPseq,'{');
    levelindex = zeros(1,length(mainstructGPseq));
    indtemp2main = strfind(mainstructGPseq,'}');
    levelindex(indtempmain) = 1;
    levelindex(indtemp2main) = -1;
    opencbrep = [];
    cbpair = [];
    ind = 1;
    while ind <= length(levelindex)
        if levelindex(ind) == 1
            opencbrep = [opencbrep;ind];
        elseif levelindex(ind) == -1
            cbpair = [cbpair;[opencbrep(end),ind]];
            opencbrep(end) = [];
        end
        ind = ind + 1;
    end
    cbpair = sortrows(cbpair,1);
    %% original ver.
    %     mainstructGPseq = regexprep(mainstructGPseq,'[{}]',' ');
    %     if strcmp(mainstructGPseq(end),' ')
    %         mainstructGPseq(end) = '';
    %     end
    %     [~,MonoIndex,~,~] = Split(mainstructGPseq);
    %     cbcontent = cell(size(cbpair,1),1);
    %     for i = 1:size(cbpair,1)
    %         cbcontent{i} = find(MonoIndex < cbpair(i,2) & MonoIndex >= cbpair(i,1));
    %     end
    %     for i = 1:size(cbitems,1)
    %         cbitems{i,2} = cbcontent{i};  % curly bracket coverage area, shown using monosac. index number
    %     end
    %% original ver. end
    
    %% new ver.
    % Cruly brackets SOMETIMES inteferes monosac. string position detection
    % because they are considered as a monosaccharide in func. "Split". This
    % fix aims to fix this issue
    cbpos = regexp(mainstructGPseq,'[{}]');
    mainstructGPseq = regexprep(mainstructGPseq,'[{}]','');
    if strcmp(mainstructGPseq(end),' ')
        mainstructGPseq(end) = '';
    end
    [~,MonoIndex,~,~] = Split(mainstructGPseq);
    for i = 1:length(MonoIndex)
        MonoIndex(i) = MonoIndex(i) + sum(MonoIndex(i) > cbpos);
    end
    cbcontent = cell(size(cbpair,1),1);
    for i = 1:size(cbpair,1)
        cbcontent{i} = find(MonoIndex < cbpair(i,2) & MonoIndex >= cbpair(i,1));
    end
    for i = 1:size(cbitems,1)
        cbitems{i,2} = cbcontent{i};  % curly bracket coverage area, shown using monosac. index number
    end
    
    
    %% new ver. end
    
    specialoptions.CURLYBRACKET = cbitems;
    cleanGLYseq = GLYseq(sectionpos(end,1):sectionpos(end,2));
    cleanGLYseq([indtempmain,indtemp2main]) = '';
else
    cleanGLYseq = GLYseq;
end

%% other contents
[otheroptions,cleanerseq] = getseqoptions(cleanGLYseq,'G');
if ~isempty(otheroptions)  % combine option info.
    otheroptionnames = fieldnames(otheroptions);
    for i = 1:length(otheroptionnames)
        specialoptions.(otheroptionnames{i}) = otheroptions.(otheroptionnames{i});
    end
end
end
