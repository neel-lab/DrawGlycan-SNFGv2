function [gly,pep,glypos,PTMid] = distggp(seqinput,inputformat)
% DISTGGP: Parse the input IUPAC sequence to get glycan(s) and peptide
%
% Syntax:
% [gly,pep,glypos,PTMid] = distggp(seqinput,inputformat)
%
% Input:
% seqinput: string, the sequence of the % glycan/glycopeptide to be identified.
% inputformat: string, the format of input sequence.
%
% Output:
% gly: n x 1 cell array of strings, the glycan sequence of glycan/glycopeptide.
% pep: string, the peptide sequence of glycopeptide, this output is blank
% when input is glycan only.
% glypos: n x 1 numerical array, position where glycan is attached to peptide
% backbone, counting from N-terminus.
% PTMid: n x 1 numerical array, 0 means non-glycan PTM, 1 means glycan.
% This applies to SGP1.0 format only, this value is 1 for all other
% formats.
%
% Note:
% Input must be in a string.
% Input must be consists of pure glycan(s) or contains only one
% glycopeptide. Mixture of glycan and glycopeptide is prohibited.
% If there is any, any lowercase letter in sequence outside
% brackets, the whole sequence will be treated as a glycan.
% Support for more input formats will be added in the future.
%
% Example:
% 1. IUPAC format
% [gly,pep,glypos,PTMid] = distggp('FKT[(??-?)GalNAc[(??-?)Fuc](??-?)Gal]GTK','IUPAC')
%
% gly =
%
%   1×1 cell array
%
%     {'(??-?)GalNAc[(??-?)Fuc](??-?)Gal'}
%
%
% pep =
%
%     'FKTGTK'
%
%
% glypos =
%
%      3
%
%
% PTMid =
%
%      1
%
% 2. SGP1.0 format
% [gly,pep,glypos,PTMid] = distggp('FKT{n{f}{h}}GM<o>K','SGP1')
%
% gly =
%
%   2×1 cell array
%
%     {'{n{f}{h}}'}
%     {'<o>'      }
%
%
% pep =
%
%     'FKTGMK'
%
%
% glypos =
%
%      3
%      5
%
%
% PTMid =
%
%      1
%      0
%
% Children function:
% N/A

% ver 1.1, added output "glypos", "PTMid"

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2019, Research Foundation for State University of New York. All rights reserved
%


switch upper(inputformat)
    case 'SGP1'
        gly = {};
        glypos = [];
        [p,g,m] = breakGlyPep(seqinput);  % use BREAKGLYPEP to get pep, gly and mod
        if ~isempty(g)
            gly = {g.struct};
            glypos = [g.pos];
        end
        pep = p.pep;
        if ~isempty(m)
            gly = [gly,{m.struct}];  % combine gly and mod
            glypos = [glypos,[m.pos]];
        end
        gly = gly(:);
        [glypos,ind] = sort(glypos(:));  % displayed in order of appearance in GlyPep
        gly = gly(ind);
        PTMidG = ones(length(g),1);
        PTMidNG = zeros(length(m),1);
        PTMid = [PTMidG;PTMidNG];
        PTMid = PTMid(ind);
    case 'SGP2'
        
    case 'IUPAC'
        thissgp = seqinput;
        thissgp = strtrim(thissgp);
        % locate and replace option values with blanks, in case its
        % contents intefere with regexp
        [comstart,comend] = regexp(thissgp,DrawGlycanPara.regexp_optionvalue,'start','end');
        tempthissgp = thissgp;
        %% block user-customized info using whitespaces
        for i = 1:length(comstart)
            tempthissgp(comstart(i):comend(i)) = ' ';
        end
        %% identify glycan from glypeptide by counting brackets
        opensqbr = strfind(tempthissgp,'[');
        closesqbr = strfind(tempthissgp,']');
        [sqbrpos,ind] = sort([opensqbr closesqbr]);
        sqbrcounter = [ones(size(opensqbr)),-ones(size(closesqbr))];
        sqbrcounter = sqbrcounter(ind);
        sqbrcounterind = arrayfun(@(x) sum(sqbrcounter(1:x)),1:length(sqbrcounter));
        sqbrecorder = [0 sqbrcounterind] == 0;
        sqbrecorder = [sqbrpos(sqbrecorder(1:end-1));sqbrpos(sqbrcounterind == 0)];
        % "sqbrecorder" is a numerical array, initialized with zeroes, "["
        % adds 1 to all numbers behind, "]" minuses 1 to all numbers
        % behind. The goal is to find out where number falls to zero, which
        % means a glycan sequence ends here
        glypos = zeros(size(sqbrecorder,2),1);
        %% memorizing glycan position
        pepseq = thissgp;  % now we don't know if it's a pep or gly, call it pepseq for now
        glyseq = cell(size(sqbrecorder,2),1);
        for i = size(sqbrecorder,2):-1:1
            pepseq(sqbrecorder(1,i):sqbrecorder(2,i)) = '';
            glyseq{i} = thissgp(sqbrecorder(1,i)+1:sqbrecorder(2,i)-1);
        end
        %% peptide options: e.g. fragmentation
        [pepoptstart,pepoptend] = regexp(pepseq,'\((.*?)\)|<(.*?)>','start','end');
        pepseq2 = pepseq;  % Need this duplicate for removing lowercase letters in pep. fragmentation info.
        for i = length(pepoptstart):-1:1
            pepseq2(pepoptstart(i):pepoptend(i)) = '';
        end
        if any(isstrprop(pepseq2,'lower'))  % There is at least one lowercase letter, this is a glycan. Correct input is very important for proper functioning
            gly = seqinput;
            pep = '';
        else  % input is glycopeptide
            gly = glyseq;
            pep = strtrim(pepseq);
        end
        if ~isempty(pep)
            %% Get glycosylation position after eliminating peptide fragment info
            glyseqcompen = 0;
            for i = 1:size(sqbrecorder,2)
                glypos(i) = sqbrecorder(1,i) - 1 - glyseqcompen;
                for j = length(pepoptstart):-1:1
                    if any(glypos >= pepoptend(j))
                        glypos(glypos >= pepoptend(j)) = glypos(glypos >= pepoptend(j)) - (pepoptend(j) - pepoptstart(j) + 1);
                    end
                end
                glyseqcompen = glyseqcompen + sqbrecorder(2,i) - sqbrecorder(1,i) + 1;
            end
        end
        PTMid = ones(length(glyseq),1);
end
end