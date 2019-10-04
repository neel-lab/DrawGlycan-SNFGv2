function Sgpout = IUPAC2Sgp(IUPAC,varargin)
%IUPAC2Sgp: Conver IUPAC to Sgp1.0 or Sgp2.0 format.
%
% Syntax:
% Sgpout = IUPAC2Sgp(IUPAC)
% Sgpout = IUPAC2Sgp(IUPAC,version)  default sgp version=1.0
%
% Input:
% IUPAC: a type of glycan nomenclature in linear format.
%
% Output:
% Sgpout: glycan in sgp format
%
% Example:
% M5 = 'Man(a1-4)[Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-?)';
% Sgpout = IUPAC2Sgp(M5)
%
% Sgpout =
% {n{n{h{h}{h{h}{h}}}}}
%
% Sgpout = IUPAC2Sgp(M5,2)
%
% Sgpout =
% {(b1-?)GlcNAc{(b1-4)GlcNAc{(b1-4)Man[{(a1-3)Man}]{(a1-6)Man[{(a1-6)Man}]{(a1-4)Man}}}}}
%
% Bibrabch = 'NeuAc(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-6)[NeuAc(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-?)';
% Sgpout = IUPAC2Sgp(Bibrabch,1)
%
% Sgpout =
% {n{f}{n{h{h{n{h{s}}}}{h{n{f}{h{s}}}}}}}
%
% Date Lastly Updated: 08/30/19

% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2019, Research Foundation for State University of New York. All rights reserved
%


% Input handle
narginchk(1,2);
version = 1;

if(isempty(IUPAC))
    Sgpout = [];
    return
end

if(~isempty(varargin))
    version = varargin{1};
end

% Split the whole sequence into several cells, each cell contains one
% monosaccharide or bracket
[MonoStr,MonoIndex,BracStr,BracIndex] = Split(IUPAC);

% Convert to Sgp format
if(version==1)
    Sgpout = Convert2Sgp(MonoStr,MonoIndex,BracStr,BracIndex);
elseif(version==2)
    Sgpout = Convert2Sgp2(MonoStr,MonoIndex,BracStr,BracIndex);
    %continue
end
end

function Sgpout = Convert2Sgp(MonoStr,MonoIndex,BracStr,BracIndex)
Hex = {'Glc[^AN]','Man[^AN]','Gal[^AN]','Gul[^AN]','Alt[^AN]','All[^AN]',...
    'Tal[^AN]','Ido[^AN]'};
HexNAc = {'GlcNAc','ManNAc','GalNAc','GulNAc','AltNAc','AllNAc',...
    'TalNAc','IdoNAc'};
HexN = {'GlcN[^A]','ManN[^A]','GalN[^A]','GulN[^A]','AltN[^A]',...
    'AllN[^A]','TalN[^A]','IdoN[^A]'};
HexA = {'GlcA','ManA','GalA','GulA','AltA','AllA','TalA','IdoA'};
DHexNAc = {'QuiNAc','RhaNAc','FucNAc'};
DHex = {'Qui[^N]','Rha[^N]','6dAlt','6dTal','Fuc[^N]'};
DiDHex = {'Oli','Tyv','Abe','Par','Dig','Col'};
Pen = {'Ara','Lyx','Xyl','Rib'};
Othersin3 = {'Kdn','NeuGc','NeuAc','Neu','Bac','LDMan',...
    'Kdo','Dha','DDMan','MurNAc','MurNGc','Mur','Api','Fru',...
    'Tag','Sor','Psi'};
Othersin1 = {'k','g','s','u','b','L',...
    'o','y','D','A','M','m','i','c','t','r','i'};
Sgpout     = [];
WholeSeq   = sort(cat(2,MonoIndex,BracIndex));
Startpoint = length(WholeSeq);
BracLvl    = 0;
Brackpos   = zeros(length(BracIndex)/2,1);
for i = 1:length(WholeSeq)
    ithIndex = WholeSeq(Startpoint);
    if(~isempty(find(ithIndex-MonoIndex==0,1)))
        ithCellStr = MonoStr{find(ithIndex-MonoIndex==0,1)};
        if(sum(cellfun(@(x)matchMono(x,ithCellStr),Hex)==1))
            Sgpunit = '{h';
            Sgpout  = [Sgpout Sgpunit];
        elseif(sum(cellfun(@(x)matchMono(x,ithCellStr),HexNAc)==1))
            Sgpunit = '{n';
            Sgpout  = [Sgpout Sgpunit];
        elseif(sum(cellfun(@(x)matchMono(x,ithCellStr),HexN)==1))
            Sgpunit = '{e';
            Sgpout  = [Sgpout Sgpunit];
        elseif(sum(cellfun(@(x)matchMono(x,ithCellStr),HexA)==1))
            Sgpunit = '{a';
            Sgpout  = [Sgpout Sgpunit];
        elseif(sum(cellfun(@(x)matchMono(x,ithCellStr),DHexNAc)==1))
            Sgpunit = '{N';
            Sgpout  = [Sgpout Sgpunit];
        elseif(sum(cellfun(@(x)matchMono(x,ithCellStr),DHex)==1))
            Sgpunit = '{f';
            Sgpout  = [Sgpout Sgpunit];
        elseif(sum(cellfun(@(x)matchMono(x,ithCellStr),DiDHex)==1))
            Sgpunit = '{d';
            Sgpout  = [Sgpout Sgpunit];
        elseif(sum(cellfun(@(x)matchMono(x,ithCellStr),Pen)==1))
            Sgpunit = '{x';
            Sgpout  = [Sgpout Sgpunit];
        elseif(sum(cellfun(@(x)matchMono(x,ithCellStr),Othersin3)==1))
            MonoId = find(cellfun(@(x)matchMono(x,ithCellStr),Othersin3));
            SgpSymbol = Othersin1{MonoId};
            Sgpunit   = ['{' SgpSymbol];
            Sgpout    = [Sgpout Sgpunit];
        end
    else
        BracPos   = find(ithIndex-BracIndex==0,1);
        BracSym   = BracStr(BracPos);
        if(strcmp(BracSym,']'))
            startpos = find(ithIndex-WholeSeq==0,1);
            BracLvl  = BracLvl+1;
            Brackpos(BracLvl) = length(Sgpout)+1;
        else
            StrLvl       = Sgpout(Brackpos(BracLvl):end);
            MonoStart    = length(strfind(StrLvl,'{'));
            Monoend      = length(strfind(StrLvl,'}'));
            Mononum      = MonoStart - Monoend;
            for j = 1 : Mononum
                Sgpout       = [Sgpout '}'];
            end
            BracLvl      = BracLvl-1;
            if(BracLvl==0)
                Brackpos   = zeros(length(BracIndex)/2,1);
            end
        end
    end
    Startpoint = Startpoint-1;
end
MonoStart    = length(strfind(Sgpout,'{'));
Monoend      = length(strfind(Sgpout,'}'));
for i = 1 : MonoStart-Monoend
    Sgpout  = [Sgpout '}'];
end
end

function Sgpout = Convert2Sgp2(MonoStr,MonoIndex,BracStr,BracIndex)
Sgpout     = [];
WholeSeq   = sort(cat(2,MonoIndex,BracIndex));
Startpoint = length(WholeSeq);
BracLvl    = 0;
Brackpos   = zeros(length(BracIndex)/2,1);
for i = 1:length(WholeSeq)
    ithIndex = WholeSeq(Startpoint);
    if(any(ithIndex==MonoIndex))
        ithCellStr = MonoStr{find(ithIndex-MonoIndex==0,1)};
        % Check linkage information
        ithCellStr = ChkLinkInfo(ithCellStr);
        Sgpunit = ['{' ithCellStr];
        Sgpout  = [Sgpout Sgpunit];
    else
        BracPos   = find(ithIndex-BracIndex==0,1);
        BracSym   = BracStr(BracPos);
        if(strcmp(BracSym,']'))
            startpos = find(ithIndex-WholeSeq==0,1);
            BracLvl  = BracLvl+1;
            Sgpout   = [Sgpout '['];
            Brackpos(BracLvl) = length(Sgpout);
        else
            StrLvl       = Sgpout(Brackpos(BracLvl):end);
            MonoStart    = length(strfind(StrLvl,'{'));
            Monoend      = length(strfind(StrLvl,'}'));
            Mononum      = MonoStart - Monoend;
            for j = 1 : Mononum
                Sgpout       = [Sgpout '}'];
            end
            Sgpout           = [Sgpout ']'];
            BracLvl      = BracLvl-1;
            if(BracLvl==0)
                Brackpos   = zeros(length(BracIndex)/2,1);
            end
        end
    end
    Startpoint = Startpoint-1;
end
MonoStart    = length(strfind(Sgpout,'{'));
Monoend      = length(strfind(Sgpout,'}'));
for i = 1 : MonoStart-Monoend
    Sgpout  = [Sgpout '}'];
end
end

function isMatch = matchMono(Mono,Monolist)
isMatch = 0;
if(~isempty(regexp(Monolist,Mono,'match')))
    isMatch = 1;
end
end

function MonoBlock = ChkLinkInfo(MSinput)
LinkExp = '\(.*\)';
LinkInfo = regexp(MSinput,LinkExp,'Match');
if(~isempty(LinkInfo))
    MonoBlock = regexprep(MSinput,LinkExp,'');
    StdLinkExp = '\([ab][\d\?][-,][\d\?].*\)';
    if(~isempty(regexp(LinkInfo{1},StdLinkExp,'Once')))  % if link info is in standard format
        newLinkInfo = LinkInfo{1}(1:end-1);
        [Contstart,Contend] = regexp(newLinkInfo,'["''].*?["'']','start','end');
        newLinkInfo(Contstart) = '"';
        newLinkInfo(Contend) = '"';
        if ~isempty(Contstart)
            partA = strrep(newLinkInfo(1:min(Contstart) - 1),',','-');
            newLinkInfo = [partA,newLinkInfo(min(Contstart):end)];
        else
            newLinkInfo = strrep(newLinkInfo,',','-');
        end
    else  % unorthodox link format
        newLinkInfo = linkrepair(LinkInfo);
    end
    LinkInfo = [newLinkInfo,')'];
else  % no link info or damaged parenthesis
    openparenthesis = strfind(MSinput,'(');
    if any(openparenthesis)
        MonoBlock = MSinput(1:openparenthesis(1)-1);
        LinkInfo = MSinput(openparenthesis(1):end);
        if strcmp(LinkInfo(end),'-') || strcmp(LinkInfo(end),',')
            LinkInfo = [LinkInfo,'?)'];
        else
            LinkInfo = [LinkInfo,')'];
        end
    else
        MonoBlock = MSinput;
        LinkInfo = '(??-?';
    end
    LinkInfo = linkrepair({LinkInfo});
    LinkInfo = [LinkInfo,')'];
end
MonoBlock = [LinkInfo,MonoBlock];
end

function newLinkInfo = linkrepair(LinkInfo)
LinkInfo = LinkInfo{1}(2:end-1);
if isempty(LinkInfo)  % no link/frag info
    newLinkInfo = '(??-?';
    fraginfo = '';
else  % have link/frag info
    fragstart = regexp(LinkInfo,'-R|-NR|-U|-D','start');
    if isempty(fragstart)
        fraginfo = '';
        fragstart = length(LinkInfo) + 1;
    else
        fraginfo = [' ',LinkInfo(fragstart:end)];
    end
    glybond = LinkInfo(1:min(fragstart)-1);
    glybond = strrep(glybond,',','-');
    if isempty(glybond)
        newLinkInfo = '(??-?';
    elseif length(glybond) == 1  % bond info is question mark/1 letter/1 digit
        glybondascii = double(LinkInfo(1:min(fragstart)-1));
        if glybondascii == 63 %% question mark
            newLinkInfo = '(??-?';
        elseif glybondascii >47 && glybondascii < 58  % 0-9
            newLinkInfo = ['(??-',LinkInfo(1:min(fragstart)-1)];
        elseif (glybondascii >64 && glybondascii < 91) || (glybondascii >96 && glybondascii < 123)  % a-z and A-Z
            newLinkInfo = ['(',LinkInfo(1:min(fragstart)-1),'?-?'];
        else  % unrecognizable
            newLinkInfo = '(??-?';
        end
    else  % bond info is longer
        carbon = regexp(glybond,'[a-z]','once');
        qmark = regexp(glybond,'?','once');
        pairnum = regexp(glybond,'[\d/\\]+\-[?\d/\\]+','once','match');
        [anynum,anynumstart] = regexp(glybond,'[\d/\\]+','once','match','start');
        if isempty(carbon)  % long bond info, no letter, question mark/number only
            if isempty(qmark)  % long bond info, no letter, no question mark
                if isempty(pairnum)  %  long bond info, no letter, no question mark, number not in pair
                    newLinkInfo = ['(??-',glybond];
                else  %  long bond info, no letter, no question mark, number in pair
                    newLinkInfo = ['(?',glybond];
                end
            elseif strcmp(glybond,'??-?')
                newLinkInfo = '(??-?';
            elseif ~isempty(anynum) && ~isempty(qmark)  % long bond info, no letter, question mark/number both exist
                if isempty(pairnum)  % long bond info, no letter, question mark exist, number not in pair
                    if qmark < anynumstart  % (?#)
                        newLinkInfo = ['(??-',anynum];
                    else  % (#?)
                        newLinkInfo = ['(?',anynum,'-?'];
                    end
                else % long bond info, no letter, question mark exist, number in pair
                    newLinkInfo = ['(?',pairnum];
                end
            else  % long bond info, unrecognizable
                newLinkInfo = '(??-?';
            end
        else  % long bond info, have carbon letter
            if isempty(qmark)  % long bond info, have carbon letter, no question mark
                if isempty(pairnum)  % long bond info, have carbon letter, no question mark, number not in pair
                    newLinkInfo = ['(',glybond(carbon),'?-',anynum];
                else  % long bond info, have carbon letter, no question mark, number in pair
                    newLinkInfo = ['(',glybond(carbon),pairnum];
                end
            else  % long bond info, have carbon letter, have question mark
                newLinkInfo = ['(',glybond(carbon),'?-?'];
            end
        end
    end
end
newLinkInfo = [newLinkInfo,fraginfo];
end