function [pepMat,glyMat,modMat]=breakGlyPep(SmallGlyPep)
%BREAKGLYPEP: Parse glycopeptide (SmallGlyPep) to return component peptide 
% (pepMat), glycan (in glyMat) and other PTMs (in modMat) in structure format
% 
% Syntax:
%     [pepMat,glyMat,modMat]=breakGlyPep(SmallGlyPep)
%  
% Input: SmallGlyPep (Can be glycopeptide, peptide or glycan)
% 
% Output: Structures pepMat (peptide), glyMat (glycan) and modMat (PTMs) 
%   -- pepMat stores position (.pos) data and the peptide sequence (.pep)
%   -- glyMat stores position (.pos) on base peptide, glycan structure (.struct), 
%      and glycan length (.len)
%   -- modMat stores PTM modification (.struct) and position (.pos) on base peptide
% 
% Examples:
% Example 1: for glycopeptide
% >> SmallGlyPep='GYLN{n{n{h{h{h{h}}}{h{h{h}}{h{h}}}}}}CT{n{h{s}}{n{h{s}{f}}}}R';
%    [pepMat,glyMat,modMat]=breakGlyPep(SmallGlyPep)
% Answer:
%       pepMat =    pep: 'GYLNCTR'
%                   pos: [1 2 3 4 38 39 61]
%       glyMat =    1x2 struct array with fields:
%                       pos
%                       struct
%                       len
%       modMat =   []
%    >> glyMat(1).pos
%       ans =   4
%    >> glyMat(1).struct
%       ans =   {n{n{h{h{h{h}}}{h{h{h}}{h{h}}}}}}
%    >> glyMat(1).len
%       ans =   11
%
% Example 2: for glycopeptide (with glycan and another modification on the
% same amino acid)
% >> SmallGlyPep='GYLN{n{n{h{h{h{h}}}{h{h{h}}{h{h}}}}}}CT<s>{n{h{s}}{n{h{s}{f}}}}R';
%    [pepMat,glyMat,modMat]=breakGlyPep(SmallGlyPep)
% Answer:
%       pepMat =    pep: 'GYLNCTR'
%                   pos: [1 2 3 4 38 39 64]
%       glyMat =    1x2 struct array with fields:
%                       pos
%                       struct
%                       len
%       modMat =    pos: 6
%                   struct: '<s>'
%    >> glyMat(2).pos
%       ans =   6
%    >> glyMat(2).struct
%       ans =   {n{h{s}}{n{h{s}{f}}}}
%    >> glyMat(2).len
%       ans =   7
%
% Example 3: glycan only
% >> SmallGlyPep='{n{h{s}}{n{h{s}{f}}}}';
%    [pepMat,glyMat,modMat]=breakGlyPep(SmallGlyPep)
% Answer:
%       pepMat =    pep: ''
%                   pos: 0
%       glyMat =    pos: 0
%                   struct: '{n{h{s}}{n{h{s}{f}}}}'
%                   len: 7
%       modMat =    []
%
% Example 4: peptide only
% >> SmallGlyPep='GYM<o>KNCT<s>';
%    [pepMat,glyMat,modMat]=breakGlyPep(SmallGlyPep)
% Answer:
%       pepMat =    pep: 'GYMKNCT'
%                   pos: [1 2 3 7 8 9 10]
%       glyMat =    []
%       modMat =    1x2 struct array with fields:
%                       pos
%                       struct
%    >> modMat(1).pos
%       ans =   3
%    >> modMat(1).struct
%       ans =   <o>
%
% Example 5: for complete glycopeptide (for custom glycan and non-glycan fragmentation)
% >> SmallGlyPep='GYLN{n{n{h{h{h{162.1}}}{h{h{h}}{h{h}}}}}}CT<+96>{n{h{s}}{n{h{s}{f}}}}R';
%    [pepMat,glyMat,modMat]=breakGlyPep(SmallGlyPep)
% Answer:
%       pepMat =    pep: 'GYLNCTR'
%                   pos: [1 2 3 4 42 43 70]
%       glyMat =    1x2 struct array with fields:
%                       pos
%                       struct
%                       len
%       modMat =    pos: 6
%                   struct: '<+96>'
%    >> glyMat(1).pos
%       ans =   4
%    >> glyMat(1).struct
%       ans =   {n{n{h{h{h{162.1}}}{h{h{h}}{h{h}}}}}}
%    >> glyMat(1).len
%       ans =   11
% 
%See also joinGlyPep, compileFrags, glycanFrag, multiSGPFrag, UQFragIon. 

% Author: Sriram Neelamegham
% Date Lasty updated: 8/11/14 by Gang Liu

%
% DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
%(c) 2019, Research Foundation for State University of New York. All rights reserved
%


 % create a data structure glyMat, pepMat and modMat to store related info
glyMat=[];        
pepMat=[];
modMat=[];

% protein sequence
% aa1letcharexpr = '[ARNDCQEGHILKMFPSTWYV]';
pepPos=regexp(SmallGlyPep, '[ARNDCQEGHILKMFPSTWYV]');
if(isempty(pepPos))
    pepPos = 0;
    pep    =  '';
else
    pep = SmallGlyPep(pepPos);
end
pepMat.pep = pep;
pepMat.pos = pepPos;

% modification 
[modpos,modendingpos] = regexp(SmallGlyPep,'<[+-\d\.a-zA-Z_0-9]+>');
for i = 1 : length(modpos)
    modMat(i).struct = SmallGlyPep(modpos(i):modendingpos(i));
    newglypep        = SmallGlyPep(1:modpos(i)-1);
    newpepPos        = regexp(newglypep, '[ARNDCQEGHILKMFPSTWYV]');
    modMat(i).pos    = length(newpepPos);
end

% glycan
if(~isempty(pep))
    [glypos,glyendingpos] = regexp(SmallGlyPep,'{[+\d.a-z{}]+}(?=[A-Z][<[+-\da-z]+>]*)|(?<=[A-Za-z<>\[\]+-\d]*){[+\d.a-z{}]+}');
    for i = 1 : length(glypos)
        glyMat(i).struct = SmallGlyPep(glypos(i):glyendingpos(i));
        newglypep        = SmallGlyPep(1:glypos(i)-1);
        newpepPos        = regexp(newglypep, '[ARNDCQEGHILKMFPSTWYV]');
        glyMat(i).pos    = length(newpepPos);
        glyMat(i).len    = length(strfind(glyMat(i).struct,'{'));    
   end
else
    glyMat(1).struct = SmallGlyPep;
    glyMat(1).pos    = 0;
    glyMat(1).len    = length(strfind(SmallGlyPep,'{'));        
end

end