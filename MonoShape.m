classdef MonoShape
    %% MONOSHAPE: shape definition for monosaccharides
    % Note: 6dAlt and 6dTal are presented as dAlt6 and dTal6 to avoid variable
    % naming issues.
    
    %
    % DrawGlycan authors: Kai Cheng, Yusen Zhou, Sriram Neelamegham
    %(c) 2019, Research Foundation for State University of New York. All rights reserved
    %
    
    
    properties (Constant)
        % Hexose
        Hex = 'Filled Circle';
        Hexose = 'Filled Circle';
        Glc = 'Filled Circle';
        Man = 'Filled Circle';
        Gal = 'Filled Circle';
        Gul = 'Filled Circle';
        Alt = 'Filled Circle';
        All = 'Filled Circle';
        Tal = 'Filled Circle';
        Ido = 'Filled Circle';
        
        % HexNAc
        HexNAc = 'Filled Square';
        GlcNAc = 'Filled Square';
        ManNAc = 'Filled Square';
        GalNAc = 'Filled Square';
        GulNAc = 'Filled Square';
        AltNAc = 'Filled Square';
        AllNAc = 'Filled Square';
        TalNAc = 'Filled Square';
        IdoNAc = 'Filled Square';
        
        % HexN
        Hexosamine = 'Crossed Square';
        GlcN = 'Crossed Square';
        ManN = 'Crossed Square';
        GalN = 'Crossed Square';
        GulN = 'Crossed Square';
        AltN = 'Crossed Square';
        AllN = 'Crossed Square';
        TalN = 'Crossed Square';
        IdoN = 'Crossed Square';
        
        % HexA
        Hexuronate = 'Divided Diamond';
        GlcA = 'Divided Diamond';
        ManA = 'Divided Diamond';
        GalA = 'Divided Diamond';
        GulA = 'Divided Diamond';
        AltA = 'Divided Diamond Inverted';
        AllA = 'Divided Diamond';
        TalA = 'Divided Diamond';
        IdoA = 'Divided Diamond Inverted';
        
        % Deoxyhexose
        Deoxyhexose = 'Filled Triangle';
        Qui   = 'Filled Triangle';
        Rha   = 'Filled Triangle';
        Fuc   = 'Filled Triangle';
        
        % DeoxyhexNAc
        DeoxyhexNAc = 'Divided Triangle';
        QuiNAc   = 'Divided Triangle';
        RhaNAc   = 'Divided Triangle';
        FucNAc   = 'Divided Triangle';
        
        % Di-Deoxyhexose
        Dideoxyhexose = 'Flat Rectangle';
        Oli   = 'Flat Rectangle';
        Tyv   = 'Flat Rectangle';
        Abe   = 'Flat Rectangle';
        Par   = 'Flat Rectangle';
        Dig   = 'Flat Rectangle';
        Col   = 'Flat Rectangle';
        
        % Pentose
        Pentose = 'Filled Star';
        Ara   = 'Filled Star';
        Lyx   = 'Filled Star';
        Xyl   = 'Filled Star';
        Rib   = 'Filled Star';
        
        % Deoxynonulosonate
        Deoxynonulosonate = 'Filled Diamond';
        Kdn    = 'Filled Diamond';
        Neu5Ac = 'Filled Diamond';
        Neu5Gc = 'Filled Diamond';
        Neu    = 'Filled Diamond';
        Sia      = 'Filled Diamond';
        
        % Dideoxynonulosonate
        Dideoxynonulosonate = 'Flattened Diamond';
        Pse = 'Flattened Diamond';
        Leg = 'Flattened Diamond';
        Aci = 'Flattened Diamond';
        
        % Unknown
        Unknown = 'Flat Hexagon';
        Bac      = 'Flat Hexagon';
        LDManHep = 'Flat Hexagon';
        Kdo      = 'Flat Hexagon';
        Dha      = 'Flat Hexagon';
        DDManHep = 'Flat Hexagon';
        MurNAc   = 'Flat Hexagon';
        MurNGc   = 'Flat Hexagon';
        Mur      = 'Flat Hexagon';
        
        % Assigned
        Assigned = 'Pentagon';
        Api = 'Pentagon';
        Fru = 'Pentagon';
        Tag = 'Pentagon';
        Sor = 'Pentagon';
        Psi = 'Pentagon';
        NonGly = 'Pentagon';
        Blank = 'None';
    end
end