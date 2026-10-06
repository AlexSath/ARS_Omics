from .color_palettes import black_red_pastel_palette
from .color_palettes import black_blue_pastel_palette
from .color_palettes import black_green_pastel_palette
from .color_palettes import black_yellow_pastel_palette

GENES_OF_INTEREST = [
    "DNMT3A",
    "DNMT3AP1",
    "DNMT3L",
    "DNMT3B",
    "CD47",
    "CD55",
    "CD81",
    "CD151",
    "B2M",
    "KIT"
]

SCREEN_GENES = {
    "LX2 hits": [
        "SMAD3", # Daniel LX-2 screen hit
        "NFX1", # Daniel LX-2 screen hit
        "GPX4", # Daniel LX-2 screen hit
        "ESYT2", # Daniel LX-2 screen hit
        "LAMTOR1", # Daniel LX-2 screen hit
        "ZNF259", # Daniel LX-2 screen hit
        "PFKFB3", # Daniel LX-2 screen hit
        "TGFBR2", # Daniel LX-2 screen hit
        "C16orf87", # proposed renamed HDIP for HDAC-interacting protein, Daniel LX-2 screen hit
    ],
    "LX2 below thresh": [ # Daniel LX-2 screen below threshold
        "SMAD4", 
        "FOXO3",
        "HOXA3",
        "MEOX1", 
        "BRD2", 
        "MYH9",
        "CCN2", # aka Ctgf
        "TGFB1",
        "ZSCAN21",
        "ZNF559",
        "ZNF330",
    ],
    "LX2 screen enriched": [
        "MTHFD2", # Daniel screen enriched
        "ZNF567", # Daniel screen enriched
        "NUPL1", # Daniel screen enriched
    ],
}

FIBROSIS_GENES = [
    "TILAM",
    "POSTN",
    "JUNB",
    "GATA4",
    "GATA6",
    "TP73",
    "FOXC2",
    "SERPINE1",
    "ACTA1",
    "ACTA2",
    "YAP1",
    "LTBP2",
    "WWTR1", # TAZ, another supposed mechanosensitive transcription factor
    "TGFBR1",
    "TGFBR2",
    "BRD3",
    "BRD4",
    "BRDT",
    "COL1A1",
    "SPI1", # see "PU.1 controls fibroblast..." (2019)
    "SMAD2",
]

CONSTRUCT_NAME_DICT = {
    "zim3": "dCas9-ZIM3",
    "ZIM3": "dCas9-ZIM3",
    "ARS001": "dCas9-ZIM3",
    "pARS001": "dCas9-ZIM3",
    "pARS1": "dCas9-ZIM3",
    "pCH45": "multiAsCas12a-KOX1",
    "pARS002": "multiAsCas12a-ZIM3",
    "pRA2": "hyperLbCas12a-KOX1",
    "pCH95": "H3t-D3L-dCas9",
    "pCH96": "H3t-D3L-dCas9-KOX1",
    "pCH97": "H3t-KOX1-D3L-dCas9",
    "HEK-CLTA": "HEK",
    "sgCLTA-1X": "sgCLTA",
    "sgCLTA-20X": "sgCLTA",
    "Mock": "LX2",
    "Ctrl": "Unstained",
    "LX-2": "LX2",
    "LX2": "LX2",
}

GUIDE_NAME_DICT = {
    "_HEK_" : "CLTA-",
    "_HEK-CLTA_" : "CLTA+",
    "-1X" : "1X",
    "-20X" : "20X",
    "Mock" : "LX-2",
    "Parental" : "",
    "Unstained" : "Ab-",
    "Stained" : "Ab+",
    "sgNT": "NT",
    "gLG20": "NT",
    "no1ary" : "1°Ab-",
    "no1ab" : "1°Ab-",
    "w1ary" : "1°Ab+",
    "w1ab": "1°Ab+",
    "Stained" : "1°Ab+",
    "DC01" : "CD81",
    "sgDC01" : "CD81",
    "DC02" : "B2M",
    "sgDC02" : "B2M",
    "DC03": "Col1a1",
    "sgDC03": "Col1a1",
    "gDC03": "Col1a1" 
}


