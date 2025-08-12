graph [
  directed 1
  node [
    id 0
    label "sucr_e"
    LabelGraphics [
      fontSize 8
    ]
    graphics [
      fill "#FFA500"
    ]
  ]
  node [
    id 1
    label "glc__D_e"
    LabelGraphics [
      fontSize 8
    ]
    graphics [
      fill "#ffff99"
    ]
  ]
  node [
    id 2
    label "glc__D_c"
    LabelGraphics [
      fontSize 8
    ]
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 3
    label "g6p_c"
    LabelGraphics [
      fontSize 8
    ]
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 4
    label "6pgl_c"
    LabelGraphics [
      fontSize 8
    ]
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 5
    label "6pgc_c"
    LabelGraphics [
      fontSize 8
    ]
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 6
    label "co2_c"
    LabelGraphics [
      fontSize 8
    ]
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 7
    label "hco3_c"
    LabelGraphics [
      fontSize 8
    ]
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 8
    label "malcoa_c"
    LabelGraphics [
      fontSize 8
    ]
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 9
    label "coa_c"
    LabelGraphics [
      fontSize 8
    ]
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 10
    label "pmtcoa_c"
    LabelGraphics [
      fontSize 8
    ]
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 11
    label "pmtcoa_r"
    LabelGraphics [
      fontSize 8
    ]
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 12
    label "pmtalc_r"
    LabelGraphics [
      fontSize 8
    ]
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 13
    label "pmtalc_e"
    LabelGraphics [
      fontSize 8
    ]
    graphics [
      fill "#90EE90"
    ]
  ]
  edge [
    source 0
    target 1
    id "SUCR_e"
    label "SUCR_e&#10;(0.3329)"
    rounded_flux 0.3329
    flux "0.3329249999999999 0.332925"
    width 0.332925
    weight 11.935311280597027
  ]
  edge [
    source 1
    target 2
    id "GLCt_c_e"
    label "GLCt_c_e&#10;(0.3329)"
    rounded_flux 0.3329
    flux "0.3329249999999999 0.332925"
    width 0.332925
    weight 11.935311280597027
  ]
  edge [
    source 2
    target 3
    id "HEX1_c"
    label "HEX1_c&#10;(0.3329)"
    rounded_flux 0.3329
    flux "0.332925 0.3329249999999999"
    width 0.3329249999999999
    weight 11.93531128059703
  ]
  edge [
    source 3
    target 4
    id "G6PDH2i_c"
    label "G6PDH2i_c&#10;(0.8268 0.6808)"
    rounded_flux "0.8268 0.6808"
    flux "0.8267944705917778 0.680786390787768"
    width 0.680786390787768
    weight 5.836725824520052
  ]
  edge [
    source 4
    target 5
    id "PGL_c"
    label "PGL_c&#10;(0.8268 0.6808)"
    rounded_flux "0.8268 0.6808"
    flux "0.8267944705917778 0.680786390787768"
    width 0.680786390787768
    weight 5.836725824520052
  ]
  edge [
    source 5
    target 6
    id "GND_c"
    label "GND_c&#10;(0.8268 0.6808)"
    rounded_flux "0.8268 0.6808"
    flux "0.8267944705917778 0.680786390787768"
    width 0.680786390787768
    weight 5.836725824520052
  ]
  edge [
    source 6
    target 7
    id "HCO3E_c"
    label "HCO3E_c&#10;(0.8152 0.8522)"
    rounded_flux "0.8152 0.8522"
    flux "0.8151703024022106 0.8521934513928917"
    width 0.8521934513928917
    weight 4.6627482311651915
  ]
  edge [
    source 7
    target 8
    id "ACCOAC_c"
    label "ACCOAC_c&#10;(0.7615 0.8147)"
    rounded_flux "0.7615 0.8147"
    flux "0.76147935404601 0.8147233180595586"
    width 0.8147233180595586
    weight 4.877193790840153
  ]
  edge [
    source 8
    target 9
    id "MCOATA_c"
    label "MCOATA_c&#10;(0.761)"
    rounded_flux 0.761
    flux 0.7610437469499546
    width 0.7610437469499546
    weight 8.704902484142611
    graphics [
      fill "#ff0000"
      arrow "last"
    ]
    LabelGraphics [
      color "#ff0000"
    ]
  ]
  edge [
    source 8
    target 10
    id "lumpFACS160_c"
    label "lumpFACS160_c&#10;(0.8099)"
    rounded_flux 0.8099
    flux 0.8099261163929985
    width 0.8099261163929985
    weight 4.90608146554035
    graphics [
      fill "#0000ff"
      arrow "last"
    ]
    LabelGraphics [
      color "#0000ff"
    ]
  ]
  edge [
    source 9
    target 10
    id "PCOATA160_c"
    label "PCOATA160_c&#10;(0.1069)"
    rounded_flux 0.1069
    flux 0.1069070923353133
    width 0.1069070923353133
    weight 61.96793363893191
    graphics [
      fill "#ff0000"
      arrow "last"
    ]
    LabelGraphics [
      color "#ff0000"
    ]
  ]
  edge [
    source 10
    target 11
    id "PMTCOAt_c_r"
    label "PMTCOAt_c_r&#10;(0.1066 0.1142)"
    rounded_flux "0.1066 0.1142"
    flux "0.1065627880563086 0.1142344227882721"
    width 0.1142344227882721
    weight 34.78429190698123
  ]
  edge [
    source 11
    target 12
    id "FARpmtalc_r"
    label "FARpmtalc_r&#10;(0.1066 0.1142)"
    rounded_flux "0.1066 0.1142"
    flux "0.1065627880563086 0.1141710361216055"
    width 0.1141710361216055
    weight 34.80360381297106
  ]
  edge [
    source 12
    target 13
    id "TPS_pmtalc_r"
    label "TPS_pmtalc_r&#10;(0.1066 0.1142)"
    rounded_flux "0.1066 0.1142"
    flux "0.1065627880563087 0.1141710361216055"
    width 0.1141710361216055
    weight 34.80360381297106
  ]
]
