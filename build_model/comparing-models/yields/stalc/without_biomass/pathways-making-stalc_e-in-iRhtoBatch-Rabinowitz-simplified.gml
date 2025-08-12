graph [
  directed 1
  node [
    id 0
    label "glc__D_e"
    graphics [
      fill "#FFA500"
    ]
  ]
  node [
    id 1
    label "glc__D_c"
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 2
    label "g6p_c"
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 3
    label "6pgl_c"
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 4
    label "6pgc_c"
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 5
    label "co2_c"
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 6
    label "hco3_c"
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 7
    label "malcoa_c"
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 8
    label "coa_c"
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 9
    label "stcoa_c"
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 10
    label "stcoa_rm"
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 11
    label "stalc_rm"
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 12
    label "stalc_e"
    graphics [
      fill "#90EE90"
    ]
  ]
  edge [
    source 0
    target 1
    id "GLCt_c_e"
    label "GLCt_c_e (0.3329)"
    flux 0.3329249999999993
    width 0.3329249999999993
    weight 11.939415577776252
  ]
  edge [
    source 1
    target 2
    id "HEX1_c"
    label "HEX1_c (0.3329)"
    flux 0.3329249999999993
    width 0.3329249999999993
    weight 11.939415577776252
  ]
  edge [
    source 2
    target 3
    id "G6PDH2i_c"
    label "G6PDH2i_c (1.0367)"
    flux 1.036725040892735
    width 1.036725040892735
    weight 3.834121656604624
  ]
  edge [
    source 3
    target 4
    id "PGL_c"
    label "PGL_c (1.0367)"
    flux 1.036725040892735
    width 1.036725040892735
    weight 3.834121656604624
  ]
  edge [
    source 4
    target 5
    id "GND_c"
    label "GND_c (1.0367)"
    flux 1.036725040892735
    width 1.036725040892735
    weight 3.834121656604624
  ]
  edge [
    source 5
    target 6
    id "HCO3E_c"
    label "HCO3E_c (1.0194)"
    flux 1.0194102040476716
    width 1.0194102040476716
    weight 3.899244794144975
  ]
  edge [
    source 6
    target 7
    id "ACCOAC_c"
    label "ACCOAC_c (1.0194)"
    flux 1.0194102040476716
    width 1.0194102040476716
    weight 3.899244794144975
  ]
  edge [
    source 7
    target 8
    id "MCOATA_c"
    label "MCOATA_c (1.0194)"
    flux 1.019410204047672
    width 1.019410204047672
    weight 3.8992447941449737
  ]
  edge [
    source 8
    target 9
    id "PCOATA180_c"
    label "PCOATA180_c (0.1274)"
    flux 0.127426275505959
    width 0.127426275505959
    weight 31.19395835315979
  ]
  edge [
    source 9
    target 10
    id "STCOAt_c_rm"
    label "STCOAt_c_rm (0.1274)"
    flux 0.1274262755059589
    width 0.1274262755059589
    weight 31.193958353159815
  ]
  edge [
    source 10
    target 11
    id "FARstalc_rm"
    label "FARstalc_rm (0.1274)"
    flux 0.1274262755059589
    width 0.1274262755059589
    weight 31.193958353159815
  ]
  edge [
    source 11
    target 12
    id "TPS_stalc_rm"
    label "TPS_stalc_rm (0.1274)"
    flux 0.1274262755059589
    width 0.1274262755059589
    weight 31.193958353159815
  ]
]
