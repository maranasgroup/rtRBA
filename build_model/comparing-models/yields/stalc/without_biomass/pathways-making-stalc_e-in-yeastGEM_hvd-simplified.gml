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
    label "stcoa_c"
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 9
    label "stcoa_rm"
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 10
    label "stalc_rm"
    graphics [
      fill "#D3D3D3"
    ]
  ]
  node [
    id 11
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
    flux 0.332925
    width 0.332925
    weight 7.851044494896893
  ]
  edge [
    source 1
    target 2
    id "HEX1_c"
    label "HEX1_c (0.3329)"
    flux 0.3329249999999999
    width 0.3329249999999999
    weight 7.851044494896896
  ]
  edge [
    source 2
    target 3
    id "G6PDH2i_c"
    label "G6PDH2i_c (0.5342)"
    flux 0.5341967288935562
    width 0.5341967288935562
    weight 4.892970786019872
  ]
  edge [
    source 3
    target 4
    id "PGL_c"
    label "PGL_c (0.5342)"
    flux 0.5341967288935562
    width 0.5341967288935562
    weight 4.892970786019872
  ]
  edge [
    source 4
    target 5
    id "GND_c"
    label "GND_c (0.5342)"
    flux 0.5341967288935562
    width 0.5341967288935562
    weight 4.892970786019872
  ]
  edge [
    source 5
    target 6
    id "HCO3E_c"
    label "HCO3E_c (1.0292)"
    flux 1.0292385988463555
    width 1.0292385988463555
    weight 2.5395559313392377
  ]
  edge [
    source 6
    target 7
    id "ACCOAC_c"
    label "ACCOAC_c (1.0292)"
    flux 1.0292385988463555
    width 1.0292385988463555
    weight 2.5395559313392377
  ]
  edge [
    source 7
    target 8
    id "lumpFACS180_c"
    label "lumpFACS180_c (1.0292)"
    flux 1.029238598846356
    width 1.029238598846356
    weight 2.539555931339237
  ]
  edge [
    source 8
    target 9
    id "STCOAt_c_rm"
    label "STCOAt_c_rm (0.1287)"
    flux 0.1286548248557945
    width 0.1286548248557945
    weight 20.316447450713895
  ]
  edge [
    source 9
    target 10
    id "FARstalc_rm"
    label "FARstalc_rm (0.1287)"
    flux 0.1286548248557944
    width 0.1286548248557944
    weight 20.316447450713905
  ]
  edge [
    source 10
    target 11
    id "TPS_stalc_rm"
    label "TPS_stalc_rm (0.1287)"
    flux 0.1286548248557944
    width 0.1286548248557944
    weight 20.316447450713905
  ]
]
