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
    weight 7.856847603705226
  ]
  edge [
    source 1
    target 2
    id "HEX1_c"
    label "HEX1_c (0.3329)"
    flux 0.3329249999999999
    width 0.3329249999999999
    weight 7.856847603705229
  ]
  edge [
    source 2
    target 3
    id "G6PDH2i_c"
    label "G6PDH2i_c (0.5342)"
    flux 0.5341967288935505
    width 0.5341967288935505
    weight 4.896587431153668
  ]
  edge [
    source 3
    target 4
    id "PGL_c"
    label "PGL_c (0.5342)"
    flux 0.5341967288935505
    width 0.5341967288935505
    weight 4.896587431153668
  ]
  edge [
    source 4
    target 5
    id "GND_c"
    label "GND_c (0.5342)"
    flux 0.5341967288935505
    width 0.5341967288935505
    weight 4.896587431153668
  ]
  edge [
    source 5
    target 6
    id "HCO3E_c"
    label "HCO3E_c (1.0292)"
    flux 1.0292385988463546
    width 1.0292385988463546
    weight 2.5414330471044084
  ]
  edge [
    source 6
    target 7
    id "ACCOAC_c"
    label "ACCOAC_c (1.0292)"
    flux 1.0292385988463548
    width 1.0292385988463548
    weight 2.541433047104408
  ]
  edge [
    source 7
    target 8
    id "lumpFACS180_c"
    label "lumpFACS180_c (1.0292)"
    flux 1.0292385988463544
    width 1.0292385988463544
    weight 2.541433047104409
  ]
  edge [
    source 8
    target 9
    id "STCOAt_c_rm"
    label "STCOAt_c_rm (0.1287)"
    flux 0.1286548248557943
    width 0.1286548248557943
    weight 20.33146437683527
  ]
  edge [
    source 9
    target 10
    id "FARstalc_rm"
    label "FARstalc_rm (0.1287)"
    flux 0.1286548248557943
    width 0.1286548248557943
    weight 20.33146437683527
  ]
  edge [
    source 10
    target 11
    id "TPS_stalc_rm"
    label "TPS_stalc_rm (0.1287)"
    flux 0.1286548248557943
    width 0.1286548248557943
    weight 20.33146437683527
  ]
]
