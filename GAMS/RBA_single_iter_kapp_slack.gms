* run RBA without bisection, using default settings and outputs
$setGlobal gms %system.FN%
$include %model_root_path%GAMS/defaults.gms
solve rba using lp minimizing z;
$include %model_root_path%GAMS/default_outputs.gms

* kapp_slack_ub('RXN-12AMANTF_g_FWD-rt2093') - kapp_slack_lb('RXN-12AMANTF_g_FWD-rt2093') + v('ENZLOAD-12AMANTF_g_FWD-rt2093') * kapp('RXN-12AMANTF_g_FWD-rt2093') =e= mu * v('RXN-12AMANTF_g_FWD-rt2093');