"""Confirm the J2 replay adapter changes dynamics only, with original random streams."""
from run_revision import *
import j2_revision as j
from space_cluster_scheduling.run_accuracy import chain
from space_cluster_scheduling.nonlinear_math import rk4,transition
p=j.Geometry(scene_for(2101,PROTOCOL['settings'][1]));p.x=p.two_x;p.obs=p.two_obs
plan=json.loads((ARCHIVE/'plans/P1-2101-s3101-b48-IP-AD.json').read_text())['plan'];visits=[v for v in plan if v['target']==0]
j._ns.update(truth_step=rk4,estimate_step=rk4,transition=transition)
a=chain(p,0,visits,2101,200);b=j._ns['j2_chain'](p,0,visits,2101,200)
for k in ['truth_sha256','time_position_rmse_m','time_velocity_rmse_m_s','full_interval_pass_count','seed']:assert a[k]==b[k],k
for k in ['max_position_error_per_realization_m','max_velocity_error_per_realization_m_s']:assert np.array_equal(a[k],b[k]),k
write(ROOT/'checks/j2-replay-adapter.json',dict(status='PASS',two_body_adapter_exact=True,truth_hash=a['truth_sha256'],seed=a['seed'],repetitions=200,standardized_noise_draws='Original evaluator retained verbatim; draws occur at all steps, independent of actual observation calendar.',scope='Original evaluator versus dynamics-parameterized adapter with both dynamics set to original two-body functions; scene2101,target0,IP-AD,seed3101.'))
print('J2 REPLAY ADAPTER PASS')
