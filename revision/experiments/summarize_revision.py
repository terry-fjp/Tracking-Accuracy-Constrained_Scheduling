"""Complete-matrix aggregation, paired differences, tables and manuscript inserts."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ.setdefault('MPLCONFIGDIR','/tmp/revision-20261018-mpl')
from pathlib import Path
import json,csv,hashlib,sys
import numpy as np
from run_revision import ROOT,OUT,PROTOCOL,write
R=ROOT.parent;P=R/'paper'
sys.path.insert(0,str(ROOT/'vendor/plot'))
fields=['position_rmse_m','velocity_rmse_m_s','compliance_pct','timely_targets','late_targets','unobserved_targets','laser_s','occupied_s','planning_wall_s','planning_cpu_s','evaluation_s','observations','canceled_visits']

def csvwrite(p,rows):
 p.parent.mkdir(parents=True,exist_ok=True)
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,keys);w.writeheader();w.writerows(rows)

def main():
 check=json.loads((ROOT/'checks/full-matrix.json').read_text());assert check['status']=='PASS' and check['full_default_plan_reproductions']==80
 (P/'figures').mkdir(parents=True,exist_ok=True)
 rows=[json.loads(f.read_text()) for f in sorted((OUT/'runs').glob('*.json'))]
 jrows=[json.loads(f.read_text()) for f in sorted((OUT/'j2/runs').glob('*.json'))]
 summary={};seedrows=[]
 for group in [(s['id'],m) for s in PROTOCOL['settings'] for m in s['methods']]+[(mode,m) for mode in ['J2-two-body','J2-J2'] for m in ['IP-AD','IP-F20-H']]:
  setting,method=group;rs=[r for r in rows+jrows if r.get('setting',r.get('mode'))==setting and r['method']==method];assert len(rs)==40
  d=dict(setting=setting,method=method,scenes=8,seeds=5,aggregation='mean and sample SD of five equally weighted eight-scene means')
  for field in fields:
   if not all(field in r for r in rs):continue
   vals=[]
   for seed in PROTOCOL['search_seeds']:
    a=[r[field] for r in rs if r['seed']==seed];assert len(a)==8;v=float(np.mean(a));vals.append(v);seedrows.append(dict(setting=setting,method=method,seed=seed,metric=field,eight_scene_mean=v))
   d[field]=dict(mean=float(np.mean(vals)),sd=float(np.std(vals,ddof=1)),seed_means=vals)
  summary[setting+'/'+method]=d
 flat=[]
 for d in summary.values():
  f={k:v for k,v in d.items() if not isinstance(v,dict)}
  for k,v in d.items():
   if isinstance(v,dict):f[k+'_mean']=v['mean'];f[k+'_sd']=v['sd']
  flat.append(f)
 csvwrite(OUT/'summary.csv',flat);csvwrite(OUT/'seed-means.csv',seedrows)
 csvwrite(OUT/'scene-seed-results.csv',[{k:r[k] for k in ['run_id','scene','seed','method','target_count']+fields if k in r} for r in rows+jrows])
 pairs=[]
 for s in PROTOCOL['settings']:
  if len(s['methods'])!=2:continue
  for sid in PROTOCOL['scenes']:
   for seed in PROTOCOL['search_seeds']:
    a,b=[next(r for r in rows if r['setting']==s['id'] and r['scene']==sid and r['seed']==seed and r['method']==m) for m in s['methods']]
    pairs.append(dict(setting=s['id'],scene=sid,seed=seed,contrast='IP-AD minus IP-F20-H',**{f:a[f]-b[f] for f in fields if f in a and f in b}))
 for d in summary.values():
  setting,method=d['setting'],d['method']
  if setting=='hold-200':continue
  baseline='population-24' if setting.startswith('population') else 'hold-200'
  for sid in PROTOCOL['scenes']:
   for seed in PROTOCOL['search_seeds']:
    a=next(r for r in rows+jrows if r.get('setting',r.get('mode'))==setting and r['scene']==sid and r['seed']==seed and r['method']==method)
    b=next(r for r in rows if r['setting']==baseline and r['scene']==sid and r['seed']==seed and r['method']==method)
    pairs.append(dict(setting=setting,scene=sid,seed=seed,contrast=method+' minus '+baseline,**{f:a[f]-b[f] for f in fields if f in a and f in b}))
 for sid in PROTOCOL['scenes']:
  for seed in PROTOCOL['search_seeds']:
   for method in ['IP-AD','IP-F20-H']:
    a,b=[next(r for r in jrows if r['mode']==mode and r['scene']==sid and r['seed']==seed and r['method']==method) for mode in ['J2-J2','J2-two-body']]
    pairs.append(dict(setting='J2-filter-contrast',scene=sid,seed=seed,contrast=method+' J2 estimator minus two-body estimator under J2 truth',**{f:a[f]-b[f] for f in fields if f in a and f in b}))
   for mode in ['J2-two-body','J2-J2']:
    a,b=[next(r for r in jrows if r['mode']==mode and r['scene']==sid and r['seed']==seed and r['method']==m) for m in ['IP-AD','IP-F20-H']]
    pairs.append(dict(setting=mode,scene=sid,seed=seed,contrast='IP-AD minus IP-F20-H',**{f:a[f]-b[f] for f in fields if f in a and f in b}))
 csvwrite(OUT/'paired-differences.csv',pairs);write(OUT/'summary.json',summary)
 def val(s,m,f):return summary[s+'/'+m][f]['mean']
 def fmt(s,m,f,n=2):
  d=summary[s+'/'+m][f];return f"{d['mean']:,.{n}f} ± {d['sd']:,.{n}f}"
 def texval(s,m,f,n=2):return '$'+fmt(s,m,f,n).replace(' ± ',r'\pm ').replace(',','{,}')+'$'
 def series(settings,metric,n=2,method='IP-AD'):return ', '.join(f'{val(s,method,metric):,.{n}f}' for s in settings)
 hs=[f'hold-{n}' for n in [100,200,300,500]];ps=[f'population-{n}' for n in [24,48,72,96]]
 he=['hold-200']+[s['id'] for s in PROTOCOL['settings'] if s['group']=='heuristic']
 comp=[val(s,'IP-AD','compliance_pct') for s in he];ers=[val(s,'IP-AD','position_rmse_m') for s in he]
 findings={
 'HOLD_FINDINGS':f"For H = 100/200/300/500 s, IP-AD's five-seed scene-mean position RMSEs are {series(hs,'position_rmse_m',1)} m, and full-interval compliance is {series(hs,'compliance_pct')}%, respectively. The appendix and source tables report both methods, sample standard deviations, and paired scene–seed differences.",
 'HOLD_FINDINGS_ZH':f"H=100/200/300/500 s 时，IP-AD 五种子场景均值位置 RMSE 依次为 {series(hs,'position_rmse_m',1)} m，全程达标率为 {series(hs,'compliance_pct')}%。附录和源数据表报告两种方法、样本标准差和场景—种子配对差异。",
 'HEURISTIC_FINDINGS':f"The baseline position RMSE is {val('hold-200','IP-AD','position_rmse_m'):,.1f} m with {val('hold-200','IP-AD','compliance_pct'):.2f}% compliance. Across the nine settings, IP-AD's five-seed mean position RMSE ranges from {min(ers):,.1f} to {max(ers):,.1f} m and compliance from {min(comp):.2f}% to {max(comp):.2f}%. Full per-setting results, seed variability, and resource outcomes retain the original baseline. A separately labeled warning variant is assessed on new search seeds.",
 'HEURISTIC_FINDINGS_ZH':f"基线位置 RMSE 为 {val('hold-200','IP-AD','position_rmse_m'):,.1f} m、达标率 {val('hold-200','IP-AD','compliance_pct'):.2f}%。九个设置的 IP-AD 五种子位置 RMSE 均值范围为 {min(ers):,.1f}–{max(ers):,.1f} m，达标率为 {min(comp):.2f}%–{max(comp):.2f}%。完整报告各设置、种子变异及资源指标，保留原基线；另行标注预警变体并使用新搜索种子评价。",
 'WARNING_FINDINGS':f"For warning 0.7/0.8/0.9, IP-AD's mean compliance is {series(['warning-0.7','hold-200','warning-0.9'],'compliance_pct')}%, respectively; the appendix also reports RMSE, variability, and laser use.",
 'WARNING_FINDINGS_ZH':f"预警 0.7/0.8/0.9 下 IP-AD 平均达标率依次为 {series(['warning-0.7','hold-200','warning-0.9'],'compliance_pct')}%，附录同时报告 RMSE、变异和激光占用。",
 'MULTIPLIER_FINDINGS':f"For multipliers 1.25/1.5/1.75, mean compliance is {series(['multiplier-1.25','hold-200','multiplier-1.75'],'compliance_pct')}%, respectively; matched errors and resource use are retained in the appendix.",
 'MULTIPLIER_FINDINGS_ZH':f"倍率 1.25/1.5/1.75 的平均达标率依次为 {series(['multiplier-1.25','hold-200','multiplier-1.75'],'compliance_pct')}%，附录同时保留配对误差和资源指标。",
 'POPULATION_FINDINGS':f"For 24/48/72/96 targets, IP-AD's mean compliance is {series(ps,'compliance_pct')}%, and the mean numbers of unobserved targets are {series(ps,'unobserved_targets')}, respectively. At 96 targets, IP-F20-H has mean compliance {val('population-96','IP-F20-H','compliance_pct'):.2f}%; position RMSE is {val('population-96','IP-AD','position_rmse_m'):,.1f} m for IP-AD and {val('population-96','IP-F20-H','position_rmse_m'):,.1f} m for IP-F20-H. Saved admission keys exclude a mean of 8.075 targets, while 52.500 admitted targets receive no observations, leaving 35.425 observed targets for IP-AD at this load. Admission is permission to attempt service, not a continuous-tracking guarantee. All matched error/resource metrics are reported in the appendix and source tables.",
 'POPULATION_FINDINGS_ZH':f"24/48/72/96 目标时，IP-AD 平均达标率依次为 {series(ps,'compliance_pct')}%，平均未观测目标数为 {series(ps,'unobserved_targets')}。96 目标下 IP-F20-H 达标率为 {val('population-96','IP-F20-H','compliance_pct'):.2f}%，IP-AD 与 IP-F20-H 的位置 RMSE 分别为 {val('population-96','IP-AD','position_rmse_m'):,.1f}、{val('population-96','IP-F20-H','position_rmse_m'):,.1f} m。96 目标 IP-AD 平均有 8.075 个目标被接纳键排除、52.500 个已接纳但未观测，实际观测 35.425 个；允许尝试服务不等于持续跟踪保证。附录与源数据报告全部配对误差/资源指标。"
 }
 geometry={sid:json.loads((OUT/'j2/geometry'/f'{sid}.json').read_text()) for sid in PROTOCOL['scenes']}
 maximum=lambda kind,t,f:max(g[kind][str(t)][f] for g in geometry.values())
 lost=sum(g['visibility']['lost'] for g in geometry.values());gained=sum(g['visibility']['gained'] for g in geometry.values());allvis=sum(g['visibility']['all_samples'] for g in geometry.values())
 findings['J2_GEOMETRY']=f"Across the eight scenes, the maximum deterministic target position differences are {maximum('targets',200,'position_max_m'):.2f} m at 200 s and {maximum('targets',900,'position_max_m'):.2f} m at 900 s; corresponding velocity differences are {maximum('targets',200,'velocity_max_m_s'):.3f} and {maximum('targets',900,'velocity_max_m_s'):.3f} m/s. The visibility grid has {lost:,} lost and {gained:,} gained target–observer–time samples out of {allvis:,} samples. Observer differences and per-scene values are also reported."
 findings['J2_GEOMETRY_ZH']=f"八场景中目标确定性位置差异最大值在 200 s、900 s 分别为 {maximum('targets',200,'position_max_m'):.2f}、{maximum('targets',900,'position_max_m'):.2f} m，速度分别为 {maximum('targets',200,'velocity_max_m_s'):.3f}、{maximum('targets',900,'velocity_max_m_s'):.3f} m/s。共 {allvis:,} 个目标—观测星—时间样本中，可见性损失 {lost:,} 个、新增 {gained:,} 个；另报告观测星差异及逐场景值。"
 masks=[json.loads(f.read_text()) for f in (OUT/'j2/plans').glob('*.json')];assert len(masks)==80
 canceled=sum(len(m['canceled']) for m in masks);original_visits=sum(m['original_visits'] for m in masks)
 modes=['hold-200','J2-two-body','J2-J2']
 findings['J2_FINDINGS']=f"For IP-AD in those three modes, the mean position RMSEs are {series(modes,'position_rmse_m',1)} m and compliance rates are {series(modes,'compliance_pct')}%. Across the 80 original calendars, {canceled} of {original_visits} visits are canceled under nominal J₂ execution checks. The fixed-duration comparison is retained, so the baseline-to-J₂ comparison includes geometry and execution changes as well as dynamics; the two J₂ estimator modes isolate the estimator-model change under the same executed calendar."
 findings['J2_FINDINGS_ZH']=f"IP-AD 三模式平均位置 RMSE 依次为 {series(modes,'position_rmse_m',1)} m，达标率为 {series(modes,'compliance_pct')}%。80 份原日历的 {original_visits} 次观测中有 {canceled} 次在名义 J₂ 执行检查后取消。保留固定时长对照，因此基线至 J₂ 比较同时含几何和执行变化；两种 J₂ 估计模式则在相同执行日历下比较估计器模型变化。"
 write(OUT/'response-findings.json',findings)
 # Appendix tables use a new counter namespace, preserving original Table 9.
 groups=[('Hold duration',[(s,m) for s in hs for m in ['IP-AD','IP-F20-H']]),('Local heuristics',[(s,'IP-AD') for s in he]),('Nested populations',[(s,m) for s in ps for m in ['IP-AD','IP-F20-H']]),('Fixed-calendar dynamics',[(s,m) for s in modes for m in ['IP-AD','IP-F20-H']])]
 def short(s):return s.replace('hold-','H=').replace('population-','N=').replace('warning-','w=').replace('multiplier-','m=').replace('weights-','W')
 tables=[]
 for idx,(title,comb) in enumerate(groups,1):
  tab=r'\begin{table}[htbp]\centering\small\setlength{\tabcolsep}{3pt}'+'\n'+r'\begin{revision}\caption{'+title+r'. Means $\pm$ sample standard deviations of five equally weighted eight-scene means. Compliance uses 500 m/25 m/s throughout 120--900 s.}\label{tab:supp-'+str(idx)+r'}\end{revision}'+'\n'+r'\begin{revision}\begin{tabular}{llrrr}\toprule Setting & Method & Position (m) & Velocity (m/s) & Compliance (\%)\\\midrule'+'\n'
  for setting,m in comb:tab+=' & '.join([short(setting),m,texval(setting,m,'position_rmse_m',1),texval(setting,m,'velocity_rmse_m_s',2),texval(setting,m,'compliance_pct',2)])+r'\\'+'\n'
  tab+=r'\bottomrule\end{tabular}\end{revision}\end{table}'+'\n';tables.append(tab)
 resources=r'\begin{table}[p]\centering\small\setlength{\tabcolsep}{4pt}\begin{revision}\caption{Service and computational outcomes. Values are five-seed means; full seed standard deviations, occupied time, wall time, evaluation time, and scene--seed paired differences are retained in the source tables. Planning CPU excludes preprocessing and EKF replay. $T/L/U$ count timely, late, and unobserved targets per scene. Laser time is summed over observers. A dash denotes no scheduling optimization in fixed-calendar replay.}\label{tab:supp-resource}\end{revision}\begin{revision}\begin{tabular}{llrrrrr}\toprule Setting & Method & $T$ & $L$ & $U$ & Laser (s) & CPU (s)\\\midrule'+'\n'
 for _,comb in groups:
  for setting,m in comb:
   values=[f'{val(setting,m,f):.2f}' for f in ['timely_targets','late_targets','unobserved_targets','laser_s']]
   values+=[f'{val(setting,m,"planning_cpu_s"):.2f}' if 'planning_cpu_s' in summary[setting+'/'+m] else '--']
   resources+=' & '.join([short(setting),m]+values)+r'\\'+'\n'
 resources+=r'\bottomrule\end{tabular}\end{revision}\end{table}'
 appendix=r'''\begin{revision}\section{Supplementary Sensitivity Experiments}\label{supplementary-sensitivity}\end{revision}
\setcounter{table}{0}\renewcommand{\thetable}{C\arabic{table}}\renewcommand{\theHtable}{C.\arabic{table}}
\setcounter{figure}{0}\renewcommand{\thefigure}{C\arabic{figure}}\renewcommand{\theHfigure}{C.\arabic{figure}}
\begin{revision}
These experiments retain reference P1: planning and actual-error thresholds 500 m/25 m/s, noise 100 $\mu$rad/50 m, step 0.2 s, horizon 0--900 s, and evaluation 120--900 s. Each scheduling run uses 48 complete-plan evaluations and search seeds 3101--3105; every target has 200 nonlinear realizations. Means and sample standard deviations use five equal eight-scene means, not 200 independent scene replicates. All targets, failed visits, and poor outcomes remain included. The archived default plans and errors are reused only after exact plan reproduction for all 80 default runs and independent nonlinear spot checks. Default planning times are remeasured in this revision. The original baseline is retained; a post-review warning variant is assessed separately in Section~\ref{supp-warning-validation}.

The hold study replans IP-AD and IP-F20-H separately at $H=100,200,300,500$ s. The heuristic study varies one factor at a time in IP-AD: warning 0.7/0.8/0.9; duration multiplier 1.25/1.5/1.75; weights (0.5,0.3,0.2), (0.6,0.2,0.2), (0.4,0.4,0.2), (0.6,0.3,0.1), and (0.4,0.3,0.3). The shared baseline gives nine distinct heuristic settings. Earlier weight calibration used different scores; these new comparisons do not establish the original parameter-selection rationale.

For population scaling, each seed 2101--2108 generates a 96-target mother scene before taking prefixes of 24/48/72/96. Shared targets preserve initial states, priors, and standardized random inputs; observers, resources, and noise remain fixed within each family. New identifiers 5101--5108 distinguish these families from the original scenes. The same total constellation has 72 observers, but local geometric service availability is much smaller. This is a workload study within synthetic orbital families, not a representative catalog or a capacity bound. The 96-target case is an overload boundary, not a nominal service requirement. Admission keys allow attempts to serve a subset, but rejection or lack of observations does not remove a target from the evaluation denominator. At this boundary, IP-AD has higher compliance but worse position RMSE than the fixed-duration control (Table~\ref{tab:supp-3}); all results are retained. IP-AD admits a mean of 87.925 targets but observes only 35.425; 8.075 targets are excluded by admission and 52.500 are admitted without receiving an observation. Admission therefore does not certify sufficient service.

The $J_2$ experiment applies $J_2=1.08262668\times10^{-3}$ to targets and observers with a fixed inertial $z$ pole, the same $\mu$ and $R_E$, and RK4 integration. Its perturbing acceleration is
\[
\boldsymbol a_{J_2}=\frac{3J_2\mu R_E^2}{2r^5}\begin{bmatrix}x(5z^2/r^2-1)\\y(5z^2/r^2-1)\\z(5z^2/r^2-3)\end{bmatrix}.
\]
The analytic acceleration Jacobian updates the matched estimator transition; process noise remains unchanged. The potential-gradient, finite-difference Jacobian, zero-$J_2$, and step-halving checks are saved with the reproducibility materials. Observer ephemerides are treated as exact. Revalidation uses nominal perturbed trajectories, not realization-dependent acquisition or return detection. Both propagations start from the same Cartesian states obtained from the original circular-orbit elements. The original $H=200$ s calendars are frozen, including setup times; a new slew must fit the reserved setup interval. Entire visits failing acquisition, geometric visibility, shared-field, or transition checks are canceled without repair. Both estimators under $J_2$ truth use the same execution mask; changing only the estimator dynamics then isolates that comparison. The unperturbed mask check retains all visits in all 80 original calendars. This does not optimize new $J_2$ schedules or validate atmospheric drag.

Physical conditions share random initial deviations, process-noise draws, and standardized measurement-noise draws. Random draws occur on every grid step regardless of observations. Different dynamics share random inputs but need not share trajectory hashes. Cache keys include physical inputs, actual measurement calendars, dynamics/evaluator versions, protocol, and realization count. No incomplete run enters these summaries. Four processes use one BLAS thread each; planning CPU and wall times are recorded separately. Replay timings include cache reuse and are not standalone estimator benchmarks.
\end{revision}
'''
 appendix=appendix.replace('\n\n','\n\n\\end{revision}\n\\begin{revision}\n\n')
 for t in tables:appendix+=t+'\n'
 appendix+='\\clearpage\n'+resources+'\n\\clearpage\n\\begin{revision}\n'+findings['J2_GEOMETRY'].replace('J₂',r'$J_2$').replace('–','--')+'\n\\end{revision}\n'
 geomtab=r'\begin{table}[htbp]\centering\small\begin{revision}\caption{Maximum deterministic differences from two-body propagation. Columns report maxima over all eight scenes; visibility changes are reported in the text.}\label{tab:j2-geometry}\end{revision}\begin{revision}\begin{tabular}{llrr}\toprule Bodies & Time & Position (m) & Velocity (m/s)\\\midrule'+'\n'
 for kind in ['targets','observers']:
  for t in [200,900]:geomtab+=f'{kind.capitalize()} & {t} s & {maximum(kind,t,"position_max_m"):.2f} & {maximum(kind,t,"velocity_max_m_s"):.3f}'+r'\\'+'\n'
  geomtab+=f'{kind.capitalize()} & Full horizon & {max(g[kind]["full_horizon"]["position_max_m"] for g in geometry.values()):.2f} & {max(g[kind]["full_horizon"]["velocity_max_m_s"] for g in geometry.values()):.3f}'+r'\\'+'\n'
 appendix+=geomtab+r'\bottomrule\end{tabular}\end{revision}\end{table}'+'\n\\clearpage\n'
 # Reproducible figures; uncertainty bars are search-seed SD, not confidence intervals.
 import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
 for title,comb in groups:
  settings=list(dict.fromkeys(s for s,m in comb));xvals=[int(s.split('-')[-1]) for s in settings] if title in ['Hold duration','Nested populations'] else list(range(len(settings)));methods=list(dict.fromkeys(m for s,m in comb));fig,axes=plt.subplots(2,2,figsize=(8,5.4),layout='constrained')
  for ax,f,label in zip(axes.flat,['position_rmse_m','velocity_rmse_m_s','compliance_pct','laser_s'],['Position RMSE (m)','Velocity RMSE (m/s)','Full-interval compliance (%)','Laser occupancy (s)']):
   for m in methods:
    y=[summary[s+'/'+m][f]['mean'] for s in settings];err=[summary[s+'/'+m][f]['sd'] for s in settings]
    ax.errorbar(xvals,y,yerr=err,marker='o',linestyle='-' if title in ['Hold duration','Nested populations'] else 'none',capsize=2,label=m)
   ax.set_xticks(xvals,['Ref' if title=='Local heuristics' and s=='hold-200' else short(s) for s in settings],rotation=35 if len(settings)>4 else 0,ha='right' if len(settings)>4 else 'center');ax.set_ylabel(label);ax.grid(alpha=.2)
  axes[0,0].legend(frameon=False);file='supp-'+title.lower().replace(' ','-');fig.savefig(P/'figures'/f'{file}.pdf');fig.savefig(OUT/f'{file}.png',dpi=170);plt.close(fig)
  appendix+=r'\begin{figure}[p]\begin{adjustwidth}{-\extralength}{}\centering\includegraphics[width=.96\fulllength]{figures/'+file+r'.pdf}\end{adjustwidth}\begin{revision}\caption{'+title+(r'. Ref is the default setting; $w$ is the warning threshold, $m$ the duration multiplier, and Wabc denotes weights $(a/10,b/10,c/10)$' if title=='Local heuristics' else '')+r'. Points show five-seed means and bars show sample standard deviations of eight-scene means. All targets remain included; error bars are not confidence intervals.}\end{revision}\end{figure}'+'\n'
 appendix+='\\FloatBarrier\n\\input{parameter-validation-appendix.tex}\n';(P/'supplementary-appendix.tex').write_text(appendix)
 hold_errors=[val(s,'IP-AD','position_rmse_m') for s in hs]
 hold_pass=[val(s,'IP-AD','compliance_pct') for s in hs]
 paragraph='\n\n'.join([
  f"The supplementary experiments retain reference P1. With independent replanning at H = 100, 200, 300, and 500 s, IP-AD's mean position RMSE ranges from {min(hold_errors):,.1f} to {max(hold_errors):,.1f} m and full-interval compliance from {min(hold_pass):.2f}% to {max(hold_pass):.2f}%. Changing a visit-level hold requirement therefore needs an actual-error evaluation of the resulting calendar, rather than an inference from covariance feasibility alone.",
  f"Across nine local heuristic settings, mean IP-AD position RMSE ranges from {min(ers):,.1f} to {max(ers):,.1f} m and compliance from {min(comp):.2f}% to {max(comp):.2f}%. These one-factor comparisons describe local sensitivity; they do not establish the original selection rationale. The original baseline is retained while a selected warning variant is assessed separately.",
  "The nested 24–96-target study exposes overload under fixed local observing resources. Admission keys permit attempts to serve a subset, but all targets remain in the accuracy denominator. The 96-target boundary, including the opposing RMSE and compliance rankings of the two methods, is retained in Appendix C; it does not establish a constellation capacity limit.",
  f"For the original fixed calendars, nominal J₂ execution checks cancel {canceled} of {original_visits} visits. With J₂ truth, IP-AD mean position RMSE is {val('J2-two-body','IP-AD','position_rmse_m'):,.1f} m using the two-body estimator and {val('J2-J2','IP-AD','position_rmse_m'):,.1f} m using the J₂ estimator, compared with {val('hold-200','IP-AD','position_rmse_m'):,.1f} m in the two-body control. Full-interval compliance is {val('J2-two-body','IP-AD','compliance_pct'):.2f}% with the mismatched two-body estimator and {val('J2-J2','IP-AD','compliance_pct'):.2f}% with the matched estimator, compared with {val('hold-200','IP-AD','compliance_pct'):.2f}% in the original control. The two estimators under J₂ truth share the execution mask; comparison with the original control also includes geometry and execution changes. No J₂ scheduling optimization is performed."
 ])
 paragraph=paragraph.replace('%',r'\%').replace('J₂',r'$J_2$').replace('–','--').replace('IP-AD\'s','IP-AD\'s')
 text=r'\begin{revision}\subsection{Hold, Heuristic, Dynamics, and Workload Sensitivity}\label{supplementary-main}\end{revision}'+'\n\\begin{revision}\n'+paragraph+'\n\n'+r'These supplementary comparisons use reference P1 and retain unfavorable outcomes. Appendix~\ref{supplementary-sensitivity} gives complete protocols, both methods where applicable, uncertainty, and service/resource outcomes. The tests do not establish universal robustness, a historically optimized parameter choice, or operational validity of the sensing assumptions.'+'\n\\end{revision}\n'
 (P/'supplementary-main.tex').write_text(text+'\\input{parameter-validation-main.tex}\n')
 print('SUMMARY COMPLETE',len(rows),len(jrows),len(summary),flush=True)
if __name__=='__main__':main()
