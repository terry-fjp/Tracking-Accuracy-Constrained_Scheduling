"""Audit the complete prospective validation, then generate reports and TeX.

This reader does not change the frozen candidate, runner, or experiment outputs.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import sys,json,hashlib,csv
from pathlib import Path
from types import SimpleNamespace
from dataclasses import asdict
import numpy as np
D=Path(__file__).resolve().parent;sys.path.insert(0,str(D.parent))
import run_revision as r
from space_cluster_scheduling.model import Config
from revised_algorithm import revised_scene
R=D.parents[1];P=json.loads((D/'protocol.json').read_text());OUT=D/'results'
FIELDS=['position_rmse_m','velocity_rmse_m_s','compliance_pct','timely_targets','late_targets','unobserved_targets','laser_s','occupied_s','observations','planning_cpu_s','planning_wall_s','preprocessing_s','evaluation_s']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def csvwrite(p,rows):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
 for n,h in P['code_hashes'].items():assert sha(D/n)==h,n
 assert sha(D.parent/'results/summary.json')==P['selection_data_sha256']
 assert r.PROTOCOL['fingerprint']==P['parent_physical_protocol']
 for n,h in r.PROTOCOL['source_hashes'].items():assert sha(D.parent/n)==h,n
 rows=[json.loads(p.read_text()) for p in sorted((OUT/'runs').glob('*.json'))]
 expected={(s,k,v) for s in P['scenes'] for k in P['search_seeds'] for v in P['variants']}
 assert len(rows)==len(expected)==120,len(rows)
 assert {(x['scene'],x['seed'],x['variant']) for x in rows}==expected
 truths={};chains={};plans={}
 for x in rows:
  sid=x['scene'];v=x['variant'];n=x['target_count'];name=x['run_id']
  assert x['protocol']==P['fingerprint'] and x['physical_protocol']==r.PROTOCOL['fingerprint']
  sc=r.scene_for(sid,dict(group='validation',overrides={}))
  if v=='IP-AD-W07':sc=revised_scene(sc)
  stub=SimpleNamespace(c=Config(**sc['config']),scene=sc)
  assert x['configuration']==json.loads(json.dumps(asdict(stub.c)))
  assert n==len(sc['targets'])==x['timely_targets']+x['late_targets']+x['unobserved_targets']
  pd=json.loads((OUT/'plans'/f'{name}.json').read_text());plan=pd['plan'];plans[(sid,x['seed'],v)]=plan
  assert pd['protocol']==P['fingerprint'] and pd['search']['evaluations']==48 and pd['search']['seed']==x['seed']
  assert pd['metrics']['resource_violations']==pd['metrics']['hold_violations']==0
  assert len(plan)==x['observations'] and len(x['target_results'])==len(x['truth_hashes'])==n
  first={i:min(z['end'] for z in plan if z['target']==i) for i in {z['target'] for z in plan}}
  assert x['timely_targets']==sum(t<=600 for t in first.values())
  assert x['late_targets']==sum(t>600 for t in first.values())
  assert x['unobserved_targets']==n-len(first)
  assert np.isclose(sum((z['end']-z['start'])*.2 for z in plan),x['laser_s'])
  rr=[]
  for i,key in enumerate(x['target_results']):
   if key not in chains:chains[key]=json.loads((OUT/'chains'/f'{key}.json').read_text())
   c=chains[key];rr.append(c)
   calendar=[{k:z[k] for k in ['start','end','satellite']} for z in plan if z['target']==i]
   oldsig=r.digest(dict(physical=r.physical_target(stub,i,sid),measurements=calendar))
   assert c['physical_calendar_signature']==oldsig
   assert key==c['signature']==r.digest(dict(physical_calendar_signature=oldsig,validation_protocol=P['fingerprint']))
   assert c['validation_protocol']==P['fingerprint'] and c['measurement_schedule']==calendar
   assert c['repetitions']==200 and c['truth_sha256']==x['truth_hashes'][i]
   assert c['full_interval_pass_count']==int(np.sum((np.array(c['max_position_error_per_realization_m'])<=500)&(np.array(c['max_velocity_error_per_realization_m_s'])<=25)))
   if (sid,i) in truths:assert truths[sid,i]==c['truth_sha256']
   truths[sid,i]=c['truth_sha256']
  for field,target in [('position_rmse_m','time_position_rmse_m'),('velocity_rmse_m_s','time_velocity_rmse_m_s')]:
   assert np.isclose(x[field],np.sqrt(np.mean([c[target]**2 for c in rr])),rtol=1e-13)
  assert abs(x['compliance_pct']-100*sum(c['full_interval_pass_count'] for c in rr)/(n*200))<1e-10
  assert all(np.isfinite(x[k]) for k in FIELDS)
 # Verify the explicitly shared MC streams against the old search-seed baseline.
 for sid in P['scenes']:
  old=json.loads((r.OUT/'runs'/f'hold-200-{sid}-s3101-IP-AD.json').read_text())
  assert all(truths[sid,i]==h for i,h in enumerate(old['truth_hashes']))
 changed=sum(plans[s,k,'IP-AD']!=plans[s,k,'IP-AD-W07'] for s in P['scenes'] for k in P['search_seeds'])
 assert changed>0
 summary={};seedrows=[]
 for v in P['variants']:
  summary[v]={}
  for f in FIELDS:
   means=[]
   for k in P['search_seeds']:
    values=[x[f] for x in rows if x['variant']==v and x['seed']==k];assert len(values)==8
    mean=float(np.mean(values));means.append(mean);seedrows.append(dict(variant=v,seed=k,metric=f,eight_scene_mean=mean))
   summary[v][f]=dict(mean=float(np.mean(means)),sd=float(np.std(means,ddof=1)),seed_means=means)
 lookup={(x['scene'],x['seed'],x['variant']):x for x in rows};pairs=[];contrasts={}
 for a,b in [('IP-AD-W07','IP-AD'),('IP-AD-W07','IP-F20-H'),('IP-AD','IP-F20-H')]:
  these=[]
  for s in P['scenes']:
   for k in P['search_seeds']:
    d=dict(contrast=a+' minus '+b,scene=s,seed=k,**{f:lookup[s,k,a][f]-lookup[s,k,b][f] for f in FIELDS});these.append(d)
  pairs+=these;contrast={}
  for f in FIELDS:
   means=[float(np.mean([x[f] for x in these if x['seed']==k])) for k in P['search_seeds']]
   vals=[x[f] for x in these];contrast[f]=dict(mean=float(np.mean(means)),sd=float(np.std(means,ddof=1)),negative_pairs=sum(z<0 for z in vals),positive_pairs=sum(z>0 for z in vals),equal_pairs=sum(z==0 for z in vals))
  contrasts[a+' minus '+b]=contrast
 r.write(OUT/'summary.json',summary);r.write(OUT/'contrasts.json',contrasts)
 csvwrite(OUT/'summary.csv',[dict(variant=v,metric=f,mean=d['mean'],sample_sd=d['sd']) for v,fs in summary.items() for f,d in fs.items()])
 csvwrite(OUT/'seed-means.csv',seedrows);csvwrite(OUT/'paired-differences.csv',pairs)
 csvwrite(OUT/'scene-seed-results.csv',[{k:x[k] for k in ['run_id','scene','seed','variant','target_count']+FIELDS} for x in rows])
 check=dict(status='PASS',runs=len(rows),unique_chains=len(chains),protocol=P['fingerprint'],code_hashes_match=True,selection_data_unchanged=True,all_targets_included=True,metric_reaggregation=True,calendar_and_cache_signatures=True,shared_truth_hashes=True,matching_old_truth_streams=True,warning_changes_calendar_pairs=changed,pairs_per_contrast=40,aggregation='mean and sample SD of five equally weighted eight-scene means',summary_sha256=sha(OUT/'summary.json'))
 r.write(D/'checks.json',check)
 def m(v,f):return summary[v][f]['mean']
 def value(v,f,d=2):
  x=summary[v][f];return f"{x['mean']:,.{d}f} ± {x['sd']:,.{d}f}"
 a='IP-AD-W07';b='IP-AD';c='IP-F20-H'
 en=(f"The local sensitivity results were used to select one new variant, IP-AD-W07, with warning 0.7; hold 200 s, multiplier 1.5, and weights (0.5, 0.3, 0.2) were retained. We then froze this choice before assessing all three methods on new search seeds 3106–3110. Their mean position RMSEs were {m(a,'position_rmse_m'):,.1f} m (IP-AD-W07), {m(b,'position_rmse_m'):,.1f} m (original IP-AD), and {m(c,'position_rmse_m'):,.1f} m (IP-F20-H); corresponding compliance rates were {m(a,'compliance_pct'):.2f}%, {m(b,'compliance_pct'):.2f}%, and {m(c,'compliance_pct'):.2f}%. The appendix reports all 120 runs through aggregate tables and archived scene–seed results, including both favorable and unfavorable paired differences. These are new optimizer seeds on the same eight geometries and shared Monte Carlo inputs, not independent operational or geometric validation. Original 0.8 results remain the reference baseline, and the new variant is explicitly a post-review selection.")
 zh=(f"依据本轮局部敏感性结果，仅选择预警 0.7 的新变体 IP-AD-W07；保持 200 s、倍率 1.5、权重 (0.5,0.3,0.2) 不变。选择冻结后使用新搜索种子 3106–3110 比较三方法。IP-AD-W07、原 IP-AD、IP-F20-H 平均位置 RMSE 分别为 {m(a,'position_rmse_m'):,.1f}、{m(b,'position_rmse_m'):,.1f}、{m(c,'position_rmse_m'):,.1f} m，达标率分别为 {m(a,'compliance_pct'):.2f}%、{m(b,'compliance_pct'):.2f}%、{m(c,'compliance_pct'):.2f}%。附录统计及归档场景—种子结果保留全部 120 次运行和有利/不利配对差异。这只是相同八场景、共享蒙特卡洛输入下的新优化搜索种子验证，不是独立几何或运行验证。原 0.8 结果继续作为基线，新变体明确属于修回阶段选择。")
 r.write(OUT/'response-findings.json',dict(VALIDATION_FINDINGS=en,VALIDATION_FINDINGS_ZH=zh))
 def tex(s):return s.replace('%',r'\%').replace('±',r'$\pm$').replace('–','--')
 intro=r'''\begin{revision}\subsection{New Search Seeds for a Selected Warning Variant}\label{supp-warning-validation}
The nine-setting sensitivity study was used as selection data. Its warning 0.7 setting reduced mean position and velocity RMSE and increased compliance at similar laser use relative to warning 0.8. Only this warning setting was selected for IP-AD-W07; the hold remains 200 s, the multiplier 1.5, and the weights (0.5,0.3,0.2). The candidate and reporting protocol were frozen before evaluating search seeds 3106--3110. The original IP-AD and IP-F20-H controls were independently replanned with these same seeds, yielding 120 runs at 48 complete-plan evaluations each and 200 realizations per target. All eight original scenes and Monte Carlo inputs are shared with the selection study. Thus this check assesses variation across new optimization seeds; it does not validate transfer to unseen geometries or independent stochastic truth samples. All outcomes are retained. The original warning 0.8 remains the reference for previously reported experiments.

\end{revision}
'''
 tab=r'\begin{table}[htbp]\centering\small\setlength{\tabcolsep}{3pt}\begin{revision}\caption{Selected warning variant and matched controls on search seeds 3106--3110. Means $\pm$ sample standard deviations of five eight-scene means; the acceptance criterion remains 500 m/25 m/s throughout 120--900 s.}\label{tab:warning-validation}\end{revision}\begin{revision}\begin{tabular}{lrrr}\toprule Method & Position (m) & Velocity (m/s) & Compliance (\%)\\\midrule'+'\n'
 for v in [b,c,a]:tab+=' & '.join([v,tex(value(v,'position_rmse_m',1)),tex(value(v,'velocity_rmse_m_s')),tex(value(v,'compliance_pct'))])+r'\\'+'\n'
 tab+=r'\bottomrule\end{tabular}\end{revision}\end{table}'+'\n'
 tab+=r'\begin{table}[htbp]\centering\small\setlength{\tabcolsep}{3pt}\begin{revision}\caption{Service and resource outcomes of the new-seed comparison. Means over all 40 scene--seed pairs; source tables also report sample standard deviations. $T/L/U$ denote timely, late, and unobserved targets.}\label{tab:warning-resources}\end{revision}\begin{revision}\begin{tabular}{lrrrrr}\toprule Method & $T$ & $L$ & $U$ & Laser (s) & CPU (s)\\\midrule'+'\n'
 for v in [b,c,a]:tab+=' & '.join([v]+[f'{m(v,f):.2f}' for f in ['timely_targets','late_targets','unobserved_targets','laser_s','planning_cpu_s']])+r'\\'+'\n'
 tab+=r'\bottomrule\end{tabular}\end{revision}\end{table}'+'\n'
 d=contrasts[a+' minus '+b]
 pairtext=f"Against original IP-AD, IP-AD-W07 has lower position RMSE in {d['position_rmse_m']['negative_pairs']} of 40 pairs and higher compliance in {d['compliance_pct']['positive_pairs']} of 40 pairs. Mean paired differences are {d['position_rmse_m']['mean']:,.1f} m in position RMSE and {d['compliance_pct']['mean']:.2f} percentage points in compliance. Pair counts are descriptive and are not independent statistical tests across geometries."
 (R/'paper/parameter-validation-appendix.tex').write_text(intro+tab+'\n\\par\n\\begin{revision}\n'+tex(pairtext)+'\n\\end{revision}\n\\FloatBarrier\n')
 maintext=f"A warning-0.7 variant, selected from the local sensitivity study, was assessed separately on new search seeds 3106–3110. IP-AD-W07, original IP-AD, and IP-F20-H attained mean position RMSEs of {m(a,'position_rmse_m'):,.1f}, {m(b,'position_rmse_m'):,.1f}, and {m(c,'position_rmse_m'):,.1f} m and compliance of {m(a,'compliance_pct'):.2f}%, {m(b,'compliance_pct'):.2f}%, and {m(c,'compliance_pct'):.2f}%, respectively. The original baseline is retained; this post-review comparison shares the original geometries and Monte Carlo inputs."
 (R/'paper/parameter-validation-main.tex').write_text('\\begin{revision}\n'+tex(maintext)+r' Appendix~\ref{supp-warning-validation} reports the protocol and complete comparisons.'+'\n\\end{revision}\n')
 report=['# 补充实验对结论的影响与参数更新','本报告供作者内部核查，所有达标率保留全部目标与 200 次实现。','## 参数决定',zh,'200 s 是跟踪要求，不能为了增益而改成更宽松的评价任务。1.5 与当前权重保持不变。本轮不联合寻找最优参数，不重定义原基线。','## 新搜索种子对照','| 方法 | 位置 RMSE (m) | 速度 RMSE (m/s) | 全程达标率 (%) | 激光占用 (s) |','|---|---|---|---|---|']
 for v in [b,c,a]:report.append('| '+' | '.join([v,value(v,'position_rmse_m',1),value(v,'velocity_rmse_m_s'),value(v,'compliance_pct'),value(v,'laser_s')])+' |')
 report+=['均值 ± 样本标准差按五个搜索种子的八场景均值计算；不是 200 个独立场景的标准误。',pairtext,'## 对核心结论的影响','原文关于自适应时长的条件性优势继续使用原对照结果；新变体只能支持本轮指定条件下的改善，不支持所有指标、场景或噪声条件上的普遍优势。全程达标率仍远低于 100%，不能声称完成全目标持续跟踪任务。','保持扫描的收益非单调。H=100/200/300/500 s 时原 IP-AD 达标率为 30.73/31.43/33.36/29.40%，较长的单次保持可行性不自动改善完整日历的实际误差。','J₂ 真值配两体估计器的全程达标率为 0.00%，但平均 RMSE 变化较小。匹配 J₂ 估计器后达标率为 31.41%。这加强了“协方差可行性不等于实际误差保证”的结论，也说明模型一致性对持续验收特别重要；不是阻力、全部摄动或飞行验证。','## 96 目标应如何呈现','算法包含目标接纳键，允许仅尝试服务部分目标，实际安排还受几何与资源可行性限制。这不是已验证的最优目标子集选择器。晚初始化、未观测、拒绝服务目标都保留在总体分母。','保存计划显示，96 目标 IP-AD 平均接纳 87.925 个、观测 35.425 个；未观测可分为接纳键排除 8.075 个及已接纳但未排到观测 52.500 个。平均未观测 60.58 个，达标率 4.52%，固定对照为 1.44%；IP-AD 的位置 RMSE 86,605.8 m 劣于固定对照 81,214.3 m。服务不足明显，但尚不能分离几何、资源预算、搜索预算、启发式策略等原因，因此不能推导星座绝对容量。','按用户最新方向压缩主文，将 96 目标作为过载边界，完整数据及两种方法的相反指标排序保留附录和第三位审稿人回复。没有删除原始结果，也不新增 96 目标实验。若未来改成任务级准入控制，应同时报告全体覆盖率及已接纳目标表现，不能仅更换达标率分母。','## 复现与验证','冻结协议见 parameter-update/protocol.json；检查见 parameter-update/checks.json；全部 120 个运行、40 对/对照的配对差异和五种子统计均保留。新搜索种子不等于新几何或新蒙特卡洛真值流。']
 report.insert(report.index('## 对核心结论的影响'),f"新变体相对原 IP-AD 的位置 RMSE 降低 {(1-m(a,'position_rmse_m')/m(b,'position_rmse_m'))*100:.2f}%，速度 RMSE 降低 {(1-m(a,'velocity_rmse_m_s')/m(b,'velocity_rmse_m_s'))*100:.2f}%。平均总占用从 {m(b,'occupied_s'):.2f} s 增至 {m(a,'occupied_s'):.2f} s，规划 CPU 从 {m(b,'planning_cpu_s'):.2f} s 增至 {m(a,'planning_cpu_s'):.2f} s，平均未观测目标从 {m(b,'unobserved_targets'):.2f} 增至 {m(a,'unobserved_targets'):.2f}。激光占用仅为量测时间，不能代替包括机动和捕获等开销的总占用。")
 (R/'audit/conclusion-impact-and-parameter-update.md').write_text('\n\n'.join(report).replace(' |\n\n|', ' |\n|').replace('|\n\n|', '|\n|')+'\n')
 print(json.dumps(check,indent=2));print(en)
if __name__=='__main__':main()
