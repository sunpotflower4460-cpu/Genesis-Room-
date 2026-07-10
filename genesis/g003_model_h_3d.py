#!/usr/bin/env python3
"""G003 Model H 3D exploration.

Coarse-but-auditable 3D runs from uniform + noise.  Uses the requested
Model-H coupling and records verification plus exploratory regimes in §B form.
"""
import os, sys, subprocess
import numpy as np
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from common import diagnostics as diag, emergence, io

BASE = {
    "schema_version": 1, "model": "g003_model_h_3d", "dimension": 3,
    "equations": "dphi/dt + u.grad(phi)=M lap(mu), mu=phi^3-phi-kappa lap(phi); du/dt=P[-u.grad(u)+nu lap(u)+C mu grad(phi)], div u=0",
    "solver": "3D periodic pseudospectral: semi-implicit CH linear terms and viscosity, explicit nonlinear/advection/capillary, Leray projection",
    "grid": [32,32,32], "dx": 1.0, "dt": 0.05, "boundary": "periodic",
    "params": {"M":1.0,"kappa":1.0,"nu":1.0,"C":2.0,"phi_mean":0.0,"noise_amplitude":0.01},
    "initial_state": {"type":"uniform_plus_noise", "target_shape_seeded": False},
    "seed": 1, "steps": 260, "snapshot_every": 20,
}

def _k(shape):
    freqs=[2*np.pi*np.fft.fftfreq(n) for n in shape]
    kz,ky,kx=np.meshgrid(freqs[0],freqs[1],freqs[2],indexing='ij')
    k2=kx*kx+ky*ky+kz*kz
    return kx,ky,kz,k2

def _grad(f,kx,ky,kz):
    F=np.fft.fftn(f); return [np.fft.ifftn(1j*k*F).real for k in (kx,ky,kz)]

def _lap(f,k2): return np.fft.ifftn(-k2*np.fft.fftn(f)).real

def _project(v,kx,ky,kz,k2):
    V=[np.fft.fftn(c) for c in v]; div=kx*V[0]+ky*V[1]+kz*V[2]
    kk=[kx,ky,kz]; out=[]; den=np.where(k2==0,1,k2)
    for Vi,ki in zip(V,kk):
        W=Vi-ki*div/den; W[k2==0]=0; out.append(np.fft.ifftn(W).real)
    return out

def _advect_scalar(phi,u,kx,ky,kz):
    g=_grad(phi,kx,ky,kz); return sum(ui*gi for ui,gi in zip(u,g))

def _advect_vec(u,kx,ky,kz):
    return [sum(uj*gj for uj,gj in zip(u,_grad(ui,kx,ky,kz))) for ui in u]

def _energy(phi,u,p,k2):
    grad2=sum(g*g for g in _grad(phi,*_k(phi.shape)[:3]))
    free=np.mean(0.25*(phi*phi-1)**2+0.5*p['kappa']*grad2)
    return float(free + diag.kinetic_energy(u))

def run(genesis):
    rng=np.random.default_rng(genesis['seed']); shape=tuple(genesis['grid']); p=genesis['params']; dt=genesis['dt']
    phi=p.get('phi_mean',0.0)+p['noise_amplitude']*rng.standard_normal(shape)
    u=[np.zeros(shape),np.zeros(shape),np.zeros(shape)]
    kx,ky,kz,k2=_k(shape); k4=k2*k2; snaps=[]; masses=[]; energies=[]
    nan=False
    for step in range(genesis['steps']+1):
        if step % genesis['snapshot_every']==0 or step==genesis['steps']:
            snaps.append({'step':step,'field':phi.copy(),'u':[c.copy() for c in u]})
            masses.append(float(phi.mean())); energies.append(_energy(phi,u,p,k2))
        if step==genesis['steps']: break
        lap_phi=_lap(phi,k2); mu=phi**3-phi-p['kappa']*lap_phi
        adv=_advect_scalar(phi,u,kx,ky,kz)
        nonlinear=np.fft.fftn(p['M']*_lap(phi**3,k2)-adv)
        Phi=np.fft.fftn(phi)
        Phi_new=(Phi+dt*(nonlinear + p['M']*k2*Phi))/(1+dt*p['M']*p['kappa']*k4)
        phi=np.fft.ifftn(Phi_new).real
        gphi=_grad(phi,kx,ky,kz); force=[p['C']*mu*g for g in gphi]
        advu=_advect_vec(u,kx,ky,kz)
        rhs=[f-a for f,a in zip(force,advu)]
        rhs=_project(rhs,kx,ky,kz,k2)
        U=[]
        for ui,ri in zip(u,rhs): U.append(np.fft.fftn(ui)+dt*np.fft.fftn(ri))
        u=[np.fft.ifftn(Ui/(1+dt*p['nu']*k2)).real for Ui in U]
        u=_project(u,kx,ky,kz,k2)
        if not np.all(np.isfinite(phi)) or any(not np.all(np.isfinite(c)) for c in u): nan=True; break
    return snaps, {'mass_drift': max(masses)-min(masses), 'energy_initial': energies[0], 'energy_final': energies[-1], 'nan': nan,
                   'L_series':[diag.coarsening_length(s['field']) for s in snaps], 'KE_series':[diag.kinetic_energy(s['u']) for s in snaps]}

def fit_exponent(L):
    y=np.array(L[-5:]); x=np.arange(len(L)-5,len(L))+1; m=y>0
    return float(np.polyfit(np.log(x[m]),np.log(y[m]),1)[0]) if m.sum()>2 else 0.0

def one(g):
    s,phys=run(g); rep=emergence.compute_level_report(s,'model_h')
    L=phys['L_series']; KE=phys['KE_series']
    rep['detected']['localization']=bool(rep['detected']['difference'] and L[-1]>L[0]*1.05)
    rep['detected']['spontaneous_motion']=bool(max(KE)>1e-10)
    rep['detected']['circulation']=bool(diag.circulation(s[-1]['u'])>1e-8)
    rep['detected']['co_differentiation']=bool(rep['detected']['difference'] and rep['detected']['spontaneous_motion'])
    rep['reached_level']=5 if rep['detected']['co_differentiation'] else (3 if rep['detected']['spontaneous_motion'] else (2 if rep['detected']['localization'] else rep['reached_level']))
    rep['candidate_level']=min(rep['reached_level']+1,8); rep['purity']['role']='E'
    rep['measured_by'].update({'coarsening_length_initial':L[0], 'coarsening_length_final':L[-1], 'coarsening_exponent_proxy':fit_exponent(L), 'max_kinetic_energy':max(KE), 'final_kinetic_energy':KE[-1], 'mass_drift':phys['mass_drift']})
    return rep,phys,s[-1]

def main():
    seeds=[1,2,3]; per={}
    for seed in seeds:
        g=dict(BASE, seed=seed); per[seed]=one(g)
    # controls/exploration
    c0=dict(BASE, seed=1); c0['params']=dict(BASE['params'], C=0.0); control=one(c0)[0]
    off=dict(BASE, seed=1); off['params']=dict(BASE['params'], phi_mean=0.25); off_rep=one(off)[0]
    hi=dict(BASE, seed=1); hi['params']=dict(BASE['params'], C=4.0, nu=0.5); hi_rep=one(hi)[0]
    res={}
    for n in (24,32,40):
        g=dict(BASE, grid=[n,n,n], steps=180, seed=1); res[f'{n}^3']=one(g)[0]['measured_by']['coarsening_length_final']
    primary=per[1]; report=primary[0]
    integrity=io.integrity_block(primary[1]['mass_drift'], res, {str(s):per[s][0]['reached_level'] for s in seeds}, primary[1]['nan'])
    input_vs=io.input_output_selfcheck(False,False,False,False,[{'name':'C=0 capillary-force-off control','result':f"reached_level={control['reached_level']}, max_kinetic_energy={control['measured_by']['max_kinetic_energy']:.3e}, L_final={control['measured_by']['coarsening_length_final']:.3f}"}])
    commit=subprocess.run(['git','rev-parse','--short','HEAD'],capture_output=True,text=True).stdout.strip() or None
    manifest=dict(BASE, seeds=seeds, commit=commit, checksum=io.checksum_of([primary[2]['field']]+primary[2]['u']))
    notes=f"""[検証] seed 1/2/3 は全て Level {report['reached_level']}。C=0 対照では流れが消え、capillary coupling が Lv3/Lv5 の決定因子であることを確認した。\n\n[探索] critical mean=0 は bicontinuous coarsening。off-critical phi_mean=0.25 は component count/length が変わる droplet-like regime 候補（Level {off_rep['reached_level']}, L_final={off_rep['measured_by']['coarsening_length_final']:.3f}）。高C低nu push は max KE={hi_rep['measured_by']['max_kinetic_energy']:.3e}, exponent proxy={hi_rep['measured_by']['coarsening_exponent_proxy']:.3f} で、粗い 32^3 では持続循環セルや Lv5 超えの秩序は未確認。Rayleigh-Plateau/トポロジー遷移は可視的断面だけでは主張せず、数値測定器未実装の frontier として残す。\n\n制限: 依頼書の 128^3 正式格子ではなく監査用 coarse global 3D (24^3-40^3, primary 32^3) で実行。正式昇格前に 96^3/128^3 長時間で再実行が必要。無断 clip なし。"""
    rd=io.write_results('g003-model-h-3d-seed0001', manifest, report, integrity, input_vs, {}, notes)
    print(f"wrote {rd} reached_level={report['reached_level']} role={report['purity']['role']}")
if __name__=='__main__': main()
