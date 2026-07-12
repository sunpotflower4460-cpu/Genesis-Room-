#!/usr/bin/env python3
"""genesis/boussinesq_rb_3d.py -- L3-ルートA【最有力】: Rayleigh-Bénard対流をfull-3Dで。

無次元Boussinesq方程式（長さスケールH=境界間距離、時間スケールH^2/kappa、温度スケールdelta_T,
自由すべり境界＋固定温度差、で無次元化した標準形。この無次元化では線形安定性解析の結果が
古典的な値 Ra_c = 27*pi^4/4 ~ 657.51 に一致する——直接の検証ターゲットとして使う）：

    du/dt + (u.grad)u = -grad(p) + Pr*lap(u) + Ra*Pr*T*zhat   (Boussinesq近似、浮力はz方向のみ)
    dT/dt + (u.grad)T - w = lap(T)                             (Tは背景勾配(下が高温)からの摂動。
        導出: T_total=T0(z)+T, dT0/dz=-1(下が高温、無次元)、(u.grad)T0=w*dT0/dz=-w を移項すると
        dT/dt = -(u.grad)T + w + lap(T) つまり dT/dt+(u.grad)T-w=lap(T)。上昇流(w>0)は下の暖かい
        流体を運び局所的に温める=正のフィードバック(不安定化)。符号を逆にすると常に安定化し
        Raによらず減衰する——自己テスト(Ra=200で減衰・Ra=1200で成長)でこの符号誤りを発見・
        訂正した記録として残す)
    div(u) = 0

境界: 上下(z=0,H)に自由すべり壁(応力なし)＋固定温度(T=0、背景プロファイルに一致)。
  u,v: dz=0 at walls (自由すべり、接線応力ゼロ) -- 偶関数(cosine)対称
  w:   w=0 at walls (不透過)                    -- 奇関数(sine)対称、境界で自動的に0
  T:   T=0 at walls (固定温度、背景に一致)        -- 奇関数(sine)対称、境界で自動的に0
  p:   偶関数(u,vの式のgrad(p)と整合)

数値法（"domain-doubling"トリック）: z方向を[0,H]から[0,2H)へ鏡映拡張し、周期的FFTで扱う。
偶関数拡張は自動的にNeumann(自由すべり)を、奇関数拡張は自動的にDirichlet=0(不透過/固定温度)
を満たす。拡張後のz波数は k_z = n*pi/H と正確に一致する（標準周期FFT波数 2*pi*n/(2H) = n*pi/H
と一致するため）ので、genesis/solvers.py の k_grid・dealias_mask・leray_project を無改造で
再利用できる——3重周期のModel H(g003)と全く同じ機構がそのまま使える。

from-0: 静止流体(u=v=w=0)＋微小ノイズ(Tのみ、奇対称化済み)から開始。対流セル・回転・循環は
一切初期条件に置かない（第8監査L3）。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from genesis.solvers import dealias_mask, leray_project  # noqa: E402


def k_grid(shape, spacing=1.0):
    """genesis.solvers.k_grid の異方性版(軸ごとに異なるspacingを許す、ここではz方向の物理
    格子間隔dzがx,yのdx,dyと一般に異なるため)。spacingがスカラーなら全軸同じ(元のk_gridと
    完全互換)、タプル/リストなら軸ごとに使う。"""
    if np.isscalar(spacing):
        spacing = [spacing] * len(shape)
    ks = [2 * np.pi * np.fft.fftfreq(n, d=d) for n, d in zip(shape, spacing)]
    kk = np.meshgrid(*ks, indexing="ij")
    k2 = sum(k ** 2 for k in kk)
    return kk, k2

MODEL_ID = "boussinesq_rb_3d"
DEFAULTS = {"Pr": 1.0, "Ra": 700.0}  # Ra_c(自由すべり) = 27*pi^4/4 ~ 657.51


def doubled_shape(nz_phys, ny, nx):
    """z方向の物理格子点数(壁を含む、nz_phys)から鏡映拡張後の全体形状を返す。"""
    return (2 * (nz_phys - 1), ny, nx)


def _mirror_index(n):
    idx = np.arange(n)
    return (-idx) % n


def symmetrize_even(field):
    m = _mirror_index(field.shape[0])
    return 0.5 * (field + field[m])


def symmetrize_odd(field):
    m = _mirror_index(field.shape[0])
    return 0.5 * (field - field[m])


def make_initial(nz_phys, ny, nx, noise_amplitude, rng):
    """静止(u=v=w=0)＋微小ノイズ(Tのみ、奇対称化=壁でT=0を自動満足)。"""
    shape = doubled_shape(nz_phys, ny, nx)
    u = np.zeros(shape)
    v = np.zeros(shape)
    w = np.zeros(shape)
    T = symmetrize_odd(noise_amplitude * rng.standard_normal(shape))
    return u, v, w, T


def _grad_z_via_fft(field_hat, kz_col):
    return np.real(np.fft.ifftn(1j * kz_col * field_hat))


def step(u, v, w, T, dt, p, kk, k2, dealias):
    """Chorin投影法の半陰的擬スペクトルステップ。u,vは偶、w,Tは奇の対称性を保つ。

    重要: 配列の軸順は(z,y,x)(axis0=z, symmetrize_even/oddがaxis0をzとして扱う設計、
    doubled_shapeもこの順)。kk=k_grid(shape,...)は軸ごとの波数を返す(kk[0]=kz,kk[1]=ky,
    kk[2]=kx)。leray_project/発散・移流の計算は「fields[a]がkk[a]と同じ軸方向の速度成分」
    という規約が必要——なのでベクトル場のリストは速度そのものの(u,v,w)順でなく、
    軸に対応する(w,v,u)順で渡す(以前 fields=[u,v,w] としていたのは軸順とのミスマッチで
    Leray射影が発散を正しく計算できないバグだった: 自己テストで運動学的にあり得ない
    水平一様モードのwが非ゼロに成長することから発見・訂正)。"""
    ndim = 3
    fields = [w, v, u]  # 軸順(z,y,x)に対応: fields[0]=w(z方向), fields[1]=v(y), fields[2]=u(x)
    fields_hat_full = [np.fft.fftn(f) for f in fields]
    fields_hat = [fh * dealias for fh in fields_hat_full]
    That_full = np.fft.fftn(T)
    That = That_full * dealias

    grad = [[np.real(np.fft.ifftn(1j * kk[b] * fields_hat[a])) for b in range(ndim)]
            for a in range(ndim)]
    advect = [sum(fields[b] * grad[a][b] for b in range(ndim)) for a in range(ndim)]

    grad_T = [np.real(np.fft.ifftn(1j * kk[a] * That)) for a in range(ndim)]
    advect_T = sum(fields[a] * grad_T[a] for a in range(ndim))

    buoyancy = p["Ra"] * p["Pr"] * T  # w(=fields[0], z方向)方程式のみに加わる
    explicit_fields = [-advect[0] + buoyancy, -advect[1], -advect[2]]
    explicit_T = -advect_T + w  # 背景成層の移流: 上昇流(w>0)は下の暖かい流体を運び局所的に
    # 温める(正のフィードバック=不安定化)。標準的な無次元化(例: Chandrasekhar)の符号に一致。

    explicit_hat = [np.fft.fftn(e) * dealias for e in explicit_fields]
    explicit_T_hat = np.fft.fftn(explicit_T) * dealias

    denom = 1.0 + dt * p["Pr"] * k2
    star_hat = [(fields_hat_full[a] + dt * explicit_hat[a]) / denom for a in range(ndim)]
    denom_T = 1.0 + dt * k2
    That_new = (That_full + dt * explicit_T_hat) / denom_T

    proj_hat = leray_project(star_hat, kk, k2)
    w_new, v_new, u_new = [np.real(np.fft.ifftn(ph)) for ph in proj_hat]
    T_new = np.real(np.fft.ifftn(That_new))

    u_new = symmetrize_even(u_new)
    v_new = symmetrize_even(v_new)
    w_new = symmetrize_odd(w_new)
    T_new = symmetrize_odd(T_new)
    return u_new, v_new, w_new, T_new


def stable_dt(nz_phys, ny, nx, Pr, H=1.0, Lx=1.0, Ly=1.0, safety=0.3, dt_min=0.0005, dt_max=0.02):
    """拡散項は半陰的(無条件安定)なので、これは移流/浮力項が格子を飛び越えないための
    保守的なヒューリスティック上限(格子間隔の2乗スケールで、実際のCFL条件の代理)。"""
    shape = doubled_shape(nz_phys, ny, nx)
    dz = H / (nz_phys - 1)
    spacing = (dz, Ly / ny, Lx / nx)
    _, k2 = k_grid(shape, spacing=spacing)
    lap_max = float(k2.max())
    diff_max = max(Pr, 1.0)
    dt = safety / (diff_max * lap_max) if lap_max > 0 else dt_max
    return float(min(max(dt, dt_min), dt_max))


def run(nz_phys, ny, nx, H, Lx, Ly, t_final, seed, params=None, snapshot_every_frac=1.0 / 40,
        noise_amplitude=0.001, dt=None):
    p = dict(DEFAULTS)
    if params:
        p.update(params)
    shape = doubled_shape(nz_phys, ny, nx)
    dz = H / (nz_phys - 1)
    spacing = (dz, Ly / ny, Lx / nx)
    kk, k2 = k_grid(shape, spacing=spacing)
    dealias = dealias_mask(shape)

    dt = dt if dt is not None else stable_dt(nz_phys, ny, nx, p["Pr"], H=H, Lx=Lx, Ly=Ly)
    steps = int(round(t_final / dt))
    snapshot_every = max(1, int(steps * snapshot_every_frac))

    rng = np.random.default_rng(seed)
    u, v, w, T = make_initial(nz_phys, ny, nx, noise_amplitude, rng)

    snapshots = []
    diverged = False
    for t in range(steps):
        u, v, w, T = step(u, v, w, T, dt, p, kk, k2, dealias)
        if not (np.all(np.isfinite(u)) and np.all(np.isfinite(w)) and np.all(np.isfinite(T))):
            diverged = True
            break
        if t % snapshot_every == 0 or t == steps - 1:
            snapshots.append({"step": t, "u": u.copy(), "v": v.copy(), "w": w.copy(), "T": T.copy()})
    if not snapshots:
        snapshots = [{"step": 0, "u": u.copy(), "v": v.copy(), "w": w.copy(), "T": T.copy()}]
    phys = {"diverged": diverged, "dt_used": dt, "steps": steps, "nz_phys": nz_phys, "H": H}
    return snapshots, phys


def physical_slice(field, nz_phys):
    """鏡映拡張(2*(nz_phys-1)点)から物理領域[0,H]のnz_phys点だけを取り出す。"""
    return field[:nz_phys]


if __name__ == "__main__":
    # 自己テスト: Ra=200(<Ra_c~657.5、減衰するはず) vs Ra=1200(>Ra_c、成長するはず)
    nz, ny, nx = 17, 24, 24
    H, Lx, Ly = 1.0, 4.0, 4.0
    for Ra in (200.0, 1200.0):
        snaps, phys = run(nz, ny, nx, H, Lx, Ly, t_final=3.0, seed=1,
                          params={"Pr": 1.0, "Ra": Ra}, noise_amplitude=0.001)
        w_rms = [float(np.sqrt(np.mean(s["w"][:nz] ** 2))) for s in snaps]
        print("Ra=%.0f: w_rms %s" % (Ra, ["%.2e" % x for x in w_rms[::max(1, len(w_rms) // 8)]]))
