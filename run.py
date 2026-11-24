"""Reproduce every result: python run.py [--size 1.0] [--out results]

--size scales the k-meshes (1.0 = full run; 0.125 = pilot).
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).parent / ".mplcache"))
sys.path.insert(0, str(Path(__file__).parent / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from tbgraphene import (Graphene, Tube, dos_analytic_nn, dos_numeric, fermi_velocity,
                        find_dirac_points, min_direct_gap, tubes_in_range, HBAR_EV_S)

T_HOP = 2.7            # eV, user-specified
TP_REP = 0.1 * T_HOP   # representative t' (Castro Neto 2009 range 0.02t..0.2t) -- assumption
DELTAS = [0.05, 0.10, 0.25]  # illustrative sublattice masses, eV -- assumption
TUBES = [(5, 5), (9, 0), (10, 0), (8, 4), (13, 0)]
METAL_TOL = 1e-5       # eV: numerical gap below this counts as zero


def hexagon(gr):
    K = gr.high_symmetry()["K"]
    r = np.linalg.norm(K)
    ang = np.arctan2(K[1], K[0]) + np.arange(7) * np.pi / 3
    return r * np.c_[np.cos(ang), np.sin(ang)]


def part_graphene(gr, size, out, S):
    hs = gr.high_symmetry()
    # --- Dirac points ----------------------------------------------------
    dps = find_dirac_points(gr)
    S["dirac_points"] = [
        dict(k_inv_angstrom=k.tolist(), E_c_minus_E_v_eV=float(np.diff(gr.bands(k))[0]),
             dist_to_analytic_K_or_Kp=float(min(np.linalg.norm(k - hs["K"]), np.linalg.norm(k - hs["Kp"]))))
        for k in dps]
    # --- Fermi velocity ----------------------------------------------------
    vf = {}
    for i, k in enumerate(dps):
        v, vs, ED = fermi_velocity(gr, k, qmax=2e-3, ndir=72)
        vf[f"dirac_{i}"] = dict(v_mean_m_s=float(v), v_min=float(vs.min()), v_max=float(vs.max()),
                                rel_anisotropy=float((vs.max() - vs.min()) / v), E_D_eV=float(ED))
    vref = 3 * gr.t * gr.acc * 1e-10 / (2 * HBAR_EV_S)
    S["fermi_velocity"] = dict(per_dirac_point=vf, analytic_3t_acc_over_2hbar_m_s=vref,
                               hbar_vF_eV_angstrom=float(vf["dirac_0"]["v_mean_m_s"] * HBAR_EV_S / 1e-10),
                               fit="E_c-E_v symmetric half, quadratic fit, |q|<=2e-3 1/A, 72 directions")
    # warping: velocity along K->Gamma vs K->M at larger q
    K = dps[0]
    warp = {}
    for q in [0.02, 0.1, 0.3]:
        e = {}
        for name, tgt in (("toward_Gamma", hs["G"]), ("toward_M", hs["M"])):
            # use the image of the target closest to this K
            u = (tgt - hs["K"]) / np.linalg.norm(tgt - hs["K"])
            e[name] = float(gr.bands(hs["K"] + q * u)[1])
        warp[f"q={q}"] = e
    S["trigonal_warping_Ec_eV"] = warp

    # --- Band surfaces -----------------------------------------------------
    n = max(41, int(301 * size))
    kx = np.linspace(-1.25 * np.linalg.norm(hs["K"]), 1.25 * np.linalg.norm(hs["K"]), n)
    KX, KY = np.meshgrid(kx, kx)
    E = gr.bands(np.stack([KX, KY], -1))
    fig = plt.figure(figsize=(13, 5.5))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    for b, cm in ((0, "Blues_r"), (1, "Reds")):
        ax.plot_surface(KX, KY, E[..., b], cmap=cm, rstride=2, cstride=2, linewidth=0, alpha=0.9)
    ax.set_xlabel(r"$k_x$ (1/Å)"); ax.set_ylabel(r"$k_y$ (1/Å)"); ax.set_zlabel("E (eV)")
    ax.set_title(f"Graphene π bands, NN TB, t = {gr.t} eV")
    ax.view_init(18, -60)
    ax2 = fig.add_subplot(1, 2, 2)
    cs = ax2.contourf(KX, KY, E[..., 1], 30, cmap="viridis")
    fig.colorbar(cs, ax=ax2, label=r"$E_c$ (eV)")
    h = hexagon(gr)
    ax2.plot(h[:, 0], h[:, 1], "w-", lw=1.2, label="1st BZ")
    for k in dps:
        ax2.plot(*k, "r*", ms=12)
    ax2.plot([], [], "r*", label="Dirac points (numerical)")
    ax2.set_aspect("equal"); ax2.set_xlabel(r"$k_x$ (1/Å)"); ax2.set_ylabel(r"$k_y$ (1/Å)")
    ax2.legend(loc="upper right", fontsize=8); ax2.set_title("Conduction band, full BZ")
    fig.tight_layout(); fig.savefig(out / "fig_band_surface.png", dpi=150); plt.close(fig)

    # --- Path G-K-M-G --------------------------------------------------------
    kp, dist, ticks = gr.path(n_per_seg=max(50, int(400 * size)))
    Ep = gr.bands(kp)
    Ep_num = gr.bands_numeric(kp)
    S["path_closed_form_vs_diag_max_diff_eV"] = float(np.abs(Ep - Ep_num).max())
    S["path_energies_eV"] = {l: gr.bands(hs[l]).tolist() for l in ("G", "K", "M")}
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(dist, Ep, "k-")
    for x in ticks:
        ax.axvline(x, color="0.7", lw=0.8)
    ax.axhline(0, color="0.7", ls=":")
    ax.set_xticks(ticks, [r"$\Gamma$", "K", "M", r"$\Gamma$"])
    ax.set_xlim(dist[0], dist[-1]); ax.set_ylabel("E (eV)")
    ax.set_title(f"Graphene, NN TB, t = {gr.t} eV")
    fig.tight_layout(); fig.savefig(out / "fig_bands_GKMG.png", dpi=150); plt.close(fig)

    # --- DOS -----------------------------------------------------------------
    nk = max(200, int(3000 * size))
    edges = np.linspace(-3.05 * gr.t, 3.05 * gr.t, 1648)
    dos = dos_numeric(gr, nk, edges)
    Ec = 0.5 * (edges[1:] + edges[:-1])
    ana = dos_analytic_nn(Ec, gr.t)
    away = (np.abs(np.abs(Ec) - gr.t) > 0.1) & (np.abs(Ec) < 2.9 * gr.t) & (np.abs(Ec) > 0.3)
    S["dos"] = dict(nk=nk, bin_eV=float(edges[1] - edges[0]),
                    integral=float((dos * np.diff(edges)).sum()),
                    max_rel_dev_from_analytic_away_from_vHs=float(np.max(np.abs(dos[away] - ana[away]) / ana[away])),
                    vhs_peak_energies_eV=[float(Ec[Ec < 0][np.argmax(dos[Ec < 0])]),
                                          float(Ec[Ec > 0][np.argmax(dos[Ec > 0])])],
                    units="states / eV / unit cell / spin (2 atoms per cell)")
    low = (np.abs(Ec) < 0.5) & (np.abs(Ec) > 0.05)
    slope = np.polyfit(np.abs(Ec[low]), dos[low], 1)[0]
    vF = S["fermi_velocity"]["per_dirac_point"]["dirac_0"]["v_mean_m_s"]
    hv = vF * HBAR_EV_S / 1e-10
    S["dos"]["low_E_slope_numeric_per_eV2"] = float(slope)
    S["dos"]["low_E_slope_from_vF_Ac_over_pi_hbar2_vF2"] = float(gr.cell_area / (np.pi * hv ** 2))
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    ax[0].plot(Ec, dos, "C0-", lw=1, label=f"k-mesh {nk}×{nk} histogram")
    ax[0].plot(Ec, ana, "k--", lw=0.8, label="analytic (Hobson–Nierenberg)")
    ax[0].set_xlabel("E (eV)"); ax[0].set_ylabel("DOS (states/eV/cell/spin)")
    ax[0].set_ylim(0, 0.6); ax[0].legend(fontsize=8); ax[0].set_title("Graphene DOS")
    z = np.abs(Ec) < 1.0
    ax[1].plot(Ec[z], dos[z], "C0.", ms=3, label="numeric")
    ax[1].plot(Ec[z], gr.cell_area * np.abs(Ec[z]) / (np.pi * hv ** 2), "r-",
               label=r"$A_c|E|/(\pi\hbar^2v_F^2)$ with fitted $v_F$")
    ax[1].set_xlabel("E (eV)"); ax[1].legend(fontsize=8); ax[1].set_title("Low energy: linear DOS")
    fig.tight_layout(); fig.savefig(out / "fig_dos.png", dpi=150); plt.close(fig)
    np.savetxt(out / "dos.csv", np.c_[Ec, dos, ana], delimiter=",",
               header="E_eV,dos_numeric,dos_analytic (states/eV/cell/spin)", comments="")

    # --- Dirac cone zoom -----------------------------------------------------
    q = np.linspace(-0.3, 0.3, 121)
    QX, QY = np.meshgrid(q, q)
    Ez = gr.bands(K + np.stack([QX, QY], -1))
    fig = plt.figure(figsize=(11, 4.5))
    a = fig.add_subplot(1, 2, 1, projection="3d")
    a.plot_surface(QX, QY, Ez[..., 0], cmap="Blues_r", lw=0)
    a.plot_surface(QX, QY, Ez[..., 1], cmap="Reds", lw=0)
    a.set_xlabel(r"$q_x$ (1/Å)"); a.set_ylabel(r"$q_y$ (1/Å)"); a.set_zlabel("E (eV)")
    a.set_title("Dirac cone at K")
    b = fig.add_subplot(1, 2, 2)
    qq = np.linspace(0, 0.6, 200)
    for name, tgt in (("K→Γ", hs["G"]), ("K→M", hs["M"])):
        u = (tgt - hs["K"]) / np.linalg.norm(tgt - hs["K"])
        b.plot(qq, gr.bands(hs["K"] + qq[:, None] * u)[:, 1], label=name)
    b.plot(qq, hv * qq, "k--", label=r"$\hbar v_F q$")
    b.set_xlabel("|q| (1/Å)"); b.set_ylabel(r"$E_c$ (eV)"); b.legend()
    b.set_title(f"$v_F$ = {vF/1e5:.4f}×10$^5$ m/s (trigonal warping at large q)")
    fig.tight_layout(); fig.savefig(out / "fig_dirac_cone.png", dpi=150); plt.close(fig)


def part_tubes(gr, size, out, S):
    nk = max(201, int(2001 * size))
    rows = []
    fig, axs = plt.subplots(1, len(TUBES), figsize=(16, 5), sharey=True)
    for ax, (n, m) in zip(axs, TUBES):
        tb = Tube(n, m, gr)
        g, info = tb.gap_cutting(nk=nk)
        tb.build_cell()
        kap = np.linspace(-np.pi / tb.Tlen, np.pi / tb.Tlen, max(41, int(201 * size)))
        e_cut = tb.bands_cutting(kap)
        e_rs = tb.bands_realspace(kap)
        # real-space route gap on the same grid (independent of cutting-line code)
        nat = e_rs.shape[1]
        g_rs = tb.gap_realspace(nk=max(101, int(801 * size)))
        metallic = abs(g) < METAL_TOL
        rows.append(dict(n=n, m=m, d_nm=tb.d / 10, theta_deg=tb.theta, N_hex=tb.N, atoms=nat,
                         T_angstrom=tb.Tlen, n_minus_m_mod3=(n - m) % 3,
                         K_on_cutting_line=bool(tb.k_on_cutting_line()),
                         metallic=bool(metallic), gap_eV=0.0 if metallic else g, gap_raw_eV=g,
                         gap_realspace_eV=g_rs,
                         realspace_vs_cutting_max_diff_eV=float(np.abs(e_rs - e_cut).max()),
                         gap_times_d_eV_nm=g * tb.d / 10,
                         ref_2acc_t_over_d_eV=2 * gr.acc * gr.t / tb.d,
                         k_cbm_over_pi_T=info["k_c"] * tb.Tlen / np.pi))
        ax.plot(kap * tb.Tlen / np.pi, e_cut, "k-", lw=0.6)
        ax.plot(kap * tb.Tlen / np.pi, e_rs, "r:", lw=0.6)
        ax.set_ylim(-3, 3); ax.axhline(0, color="0.6", ls=":")
        ax.set_title(f"({n},{m}) d={tb.d/10:.3f} nm\n" + ("metallic" if metallic else f"gap {g:.3f} eV"))
        ax.set_xlabel(r"$k\,T/\pi$")
    axs[0].set_ylabel("E (eV)")
    axs[0].plot([], [], "k-", label="cutting lines"); axs[0].plot([], [], "r:", label="real-space cell")
    axs[0].legend(fontsize=7, loc="lower left")
    fig.tight_layout(); fig.savefig(out / "fig_tube_bands.png", dpi=150); plt.close(fig)
    S["tubes"] = rows
    with open(out / "tubes.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    # --- sweep: gap vs diameter ---------------------------------------------
    sweep = []
    for (n, m) in tubes_in_range(0.5, 2.5, gr):
        tb = Tube(n, m, gr)
        g, _ = tb.gap_cutting(nk=max(101, int(801 * size)))
        sweep.append((n, m, tb.d / 10, g, (n - m) % 3, tb.k_on_cutting_line()))
    sw = np.array([(s[2], s[3], s[4]) for s in sweep])
    with open(out / "gap_sweep.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["n", "m", "d_nm", "gap_eV", "n_minus_m_mod3", "K_on_cutting_line"])
        w.writerows(sweep)
    semi = np.abs(sw[:, 1]) > METAL_TOL
    gd = sw[semi, 0] * sw[semi, 1]
    S["gap_sweep"] = dict(n_tubes=len(sweep), n_semiconducting=int(semi.sum()),
                          metallic_iff_n_minus_m_mod3_eq_0=bool(np.all((sw[:, 2] == 0) == ~semi)),
                          metallic_iff_K_on_line=bool(np.all(np.array([s[5] for s in sweep]) == ~semi)),
                          gap_times_d_eV_nm_mean=float(gd.mean()), gap_times_d_min=float(gd.min()),
                          gap_times_d_max=float(gd.max()),
                          gap_times_d_mean_d_gt_1p5nm=float(gd[sw[semi, 0] > 1.5].mean()),
                          ref_2acc_t_eV_nm=2 * gr.acc * gr.t / 10)
    fig, ax = plt.subplots(figsize=(7, 5))
    mod = sw[:, 2]
    for r, c, lab in ((1, "C0", "(n−m) mod 3 = 1"), (2, "C1", "(n−m) mod 3 = 2"), (0, "C2", "metallic, (n−m) mod 3 = 0")):
        s = mod == r
        ax.plot(sw[s, 0], sw[s, 1], "o", color=c, ms=3.5, alpha=0.6, label=lab)
    dd = np.linspace(0.45, 2.55, 200)
    ax.plot(dd, 2 * gr.acc * gr.t / (10 * dd), "k--", lw=1, label=r"$2a_{cc}t/d$")
    for r in rows:
        ax.plot(r["d_nm"], r["gap_eV"], "r*", ms=13, mec="k")
        ax.annotate(f"({r['n']},{r['m']})", (r["d_nm"], r["gap_eV"]), textcoords="offset points",
                    xytext=(6, 4) if r["gap_eV"] > 0 else ((-30, 10) if r["m"] == r["n"] else (6, 10)),
                    fontsize=9)
    ax.plot([], [], "r*", ms=10, mec="k", label="requested tubes")
    ax.set_xlabel("diameter d (nm)"); ax.set_ylabel("band gap (eV)")
    ax.set_title(f"Zone-folded SWCNT gaps, NN TB t = {gr.t} eV")
    ax.legend(fontsize=8); ax.set_xlim(0.45, 2.55); ax.set_ylim(-0.05, 1.3)
    fig.tight_layout(); fig.savefig(out / "fig_gap_vs_diameter.png", dpi=150); plt.close(fig)


def part_tp_delta(size, out, S):
    res = {}
    base = Graphene(t=T_HOP)
    hs = base.high_symmetry()
    kp, dist, ticks = base.path(n_per_seg=max(50, int(400 * size)))
    variants = [("NN only", Graphene(t=T_HOP)),
                (f"t′={TP_REP:.2f} eV", Graphene(t=T_HOP, tp=TP_REP)),
                (f"t′={2*TP_REP:.2f} eV", Graphene(t=T_HOP, tp=2 * TP_REP))] + \
               [(f"Δ={d:.2f} eV", Graphene(t=T_HOP, delta=d)) for d in DELTAS] + \
               [(f"t′={TP_REP:.2f}, Δ=0.10 eV", Graphene(t=T_HOP, tp=TP_REP, delta=0.10))]
    nk = max(200, int(1500 * size))
    edges = np.linspace(-10.5, 10.5, 1401)
    Ec = 0.5 * (edges[1:] + edges[:-1])
    fig, axs = plt.subplots(1, 3, figsize=(16, 5))
    for i, (lab, gr) in enumerate(variants):
        g, kmins = min_direct_gap(gr)
        kmin = kmins[0]
        e_at = gr.bands(kmin)
        # band edges along everything: global extrema on a mesh
        mesh = gr.bands(gr.kmesh(max(60, int(400 * size))))
        r = dict(t=gr.t, tp=gr.tp, delta=gr.delta, min_direct_gap_eV=g,
                 k_at_min=kmin.tolist(),
                 dist_kmin_to_K=float(min(np.linalg.norm(kmin - hs["K"]), np.linalg.norm(kmin - hs["Kp"]))),
                 E_at_min=e_at.tolist(), band_bottom_eV=float(mesh[:, 0].min()),
                 band_top_eV=float(mesh[:, 1].max()),
                 E_Gamma=gr.bands(hs["G"]).tolist(), E_M=gr.bands(hs["M"]).tolist())
        if gr.delta == 0:
            v, vs, ED = fermi_velocity(gr, kmin)
            r["v_F_m_s"] = float(v); r["E_D_eV"] = float(ED)
            # e-h asymmetry at q = 0.2 1/A toward Gamma
            u = (hs["G"] - hs["K"]) / np.linalg.norm(hs["G"] - hs["K"])
            e = gr.bands(hs["K"] + 0.2 * u)
            r["eh_asym_at_q0.2_toward_G_eV"] = float((e[1] - ED) - (ED - e[0]))
        else:
            # effective mass from curvature of E_c at K: E_c ~ E0 + hbar^2 q^2 / (2 m*)
            # window q <= 0.02 Delta/(hbar v) so the quadratic regime holds (sceptic finding)
            qmax_m = 0.02 * gr.delta / (1.5 * gr.t * gr.acc)
            qs = np.linspace(qmax_m / 20, qmax_m, 20)
            ec = gr.bands(hs["K"] + qs[:, None] * np.array([1.0, 0.0]))[:, 1]
            c2 = np.polyfit(qs, ec, 2)[0]                     # eV A^2
            me = 9.1093837015e-31
            mstar = (HBAR_EV_S ** 2 / (2 * c2 * 1e-20)) * 1.602176634e-19 / me
            r["m_star_over_me"] = float(mstar)
            hv = 3 * gr.t * gr.acc / 2 * 1e-10
            # analytic curvature at K: hbar^2/m* = (hbar v)^2/Delta - 2 t' (1.5 a_cc)^2
            inv = (1.5 * gr.t * gr.acc) ** 2 / gr.delta - 2 * gr.tp * (1.5 * gr.acc) ** 2  # eV A^2
            r["m_star_check_analytic_over_me"] = float(HBAR_EV_S ** 2 / (inv * 1e-20) * 1.602176634e-19 / me)
        dos = dos_numeric(gr, nk, edges)
        r["dos_at_E_D"] = float(dos[np.argmin(np.abs(Ec - np.mean(e_at)))])
        res[lab] = r
        if i in (0, 1, 2, 4, 6):
            axs[0].plot(dist, gr.bands(kp), color=f"C{i}", lw=1, label=[lab, None])
            axs[2].plot(Ec, dos, color=f"C{i}", lw=0.8, label=lab)
        # zoom near K
        if i in (0, 1, 3, 4, 5, 6):
            sel = (dist > ticks[1] - 0.35) & (dist < ticks[1] + 0.35)
            axs[1].plot(dist[sel] - ticks[1], gr.bands(kp[sel]), color=f"C{i}", lw=1, label=[lab, None])
    axs[0].set_xticks(ticks, [r"$\Gamma$", "K", "M", r"$\Gamma$"])
    for x in ticks:
        axs[0].axvline(x, color="0.8", lw=0.6)
    axs[0].set_ylabel("E (eV)"); axs[0].legend(fontsize=7); axs[0].set_title("Bands with t′ and Δ")
    axs[1].set_xlabel("distance from K along path (1/Å)"); axs[1].set_ylim(-1.5, 2.5)
    axs[1].legend(fontsize=7); axs[1].set_title("Near K: t′ shifts cone up by 3t′, Δ opens 2Δ")
    axs[2].set_xlabel("E (eV)"); axs[2].set_ylabel("DOS (states/eV/cell/spin)")
    axs[2].set_ylim(0, 0.6); axs[2].legend(fontsize=7); axs[2].set_title("DOS")
    fig.tight_layout(); fig.savefig(out / "fig_tp_delta.png", dpi=150); plt.close(fig)
    S["tp_delta"] = res

    # tubes with t' and Delta
    trows = []
    for lab, gr in [("t'=0.27", Graphene(t=T_HOP, tp=TP_REP)), ("Delta=0.10", Graphene(t=T_HOP, delta=0.10)),
                    ("t'=0.27,Delta=0.10", Graphene(t=T_HOP, tp=TP_REP, delta=0.10))]:
        for (n, m) in TUBES:
            tb = Tube(n, m, gr)
            g, info = tb.gap_cutting(nk=max(201, int(2001 * size)))
            trows.append(dict(model=lab, n=n, m=m, d_nm=float(tb.d / 10), gap_eV=float(g),
                              gap_realspace_eV=float(tb.gap_realspace(nk=max(101, int(801 * size))))))
    S["tubes_tp_delta"] = trows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--size", type=float, default=1.0)
    ap.add_argument("--out", default="results")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    S = dict(size=a.size, t_eV=T_HOP, acc_angstrom=Graphene().acc)
    gr = Graphene(t=T_HOP)
    part_graphene(gr, a.size, out, S)
    part_tubes(gr, a.size, out, S)
    part_tp_delta(a.size, out, S)
    S["runtime_s"] = time.time() - t0
    (out / "summary.json").write_text(json.dumps(S, indent=2))
    print(json.dumps({k: S[k] for k in ("fermi_velocity", "dirac_points", "dos")}, indent=1))
    for r in S["tubes"]:
        print(f"({r['n']},{r['m']}) d={r['d_nm']:.4f} nm metallic={r['metallic']} gap={r['gap_eV']:.4f} eV "
              f"rs={r['gap_realspace_eV']:.4f} gap*d={r['gap_times_d_eV_nm']:.4f}")
    print("sweep", S["gap_sweep"])
    for k, v in S["tp_delta"].items():
        print(k, {x: v[x] for x in v if x in ("min_direct_gap_eV", "E_at_min", "v_F_m_s", "m_star_over_me", "eh_asym_at_q0.2_toward_G_eV", "band_bottom_eV", "band_top_eV")})
    for r in S["tubes_tp_delta"]:
        print(r)
    print(f"runtime {S['runtime_s']:.1f} s")


if __name__ == "__main__":
    main()
