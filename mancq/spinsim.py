"""
spinsim.py - exact quantum-mechanical 1H spin-system simulation (no external deps
beyond numpy). Used to turn GISSMO spin-system matrices (chemical shifts + J
couplings, refined against experimental spectra) into stick spectra at the
acquisition field of the data being analysed, so second-order multiplet
intensities (glutamine/glutamate CH2, citrate AB, sugar ring protons) are correct.

Method: Hamiltonian in the product basis, block-diagonalised by total Mz,
H = sum_i nu_i Iz_i + sum_{i<j} J_ij I_i.I_j (Hz). Observable F- = sum_i I-_i over
1H spins only. Heteronuclear couplings (31P in phosphates) are added as a
weakly coupled extra spin (Iz.Iz term only) that is not observed.
"""
import numpy as np
from itertools import combinations


def _components(n, couplings, shifts_hz, eq_tol=0.05):
    """Connected components of the coupling graph. Couplings between spins with
    identical shift are ignored for linking (they do not change the spectrum of
    magnetically equivalent groups such as methyl protons)."""
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for (i, j), J in couplings.items():
        if abs(J) < 0.05:
            continue
        if abs(shifts_hz[i] - shifts_hz[j]) < eq_tol:
            continue
        parent[find(i)] = find(j)
    comps = {}
    for i in range(n):
        comps.setdefault(find(i), []).append(i)
    return list(comps.values())


def simulate_component(shifts_hz, couplings, hetero=None, min_int=1e-4):
    """shifts_hz: list of 1H frequencies (Hz). couplings: dict (i,j)->J with local
    indices. hetero: dict local_i -> J to one heteronucleus (weak coupling).
    Returns (freqs_hz, intensities) with intensities summing to n_spins."""
    n = len(shifts_hz)
    nh = 1 if hetero else 0
    N = n + nh
    nu = np.array(list(shifts_hz) + [0.0] * nh)
    Jm = np.zeros((N, N))
    for (i, j), J in couplings.items():
        Jm[i, j] = Jm[j, i] = J
    if hetero:
        for i, J in hetero.items():
            Jm[i, n] = Jm[n, i] = J
    # states as bit tuples: bit=1 -> alpha (m=+1/2) ; group by number of 1H alphas
    # (heteronucleus m is conserved separately since only zz coupling to it)
    freqs, ints, spins = [], [], []
    het_states = [0, 1] if nh else [None]
    for hs in het_states:
        blocks = {}
        for k in range(n + 1):
            states = []
            for comb in combinations(range(n), k):
                s = np.zeros(n, dtype=np.int8)
                s[list(comb)] = 1
                states.append(s)
            blocks[k] = np.array(states)
        eig = {}
        for k, S in blocks.items():
            m = S - 0.5  # m_i values
            dim = len(S)
            H = np.zeros((dim, dim))
            H[np.diag_indices(dim)] = m @ nu[:n] + 0.5 * np.einsum('ai,ij,aj->a', m, Jm[:n, :n], m) \
                - 0.0
            # remove double counting of i==j (Jm diagonal is zero) -> fine
            if hs is not None:
                mh = hs - 0.5
                H[np.diag_indices(dim)] += m @ Jm[:n, n] * mh
            index = {tuple(r): a for a, r in enumerate(S)}
            for a, r in enumerate(S):
                for i in range(n):
                    for j in range(i + 1, n):
                        if Jm[i, j] != 0 and r[i] != r[j]:
                            r2 = r.copy()
                            r2[i], r2[j] = r[j], r[i]
                            b = index[tuple(r2)]
                            H[a, b] += 0.5 * Jm[i, j]
            w, v = np.linalg.eigh(H)
            eig[k] = (w, v, S, index)
        # transitions k -> k-1 via F- (lowers one alpha to beta)
        for k in range(1, n + 1):
            wu, vu, Su, _ = eig[k]
            wl, vl, Sl, idxl = eig[k - 1]
            Fi = np.zeros((n, len(Sl), len(Su)))
            for a, r in enumerate(Su):
                for i in range(n):
                    if r[i] == 1:
                        r2 = r.copy()
                        r2[i] = 0
                        Fi[i, idxl[tuple(r2)], a] = 1.0
            Ti = np.einsum('lb,ibc,cu->ilu', vl.T, Fi, vu, optimize=True)
            T = Ti.sum(0)
            I = T ** 2
            f = wu[None, :] - wl[:, None]
            sel = I > min_int
            share = Ti * T[None, :, :]
            owner = np.argmax(share, axis=0)
            freqs.extend(f[sel].tolist())
            ints.extend(I[sel].tolist())
            spins.extend(owner[sel].tolist())
    freqs = np.array(freqs)
    ints = np.array(ints)
    if len(ints):
        ints = ints * (n / ints.sum())
    return freqs, ints, np.array(spins, int)


def merge_sticks(ppm, inten, group=None, tol_ppm=0.00005):
    """Merge coincident transitions (degenerate lines) of the same proton group."""
    if group is None:
        group = np.zeros(len(ppm), int)
    if len(ppm) == 0:
        return ppm, inten, group
    o = np.lexsort((ppm, group))
    ppm, inten, group = ppm[o], inten[o], group[o]
    out_p, out_i, out_g = [ppm[0]], [inten[0]], [group[0]]
    for p, i, g in zip(ppm[1:], inten[1:], group[1:]):
        if g == out_g[-1] and p - out_p[-1] < tol_ppm:
            tot = out_i[-1] + i
            out_p[-1] = (out_p[-1] * out_i[-1] + p * i) / tot
            out_i[-1] = tot
        else:
            out_p.append(p); out_i.append(i); out_g.append(g)
    return np.array(out_p), np.array(out_i), np.array(out_g)


def simulate_spin_system(shifts_ppm, couplings, sf_mhz, hetero=None, weights_by_shift=None, group_tol_ppm=0.012):
    """Full compound: split into components, simulate each, return sticks in ppm.
    weights_by_shift: {ppm: weight} -> the component containing a spin at that ppm
    gets that relative weight (used for glucose alpha/beta anomers)."""
    n = len(shifts_ppm)
    shz = [s * sf_mhz for s in shifts_ppm]
    comps = _components(n, couplings, shz)
    P, I, C = [], [], []
    # proton groups: spins of one component whose shifts are within group_tol_ppm
    group_of = {}
    gid = 0
    for comp in comps:
        order = sorted(comp, key=lambda g: shifts_ppm[g])
        prev = None
        for g in order:
            if prev is None or shifts_ppm[g] - shifts_ppm[prev] > group_tol_ppm:
                gid += 1
            group_of[g] = gid
            prev = g
    # merge groups from different components that sit at the same shift
    # (e.g. the three protons of a methyl, the two halves of citrate)
    centres = {}
    for g, gi in group_of.items():
        centres.setdefault(gi, []).append(shifts_ppm[g])
    centres = {gi: float(np.mean(v)) for gi, v in centres.items()}
    remap = {}
    for gi in sorted(centres, key=lambda q: centres[q]):
        for gj in list(remap.values()):
            if abs(centres[gj] - centres[gi]) < 5e-4:
                remap[gi] = gj
                break
        else:
            remap[gi] = gi
    group_of = {g: remap[gi] for g, gi in group_of.items()}
    for ci, comp in enumerate(comps):
        loc = {g: l for l, g in enumerate(comp)}
        cc = {(loc[i], loc[j]): J for (i, j), J in couplings.items() if i in loc and j in loc}
        het = {loc[i]: J for i, J in (hetero or {}).items() if i in loc} or None
        ref = np.mean([shz[g] for g in comp])
        f, it, own = simulate_component([shz[g] - ref for g in comp], cc, het)
        C.append(np.array([group_of[comp[o]] for o in own], int))
        w = 1.0
        if weights_by_shift:
            for sp, wt in weights_by_shift.items():
                if any(abs(shifts_ppm[g] - sp) < 1e-3 for g in comp):
                    w = wt * len(comps) / 1.0  # rescaled below
        P.append((f + ref) / sf_mhz)
        I.append(it * (w if weights_by_shift else 1.0))
    ppm = np.concatenate(P)
    inten = np.concatenate(I)
    groups = np.concatenate(C)
    if weights_by_shift:
        # normalise so that the molecule carries (number of spins per anomer) protons
        n_per = n / len(comps)
        inten = inten * (n_per / inten.sum())
    return ppm, inten, groups
