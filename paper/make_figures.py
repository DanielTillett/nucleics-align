#!/usr/bin/env python3
"""Generate publication-quality SVG figures for the Nucleics Align v1.3.0 manuscript
directly from archived pop-os benchmark measurements (zero external dependencies)."""
import math
import os

OUT = os.path.join(os.path.dirname(__file__), "figures")
os.makedirs(OUT, exist_ok=True)

FONT = "Inter, -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif"
C_NAVY = "#1b365d"
C_BLUE = "#2563eb"
C_TEAL = "#0d9488"
C_AMBER = "#d97706"
C_RED = "#dc2626"
C_PURPLE = "#7c3aed"
C_GRAY = "#64748b"
C_LIGHT = "#e2e8f0"
C_BG = "#ffffff"
C_TEXT = "#0f172a"

def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

def svg_header(w, h):
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">\n'
        f'  <defs>\n'
        f'    <marker id="arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">\n'
        f'      <path d="M 0 1.5 L 10 5 L 0 8.5 z" fill="{C_GRAY}"/>\n'
        f'    </marker>\n'
        f'  </defs>\n'
        f'  <rect width="{w}" height="{h}" fill="{C_BG}"/>\n'
        f'  <style>\n'
        f'    text {{ font-family: {FONT}; fill: {C_TEXT}; }}\n'
        f'    .title {{ font-size: 13px; font-weight: 700; }}\n'
        f'    .sub {{ font-size: 10.5px; fill: {C_GRAY}; }}\n'
        f'    .lbl {{ font-size: 10px; }}\n'
        f'    .lblb {{ font-size: 10px; font-weight: 600; }}\n'
        f'    .tick {{ font-size: 9.5px; fill: #334155; }}\n'
        f'    .small {{ font-size: 8.8px; fill: #475569; }}\n'
        f'  </style>\n'
    )

def make_fig1():
    """Figure 1: Nucleics Align v1.3.0 execution dispatch & exactness architecture."""
    w, h = 860, 350
    s = [svg_header(w, h)]
    s.append(f'  <text x="20" y="24" class="title">Nucleics Align v1.3.0 — Exact Execution Dispatch Architecture</text>')
    s.append(f'  <text x="20" y="40" class="sub">All default matrix paths preserve byte-identical scalar Gotoh tie-breaking; --linear, --wfa and --fastlocal guarantee identical optimal scores.</text>')

    def box(x, y, bw, bh, title, lines, fill="#f8fafc", stroke=C_NAVY, title_color=C_NAVY):
        s.append(f'  <rect x="{x}" y="{y}" width="{bw}" height="{bh}" rx="6" fill="{fill}" stroke="{stroke}" stroke-width="1.4"/>')
        s.append(f'  <text x="{x+bw/2}" y="{y+18}" text-anchor="middle" class="lblb" style="fill:{title_color}">{esc(title)}</text>')
        for idx, line in enumerate(lines):
            s.append(f'  <text x="{x+bw/2}" y="{y+33+idx*13}" text-anchor="middle" class="small">{esc(line)}</text>')

    def arrow(x1, y1, x2, y2, label=""):
        s.append(f'  <line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{C_GRAY}" stroke-width="1.3" marker-end="url(#arr)"/>')
        if label:
            mx, my = (x1 + x2) / 2, (y1 + y2) / 2 - 4
            s.append(f'  <text x="{mx}" y="{my}" text-anchor="middle" class="small">{esc(label)}</text>')

    # Entry box
    box(20, 140, 145, 72, "Pairwise Input", ["DNA or protein (A, B)", "Affine scoring + matrix", "Checked C11 / CLI API"], fill="#eff6ff", stroke=C_BLUE, title_color=C_BLUE)

    # 4 mode branches
    box(215, 56, 175, 58, "--wfa (BiWFA)", ["O(s) trimmed fronts, 32B AVX2", "--wfabudget (auto <= 2^28)"], fill="#ecfdf5", stroke=C_TEAL, title_color=C_TEAL)
    box(215, 130, 175, 58, "--linear (Hirschberg)", ["O(|B|) rows, 16-bit + reuse", "Myers split rows (0/-1/0/-1)"], fill="#f5f3ff", stroke=C_PURPLE, title_color=C_PURPLE)
    box(215, 204, 175, 58, "--scoreonly", ["Ukkonen-banded Myers (unit)", "or banded/MT striped SIMD"], fill="#fffbeb", stroke=C_AMBER, title_color=C_AMBER)
    box(215, 278, 175, 54, "smith_waterman --maxhits 1", ["O(|B|) scans + lane hints", "8/16b strips or opt-in --fastlocal"], fill="#fef2f2", stroke=C_RED, title_color=C_RED)

    # Default matrix pipeline on right
    box(445, 72, 185, 72, "1. Band Sizing & Certificate", ["Budgeted score-only WFA or 12-mer", "Strict S_out(b) < L certificate", "Excludes tied paths outside B(b)"], fill="#f8fafc", stroke=C_NAVY, title_color=C_NAVY)
    box(445, 170, 185, 72, "2. Striped Diff-Encoded Fill", ["1-byte/cell traceback (12x RAM cut)", "8/16-bit diff lanes + prefix-scan", "64B padded condvar MT pipeline"], fill="#f8fafc", stroke=C_NAVY, title_color=C_NAVY)
    box(665, 72, 175, 72, "Strip Pruning & Strip-Kept", ["Upper bound U(x,y) < L pruning", "Auto strip-kept if >25% RAM", "1 Mb 1% in 342 MiB (exact ties)"], fill="#f0f9ff", stroke=C_BLUE, title_color=C_BLUE)
    box(665, 170, 175, 72, "Exact Output & Fallbacks", ["Byte-identical Gotoh traceback", "Unsupported -> 3-matrix scalar", "--linearfallback on OOM"], fill="#f8fafc", stroke=C_NAVY, title_color=C_NAVY)

    # Arrows from input
    for ty in (85, 159, 233, 305):
        s.append(f'  <path d="M 165 176 L 190 176 L 190 {ty} L 213 {ty}" fill="none" stroke="{C_GRAY}" stroke-width="1.3" marker-end="url(#arr)"/>')
    # Fallback from --wfa to default matrix
    arrow(390, 85, 443, 102, "fallback")
    # Default path from input around to matrix
    s.append(f'  <path d="M 390 159 L 418 159 L 418 122 L 443 122" fill="none" stroke="{C_GRAY}" stroke-width="1.2" stroke-dasharray="4,3" marker-end="url(#arr)"/>')
    s.append(f'  <text x="418" y="152" text-anchor="middle" class="small">default</text>')
    arrow(537, 144, 537, 168)
    s.append(f'  <text x="544" y="159" text-anchor="start" class="small">band ok / full</text>')
    arrow(630, 108, 663, 108, "strips")
    arrow(630, 206, 663, 206)
    arrow(752, 144, 752, 168)

    s.append('</svg>\n')
    with open(os.path.join(OUT, "fig1_architecture.svg"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(s))

def make_fig2():
    """Figure 2: (A) Fresh-process 9-config ablation (20 kb 1% & 10%); (B) 50 kb Speed vs Peak RAM Pareto."""
    w, h = 880, 390
    s = [svg_header(w, h)]

    # --- Panel A: Ablation horizontal log-scale bar chart ---
    s.append(f'  <text x="16" y="22" class="title">(a) Fresh-process 9-configuration ablation (20 kb DNA, 5 seeds)</text>')
    s.append(f'  <text x="16" y="36" class="sub">Median library wall time (log scale, ms) with median peak RSS annotated</text>')

    rows = [
        ("1. Scalar 3-mat (1t)",      2728.8, 2828.5, "4.58 GiB"),
        ("2. Compact+32b SIMD",        443.4,  453.1, "386 MiB"),
        ("3. + 8-bit diff lanes",      298.8,  296.9, "385 MiB"),
        ("4. + Band (12-mer)",          22.1,  302.8, "17 / 385 MiB"),
        ("5. + Band (WFA probe)",       24.3,  303.4, "18 / 386 MiB"),
        ("6. + Multi-thread (8t)",      24.2,   78.9, "18 / 387 MiB"),
        ("7. --linear (1t)",            39.6,  160.9, "2.8 / 3.2 MiB"),
        ("8. --linear (8t)",            28.6,   65.9, "3.2 / 3.6 MiB"),
        ("9. --wfa BiWFA (1t)",          7.8,  332.6, "2.9 / 386 MiB"),
    ]
    ax0, ax1, ay0, ay1 = 148, 435, 58, 336
    lmin, lmax = math.log10(5.0), math.log10(5000.0)
    def lx(v):
        return ax0 + (math.log10(v) - lmin) / (lmax - lmin) * (ax1 - ax0)

    # Gridlines
    for tick in [10, 100, 1000]:
        tx = lx(tick)
        s.append(f'  <line x1="{tx:.1f}" y1="{ay0}" x2="{tx:.1f}" y2="{ay1}" stroke="{C_LIGHT}" stroke-width="1"/>')
        s.append(f'  <text x="{tx:.1f}" y="{ay1+14}" text-anchor="middle" class="tick">{tick}</text>')
    s.append(f'  <line x1="{ax0}" y1="{ay1}" x2="{ax1}" y2="{ay1}" stroke="#94a3b8" stroke-width="1"/>')
    s.append(f'  <text x="{(ax0+ax1)/2:.1f}" y="{ay1+30}" text-anchor="middle" class="lblb">Median wall time (ms, log scale)</text>')

    step = (ay1 - ay0) / len(rows)
    for i, (name, v1, v10, rss) in enumerate(rows):
        y = ay0 + i * step + 4
        s.append(f'  <text x="{ax0-6}" y="{y+13}" text-anchor="end" class="lbl">{esc(name)}</text>')
        # 1% bar (teal)
        w1 = max(2.0, lx(v1) - ax0)
        s.append(f'  <rect x="{ax0}" y="{y}" width="{w1:.1f}" height="9" rx="2" fill="{C_TEAL}"/>')
        # 10% bar (blue)
        w10 = max(2.0, lx(v10) - ax0)
        s.append(f'  <rect x="{ax0}" y="{y+10}" width="{w10:.1f}" height="9" rx="2" fill="{C_BLUE}"/>')
        # RSS label
        mx = max(lx(v1), lx(v10))
        if mx < ax1 - 55:
            s.append(f'  <text x="{mx+4:.1f}" y="{y+13}" class="small">{esc(rss)}</text>')
        else:
            s.append(f'  <text x="{mx-4:.1f}" y="{y+8}" text-anchor="end" class="small" style="fill:#fff">{esc(rss)}</text>')

    # Legend A
    s.append(f'  <rect x="{ax0}" y="42" width="10" height="8" rx="1.5" fill="{C_TEAL}"/>')
    s.append(f'  <text x="{ax0+14}" y="49" class="small">1% sub.</text>')
    s.append(f'  <rect x="{ax0+66}" y="42" width="10" height="8" rx="1.5" fill="{C_BLUE}"/>')
    s.append(f'  <text x="{ax0+80}" y="49" class="small">10% sub.</text>')

    # --- Panel B: 50 kb Comparator Speed vs Peak RSS Pareto ---
    bx0, bx1, by0, by1 = 520, 860, 58, 336
    s.append(f'  <text x="480" y="22" class="title">(b) 50 kb global comparators: wall time vs peak RSS</text>')
    s.append(f'  <text x="480" y="36" class="sub">Bottom-left is faster and smaller (5-seed medians, log-log)</text>')

    xmin, xmax = math.log10(10.0), math.log10(8000.0)
    ymin, ymax = math.log10(2.0), math.log10(15000.0)
    def px(v): return bx0 + (math.log10(v) - xmin) / (xmax - xmin) * (bx1 - bx0)
    def py(v): return by1 - (math.log10(v) - ymin) / (ymax - ymin) * (by1 - by0)

    for tx_v in [10, 100, 1000]:
        x = px(tx_v)
        s.append(f'  <line x1="{x:.1f}" y1="{by0}" x2="{x:.1f}" y2="{by1}" stroke="{C_LIGHT}" stroke-width="1"/>')
        s.append(f'  <text x="{x:.1f}" y="{by1+14}" text-anchor="middle" class="tick">{tx_v}</text>')
    for ty_v, lbl in [(4, "4"), (16, "16"), (64, "64"), (256, "256"), (1024, "1 GiB"), (4096, "4 GiB")]:
        y = py(ty_v)
        s.append(f'  <line x1="{bx0}" y1="{y:.1f}" x2="{bx1}" y2="{y:.1f}" stroke="{C_LIGHT}" stroke-width="1"/>')
        s.append(f'  <text x="{bx0-6}" y="{y+3:.1f}" text-anchor="end" class="tick">{lbl}</text>')
    s.append(f'  <rect x="{bx0}" y="{by0}" width="{bx1-bx0}" height="{by1-by0}" fill="none" stroke="#94a3b8" stroke-width="1"/>')
    s.append(f'  <text x="{(bx0+bx1)/2:.1f}" y="{by1+30}" text-anchor="middle" class="lblb">Median wall time (ms, log scale)</text>')
    s.append(f'  <text x="476" y="{(by0+by1)/2:.1f}" text-anchor="middle" class="lblb" transform="rotate(-90 476 {(by0+by1)/2:.1f})">Peak RSS (MiB, log)</text>')

    pts = [
        # (wall_ms, rss_mib, color, shape, label, dx, dy)
        (38.7,   3.4,  C_TEAL,  "circle", "NA --wfa (1%,1t)",   6, 12),
        (158.5,  4.1,  C_BLUE,  "circle", "NA --lin (1%,1t)",    7,  4),
        (835.4,  4.9,  C_BLUE,  "square", "NA --lin (10%,1t)",   7,  4),
        (224.0,  6.4,  C_PURPLE,"square", "NA --lin (10%,8t)",  7, -6),
        (18.6,   5.5,  C_AMBER, "tri",    "WFA2 BiWFA (1%)",    6, -5),
        (860.5,  18.6, C_AMBER, "tri",    "WFA2 BiWFA (10%)", -32, -9),
        (68.6,  121.0, C_RED,   "tri",    "WFA2 high (1%)",     6,  3),
        (5069.8,9225.4,C_RED,   "tri",    "WFA2 high (10%)",  -86, 14),
        (2357.9,2759.9,C_GRAY,  "diamond","ksw2 (1% & 10%)",  -42, 16),
        (429.8, 2393.3,C_NAVY,  "square", "NA full (10%,8t)", -84,  4),
    ]
    for vm, vr, col, shp, lbl, dx, dy in pts:
        x, y = px(vm), py(vr)
        if shp == "circle":
            s.append(f'  <circle cx="{x:.1f}" cy="{y:.1f}" r="4.5" fill="{col}" stroke="#fff" stroke-width="1"/>')
        elif shp == "square":
            s.append(f'  <rect x="{x-4:.1f}" y="{y-4:.1f}" width="8" height="8" fill="{col}" stroke="#fff" stroke-width="1"/>')
        elif shp == "tri":
            s.append(f'  <polygon points="{x:.1f},{y-5:.1f} {x-4.5:.1f},{y+4:.1f} {x+4.5:.1f},{y+4:.1f}" fill="{col}" stroke="#fff" stroke-width="0.8"/>')
        else:
            s.append(f'  <polygon points="{x:.1f},{y-5:.1f} {x+5:.1f},{y:.1f} {x:.1f},{y+5:.1f} {x-5:.1f},{y:.1f}" fill="{col}" stroke="#fff" stroke-width="0.8"/>')
        s.append(f'  <text x="{x+dx:.1f}" y="{y+dy:.1f}" class="small" style="fill:{col};font-weight:600">{esc(lbl)}</text>')

    s.append('</svg>\n')
    with open(os.path.join(OUT, "fig2_ablation_compare.svg"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(s))

def make_fig3():
    """Figure 3: (A) 5 kb - 1 Mb log-log scaling; (B) 5-protein BLOSUM62 benchmark vs seq-align & parasail."""
    w, h = 880, 380
    s = [svg_header(w, h)]

    # --- Panel A: 5 kb - 1 Mb scaling ---
    ax0, ax1, ay0, ay1 = 62, 420, 58, 326
    s.append(f'  <text x="16" y="22" class="title">(a) Scaling from 5 kb to 1 Mb (v1.3.0, Table 3)</text>')
    s.append(f'  <text x="16" y="36" class="sub">Wall time vs sequence length (log-log; 1% solid, 0.1% dashed)</text>')

    xmin, xmax = math.log10(4.0), math.log10(1300.0)
    ymin, ymax = math.log10(2.0), math.log10(15000.0)
    def px(v): return ax0 + (math.log10(v) - xmin) / (xmax - xmin) * (ax1 - ax0)
    def py(v): return ay1 - (math.log10(v) - ymin) / (ymax - ymin) * (ay1 - ay0)

    for xv, xl in [(5,"5k"),(20,"20k"),(50,"50k"),(100,"100k"),(300,"300k"),(1000,"1M")]:
        x = px(xv)
        s.append(f'  <line x1="{x:.1f}" y1="{ay0}" x2="{x:.1f}" y2="{ay1}" stroke="{C_LIGHT}" stroke-width="1"/>')
        s.append(f'  <text x="{x:.1f}" y="{ay1+14}" text-anchor="middle" class="tick">{xl}</text>')
    for yv, yl in [(10,"10 ms"),(100,"100 ms"),(1000,"1 s"),(10000,"10 s")]:
        y = py(yv)
        s.append(f'  <line x1="{ax0}" y1="{y:.1f}" x2="{ax1}" y2="{y:.1f}" stroke="{C_LIGHT}" stroke-width="1"/>')
        s.append(f'  <text x="{ax0-6}" y="{y+3:.1f}" text-anchor="end" class="tick">{yl}</text>')
    s.append(f'  <rect x="{ax0}" y="{ay0}" width="{ax1-ax0}" height="{ay1-ay0}" fill="none" stroke="#94a3b8" stroke-width="1"/>')
    s.append(f'  <text x="{(ax0+ax1)/2:.1f}" y="{ay1+30}" text-anchor="middle" class="lblb">Sequence length (bp, log scale)</text>')
    s.append(f'  <text x="16" y="{(ay0+ay1)/2:.1f}" text-anchor="middle" class="lblb" transform="rotate(-90 16 {(ay0+ay1)/2:.1f})">Wall time (log)</text>')

    series = [
        ("full 1% (8t)",        C_NAVY,  "",    [(5,6.3),(20,25.6),(50,125.8),(100,476.7),(300,969.8),(1000,9349.3)]),
        ("--linear 1% (8t)",    C_BLUE,  "",    [(5,7.1),(20,30.1),(50,105.4),(100,304.8),(300,1170.1),(1000,7293.5)]),
        ("--wfa 1% (1t)",       C_TEAL,  "",    [(5,2.8),(20,9.5),(50,39.4),(100,136.1)]),
        ("--wfa 0.1% (1t)",     C_TEAL,  "5,3", [(100,9.9),(300,24.1),(1000,160.7)]),
        ("WFA2 BiWFA 1% (1t)",  C_AMBER, "",    [(100,57.9),(300,390.8),(1000,4009.6)]),
        ("WFA2 BiWFA 0.1% (1t)",C_AMBER, "5,3", [(100,6.9),(300,12.1),(1000,76.1)]),
    ]
    for lbl, col, dash, pts in series:
        d = " ".join(f"{'M' if i==0 else 'L'} {px(x):.1f} {py(y):.1f}" for i,(x,y) in enumerate(pts))
        da = f' stroke-dasharray="{dash}"' if dash else ""
        s.append(f'  <path d="{d}" fill="none" stroke="{col}" stroke-width="2"{da}/>')
        for x, y in pts:
            s.append(f'  <circle cx="{px(x):.1f}" cy="{py(y):.1f}" r="3.2" fill="{col}" stroke="#fff" stroke-width="0.8"/>')

    # Legend A (top-left inside plot)
    lx0, ly0 = ax0 + 8, ay0 + 10
    s.append(f'  <rect x="{lx0-4}" y="{ly0-4}" width="142" height="90" rx="4" fill="#ffffff" fill-opacity="0.92" stroke="{C_LIGHT}"/>')
    for i, (lbl, col, dash, _) in enumerate(series):
        yy = ly0 + 8 + i * 14
        da = f' stroke-dasharray="{dash}"' if dash else ""
        s.append(f'  <line x1="{lx0}" y1="{yy}" x2="{lx0+18}" y2="{yy}" stroke="{col}" stroke-width="2"{da}/>')
        s.append(f'  <text x="{lx0+23}" y="{yy+3}" class="small">{esc(lbl)}</text>')

    # --- Panel B: 5-protein BLOSUM62 benchmark ---
    bx0, bx1, by0, by1 = 505, 860, 74, 326
    s.append(f'  <text x="465" y="22" class="title">(b) Real protein benchmark (BLOSUM62, Table 5)</text>')
    s.append(f'  <text x="465" y="36" class="sub">Median wall time (ms, log scale); 5 NCBI RefSeq orthologues</text>')

    prots = [
        ("TP53\n(393 aa)",     5.6,  2.0,  2.1,  2.0,  2.1),
        ("POLR2A\n(1,970 aa)",85.9,  6.2,  6.0,  5.9,  5.9),
        ("BRCA1\n(1,863 aa)", 85.8,  6.9,  6.2,  5.8,  5.9),
        ("MTOR\n(2,549 aa)", 152.1,  9.9,  9.5,  8.3,  6.5),
        ("ATM\n(3,056 aa)",  226.5, 13.1, 12.4, 10.4,  8.1),
    ]
    pmin, pmax = math.log10(1.0), math.log10(350.0)
    def bpy(v): return by1 - (math.log10(v) - pmin) / (pmax - pmin) * (by1 - by0)

    for yv in [1, 3, 10, 30, 100, 300]:
        y = bpy(yv)
        s.append(f'  <line x1="{bx0}" y1="{y:.1f}" x2="{bx1}" y2="{y:.1f}" stroke="{C_LIGHT}" stroke-width="1"/>')
        s.append(f'  <text x="{bx0-6}" y="{y+3:.1f}" text-anchor="end" class="tick">{yv}</text>')
    s.append(f'  <line x1="{bx0}" y1="{by1}" x2="{bx1}" y2="{by1}" stroke="#94a3b8" stroke-width="1"/>')
    s.append(f'  <text x="468" y="{(by0+by1)/2:.1f}" text-anchor="middle" class="lblb" transform="rotate(-90 468 {(by0+by1)/2:.1f})">Wall time (ms, log)</text>')

    gstep = (bx1 - bx0) / len(prots)
    colors_p = [C_GRAY, C_AMBER, C_NAVY, C_BLUE, C_TEAL]
    names_p = ["seq-align (1t)", "parasail scan (1t)", "NA full (1t)", "NA --lin (1t)", "NA --lin (8t)"]
    bw = 10
    for gi, (pname, *vals) in enumerate(prots):
        gx = bx0 + gi * gstep + (gstep - 5 * (bw + 1)) / 2
        for bi, val in enumerate(vals):
            y = bpy(val)
            x = gx + bi * (bw + 1)
            s.append(f'  <rect x="{x:.1f}" y="{y:.1f}" width="{bw}" height="{by1-y:.1f}" rx="1.5" fill="{colors_p[bi]}"/>')
        cx = bx0 + (gi + 0.5) * gstep
        for li, part in enumerate(pname.split("\n")):
            s.append(f'  <text x="{cx:.1f}" y="{by1+14+li*11}" text-anchor="middle" class="small">{esc(part)}</text>')

    # Legend B
    for bi, (nm, col) in enumerate(zip(names_p, colors_p)):
        lx = bx0 + (bi % 3) * 116
        ly = 52 if bi < 3 else 64
        s.append(f'  <rect x="{lx}" y="{ly-7}" width="9" height="7" rx="1.5" fill="{col}"/>')
        s.append(f'  <text x="{lx+12}" y="{ly}" class="small">{esc(nm)}</text>')

    s.append('</svg>\n')
    with open(os.path.join(OUT, "fig3_scaling_proteins.svg"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(s))

if __name__ == "__main__":
    make_fig1()
    make_fig2()
    make_fig3()
    print("Generated 3 SVG figures in", OUT)
