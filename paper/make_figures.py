#!/usr/bin/env python3
"""Generate publication-quality SVG figures for the Nucleics Align manuscript
directly from archived pop-os v1.3.0 benchmark measurements (zero external dependencies)."""
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
    """Figure 1: Nucleics Align execution dispatch & exactness architecture."""
    w, h = 860, 365
    s = [svg_header(w, h)]
    s.append(f'  <text x="20" y="24" class="title">Nucleics Align: Exact Execution Dispatch Architecture</text>')
    s.append(f'  <text x="20" y="40" class="sub">All default matrix paths preserve byte-identical scalar Gotoh tie-breaking; --linear, --wfa and --fastlocal guarantee identical optimal scores.</text>')

    def box(x, y, bw, bh, title, lines, fill="#f8fafc", stroke=C_NAVY, title_color=C_NAVY):
        s.append(f'  <rect x="{x}" y="{y}" width="{bw}" height="{bh}" rx="6" fill="{fill}" stroke="{stroke}" stroke-width="1.4"/>')
        s.append(f'  <text x="{x+bw/2}" y="{y+18}" text-anchor="middle" class="lblb" style="fill:{title_color}">{esc(title)}</text>')
        for idx, line in enumerate(lines):
            s.append(f'  <text x="{x+bw/2}" y="{y+33+idx*13}" text-anchor="middle" class="small">{esc(line)}</text>')

    def arrow(x1, y1, x2, y2, label="", ly_off=-4):
        s.append(f'  <line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{C_GRAY}" stroke-width="1.3" marker-end="url(#arr)"/>')
        if label:
            mx, my = (x1 + x2) / 2, (y1 + y2) / 2 + ly_off
            s.append(f'  <text x="{mx}" y="{my}" text-anchor="middle" class="small">{esc(label)}</text>')

    # Entry box
    box(16, 152, 142, 72, "Pairwise Input", ["DNA or protein (A, B)", "Affine scoring + matrix", "Checked C11 / CLI API"], fill="#eff6ff", stroke=C_BLUE, title_color=C_BLUE)

    # 4 explicit/opt-in mode branches
    box(212, 74, 178, 56, "--wfa (BiWFA)", ["O(s) trimmed fronts, 32B AVX2", "--wfabudget (default <= 2^28)"], fill="#ecfdf5", stroke=C_TEAL, title_color=C_TEAL)
    box(212, 144, 178, 56, "--linear (Hirschberg)", ["O(|B|) rows, 16-bit + reuse", "Myers split rows (0/-1/0/-1)"], fill="#f5f3ff", stroke=C_PURPLE, title_color=C_PURPLE)
    box(212, 214, 178, 56, "--scoreonly", ["Ukkonen-banded Myers (unit)", "or banded/MT striped SIMD"], fill="#fffbeb", stroke=C_AMBER, title_color=C_AMBER)
    box(212, 284, 178, 56, "smith_waterman --maxhits 1", ["O(|B|) scan memory + hints", "8/16b strips or --fastlocal"], fill="#fef2f2", stroke=C_RED, title_color=C_RED)

    # Default matrix pipeline on right
    box(448, 86, 185, 74, "1. Band Sizing & Certificate", ["Budgeted WFA or 12-mer probe", "Strict S_out(b) < L (widen on fail)", "Excludes tied paths outside B(b)"], fill="#f8fafc", stroke=C_NAVY, title_color=C_NAVY)
    box(448, 190, 185, 74, "2. Striped Diff-Encoded Fill", ["1-byte/cell (12x matrix RAM cut)", "8/16-bit diff lanes + prefix-scan", "64B padded condvar MT pipeline"], fill="#f8fafc", stroke=C_NAVY, title_color=C_NAVY)
    box(668, 86, 176, 74, "Strip Pruning & Strip-Kept", ["Upper bound U(x,y) < L pruning", "Auto strip-kept if >25% RAM", "1 Mb 1%: 342 MiB (exact ties)"], fill="#f0f9ff", stroke=C_BLUE, title_color=C_BLUE)
    box(668, 190, 176, 74, "Exact Output & Fallbacks", ["Byte-identical Gotoh traceback", "Unsupported -> 3-matrix scalar", "--linearfallback on OOM"], fill="#f8fafc", stroke=C_NAVY, title_color=C_NAVY)

    # Arrows from input dispatch junction (x=186) to the 4 explicit modes
    for ty in (102, 172, 242, 312):
        s.append(f'  <path d="M 158 188 L 186 188 L 186 {ty} L 210 {ty}" fill="none" stroke="{C_GRAY}" stroke-width="1.3" marker-end="url(#arr)"/>')

    # Default global path: routes directly from dispatch junction OVER the --wfa box to Matrix Step 1
    s.append(f'  <path d="M 186 102 L 186 58 L 540 58 L 540 84" fill="none" stroke="{C_NAVY}" stroke-width="1.3" stroke-dasharray="5,3" marker-end="url(#arr)"/>')
    s.append(f'  <text x="364" y="53" text-anchor="middle" class="small" style="fill:{C_NAVY};font-weight:600">default global alignment (or --linearabove fallback)</text>')

    # Fallback from --wfa to default matrix
    arrow(390, 102, 446, 114, "fallback", -6)

    # Matrix pipeline arrows
    arrow(540, 160, 540, 188)
    s.append(f'  <text x="547" y="177" text-anchor="start" class="small">band ok / full</text>')
    arrow(633, 123, 666, 123, "strips", -5)
    arrow(633, 227, 666, 227)
    arrow(756, 160, 756, 188)

    s.append('</svg>\n')
    with open(os.path.join(OUT, "fig1_architecture.svg"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(s))

def make_fig2():
    """Figure 2: (a) Fresh-process 9-config ablation dot-and-range plot; (b) 50 kb Speed vs Peak RSS."""
    w, h = 880, 395
    s = [svg_header(w, h)]

    # --- Panel A: Ablation horizontal log-scale dot-and-range plot ---
    s.append(f'  <text x="16" y="22" class="title">(a) Fresh-process 9-configuration ablation (20 kb DNA, 5 seeds)</text>')
    s.append(f'  <text x="16" y="36" class="sub">Median [min-max] library wall time (ms, log) with median peak RSS</text>')

    # (name, (med1, lo1, hi1), (med10, lo10, hi10), rss_label, fallback10)
    rows = [
        ("1. Scalar 3-mat (1t)",     (2728.8, 2719.0, 2736.5), (2828.5, 2819.4, 2854.0), "4.47 GiB",       False),
        ("2. Compact+32b SIMD",      ( 443.4,  440.2,  459.3), ( 453.1,  428.3,  460.0), "386 MiB",        False),
        ("3. + 8-bit diff lanes",    ( 298.8,  295.7,  299.2), ( 296.9,  295.2,  299.1), "385 MiB",        False),
        ("4. + Band (12-mer)",       (  22.1,   20.1,   23.1), ( 302.8,  302.4,  304.1), "17 / 385 MiB",   False),
        ("5. + Band (WFA probe)",    (  24.3,   22.4,   25.3), ( 303.4,  301.5,  305.6), "18 / 386 MiB",   False),
        ("6. + Multi-thread (8t)",   (  24.2,   22.1,   24.9), (  78.9,   76.1,   86.2), "18 / 387 MiB",   False),
        ("7. --linear (1t)",         (  39.6,   38.4,   40.5), ( 160.9,  142.3,  169.2), "2.8 / 3.2 MiB",  False),
        ("8. --linear (8t)",         (  28.6,   27.9,   29.5), (  65.9,   64.1,   67.9), "3.2 / 3.6 MiB",  False),
        ("9. --wfa BiWFA (1t)",      (   7.8,    7.1,    8.8), ( 332.6,  330.9,  338.8), "2.9 / 386 MiB†", True),
    ]
    ax0, ax1, ay0, ay1 = 148, 442, 58, 340
    lmin, lmax = math.log10(5.0), math.log10(6000.0)
    def lx(v):
        return ax0 + (math.log10(v) - lmin) / (lmax - lmin) * (ax1 - ax0)

    for tick in [10, 100, 1000]:
        tx = lx(tick)
        s.append(f'  <line x1="{tx:.1f}" y1="{ay0}" x2="{tx:.1f}" y2="{ay1}" stroke="{C_LIGHT}" stroke-width="1"/>')
        s.append(f'  <text x="{tx:.1f}" y="{ay1+14}" text-anchor="middle" class="tick">{tick}</text>')
    s.append(f'  <rect x="{ax0}" y="{ay0}" width="{ax1-ax0}" height="{ay1-ay0}" fill="none" stroke="#94a3b8" stroke-width="1"/>')
    s.append(f'  <text x="{(ax0+ax1)/2:.1f}" y="{ay1+30}" text-anchor="middle" class="lblb">Median wall time (ms, log scale)</text>')

    step = (ay1 - ay0) / len(rows)
    for i, (name, (m1, lo1, hi1), (m10, lo10, hi10), rss, fb10) in enumerate(rows):
        yc = ay0 + (i + 0.5) * step
        if i % 2 == 1:
            s.append(f'  <rect x="{ax0+0.5}" y="{ay0+i*step:.1f}" width="{ax1-ax0-1}" height="{step:.1f}" fill="#f8fafc"/>')
        s.append(f'  <line x1="{ax0}" y1="{yc:.1f}" x2="{ax1}" y2="{yc:.1f}" stroke="#f1f5f9" stroke-width="1"/>')
        s.append(f'  <text x="{ax0-6}" y="{yc+3.5:.1f}" text-anchor="end" class="lbl">{esc(name)}</text>')
        # 1% whisker + dot (teal, slightly above center)
        y1 = yc - 4.5
        s.append(f'  <line x1="{lx(lo1):.1f}" y1="{y1:.1f}" x2="{lx(hi1):.1f}" y2="{y1:.1f}" stroke="{C_TEAL}" stroke-width="2"/>')
        s.append(f'  <circle cx="{lx(m1):.1f}" cy="{y1:.1f}" r="4.0" fill="{C_TEAL}" stroke="#fff" stroke-width="0.9"/>')
        # 10% whisker + square (blue, slightly below center; hollow if fallback)
        y10 = yc + 4.5
        s.append(f'  <line x1="{lx(lo10):.1f}" y1="{y10:.1f}" x2="{lx(hi10):.1f}" y2="{y10:.1f}" stroke="{C_BLUE}" stroke-width="2"/>')
        if fb10:
            s.append(f'  <rect x="{lx(m10)-3.5:.1f}" y="{y10-3.5:.1f}" width="7" height="7" fill="#fff" stroke="{C_BLUE}" stroke-width="1.6"/>')
        else:
            s.append(f'  <rect x="{lx(m10)-3.5:.1f}" y="{y10-3.5:.1f}" width="7" height="7" fill="{C_BLUE}" stroke="#fff" stroke-width="0.9"/>')
        # RSS annotation
        mx = max(lx(hi1), lx(hi10))
        if mx < ax1 - 68:
            s.append(f'  <text x="{mx+6:.1f}" y="{yc+3:.1f}" class="small">{esc(rss)}</text>')
        else:
            s.append(f'  <text x="{min(lx(lo1), lx(lo10))-6:.1f}" y="{yc+3:.1f}" text-anchor="end" class="small">{esc(rss)}</text>')

    # Legend A
    s.append(f'  <circle cx="{ax0+5}" cy="46" r="3.8" fill="{C_TEAL}"/>')
    s.append(f'  <text x="{ax0+13}" y="49" class="small">1% sub.</text>')
    s.append(f'  <rect x="{ax0+60}" y="42.5" width="7" height="7" fill="{C_BLUE}"/>')
    s.append(f'  <text x="{ax0+71}" y="49" class="small">10% sub.</text>')
    s.append(f'  <rect x="{ax0+124}" y="42.5" width="7" height="7" fill="#fff" stroke="{C_BLUE}" stroke-width="1.5"/>')
    s.append(f'  <text x="{ax0+135}" y="49" class="small">† budget fallback</text>')

    # --- Panel B: 50 kb Comparator Speed vs Peak RSS ---
    bx0, bx1, by0, by1 = 524, 864, 58, 340
    s.append(f'  <text x="480" y="22" class="title">(b) 50 kb global comparators: wall time vs peak RSS</text>')
    s.append(f'  <text x="480" y="36" class="sub">Bottom-left is faster and smaller (5-seed medians; 1% circles, 10% squares)</text>')

    xmin, xmax = math.log10(10.0), math.log10(8000.0)
    ymin, ymax = math.log10(2.0), math.log10(16000.0)
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
    s.append(f'  <text x="478" y="{(by0+by1)/2:.1f}" text-anchor="middle" class="lblb" transform="rotate(-90 478 {(by0+by1)/2:.1f})">Peak RSS (MiB, log)</text>')

    # Subtle 1% (1t) Pareto line: WFA2 BiWFA (18.6, 5.5) -> Nucleics --wfa (38.7, 3.4)
    s.append(f'  <line x1="{px(18.6):.1f}" y1="{py(5.5):.1f}" x2="{px(38.7):.1f}" y2="{py(3.4):.1f}" stroke="{C_TEAL}" stroke-width="1.2" stroke-dasharray="3,2"/>')
    # Subtle 10% (all) Pareto line: Nucleics --lin 8t (224.0, 6.4) -> Nucleics --lin 1t (835.4, 4.9)
    s.append(f'  <line x1="{px(224.0):.1f}" y1="{py(6.4):.1f}" x2="{px(835.4):.1f}" y2="{py(4.9):.1f}" stroke="{C_PURPLE}" stroke-width="1.2" stroke-dasharray="3,2"/>')

    pts = [
        # (wall_ms, rss_mib, color, shape, label, dx, dy)
        (38.7,   3.4,  C_TEAL,  "circle", "Nucleics --wfa (1%,1t)",   6, 12),
        (158.5,  4.1,  C_BLUE,  "circle", "Nucleics --lin (1%,1t)",   7, 12),
        (835.4,  4.9,  C_BLUE,  "square", "Nucleics --lin (10%,1t)",  7,  4),
        (224.0,  6.4,  C_PURPLE,"square", "Nucleics --lin (10%,8t)",  7, -6),
        (18.6,   5.5,  C_AMBER, "circle", "WFA2 BiWFA (1%,1t)",       6, -6),
        (860.5,  18.6, C_AMBER, "square", "WFA2 BiWFA (10%,1t)",    -38, -9),
        (68.6,  121.0, C_RED,   "circle", "WFA2 high (1%,1t)",        6,  3),
        (5069.8,9225.4,C_RED,   "square", "WFA2 high (10%,1t)",     -92, 14),
        (2370.2,2762.0,C_GRAY,  "circle", "",                         0,  0),
        (2357.9,2759.9,C_GRAY,  "square", "ksw2 (1% & 10%,1t)",     -52, 16),
        (429.8, 2393.3,C_NAVY,  "square", "Nucleics full (10%,8t)", -100,  4),
    ]
    for vm, vr, col, shp, lbl, dx, dy in pts:
        x, y = px(vm), py(vr)
        if shp == "circle":
            s.append(f'  <circle cx="{x:.1f}" cy="{y:.1f}" r="4.3" fill="{col}" stroke="#fff" stroke-width="1"/>')
        else:
            s.append(f'  <rect x="{x-3.8:.1f}" y="{y-3.8:.1f}" width="7.6" height="7.6" fill="{col}" stroke="#fff" stroke-width="1"/>')
        if lbl:
            s.append(f'  <text x="{x+dx:.1f}" y="{y+dy:.1f}" class="small" style="fill:{col};font-weight:600">{esc(lbl)}</text>')

    s.append('</svg>\n')
    with open(os.path.join(OUT, "fig2_ablation_compare.svg"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(s))

def make_fig3():
    """Figure 3: (a) 5 kb - 1 Mb log-log scaling; (b) 5-protein BLOSUM62 dot plot."""
    w, h = 880, 385
    s = [svg_header(w, h)]

    # --- Panel A: 5 kb - 1 Mb scaling ---
    ax0, ax1, ay0, ay1 = 62, 424, 58, 326
    s.append(f'  <text x="16" y="22" class="title">(a) Scaling from 5 kb to 1 Mb (Table 3)</text>')
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
    # Divider between 3-run min (<=50k) and 1-run (>=100k)
    xdiv = px(70.7)
    s.append(f'  <line x1="{xdiv:.1f}" y1="{ay0}" x2="{xdiv:.1f}" y2="{ay1}" stroke="#94a3b8" stroke-width="1" stroke-dasharray="3,3"/>')
    s.append(f'  <text x="{xdiv-4:.1f}" y="{ay1-6}" text-anchor="end" class="small">min of 3</text>')
    s.append(f'  <text x="{xdiv+4:.1f}" y="{ay1-6}" text-anchor="start" class="small">1 run</text>')

    s.append(f'  <rect x="{ax0}" y="{ay0}" width="{ax1-ax0}" height="{ay1-ay0}" fill="none" stroke="#94a3b8" stroke-width="1"/>')
    s.append(f'  <text x="{(ax0+ax1)/2:.1f}" y="{ay1+30}" text-anchor="middle" class="lblb">Sequence length (bp, log scale)</text>')
    s.append(f'  <text x="16" y="{(ay0+ay1)/2:.1f}" text-anchor="middle" class="lblb" transform="rotate(-90 16 {(ay0+ay1)/2:.1f})">Wall time (log)</text>')

    series = [
        ("Nucleics full 1% (8t)",     C_NAVY,  "",    [(5,6.3),(20,25.6),(50,125.8),(100,476.7),(300,969.8),(1000,9349.3)]),
        ("Nucleics --linear 1% (8t)", C_BLUE,  "",    [(5,7.1),(20,30.1),(50,105.4),(100,304.8),(300,1170.1),(1000,7293.5)]),
        ("Nucleics --wfa 1% (1t)",    C_TEAL,  "",    [(5,2.8),(20,9.5),(50,39.4),(100,136.1)]),
        ("Nucleics --wfa 0.1% (1t)",  C_TEAL,  "5,3", [(100,9.9),(300,24.1),(1000,160.7)]),
        ("WFA2 BiWFA 1% (1t)",        C_AMBER, "",    [(100,57.9),(300,390.8),(1000,4009.6)]),
        ("WFA2 BiWFA 0.1% (1t)",      C_AMBER, "5,3", [(100,6.9),(300,12.1),(1000,76.1)]),
    ]
    for lbl, col, dash, pts in series:
        d = " ".join(f"{'M' if i==0 else 'L'} {px(x):.1f} {py(y):.1f}" for i,(x,y) in enumerate(pts))
        da = f' stroke-dasharray="{dash}"' if dash else ""
        s.append(f'  <path d="{d}" fill="none" stroke="{col}" stroke-width="2"{da}/>')
        for x, y in pts:
            s.append(f'  <circle cx="{px(x):.1f}" cy="{py(y):.1f}" r="3.2" fill="{col}" stroke="#fff" stroke-width="0.8"/>')

    # Annotate --wfa 1% fallback past 100k
    s.append(f'  <text x="{px(100)+6:.1f}" y="{py(136.1)-3:.1f}" text-anchor="start" class="small" style="fill:{C_TEAL};font-weight:600">† fallback (&gt;=300k)</text>')

    # Legend A (top-left inside plot)
    lx0, ly0 = ax0 + 8, ay0 + 10
    s.append(f'  <rect x="{lx0-4}" y="{ly0-4}" width="160" height="90" rx="4" fill="#ffffff" fill-opacity="0.93" stroke="{C_LIGHT}"/>')
    for i, (lbl, col, dash, _) in enumerate(series):
        yy = ly0 + 8 + i * 14
        da = f' stroke-dasharray="{dash}"' if dash else ""
        s.append(f'  <line x1="{lx0}" y1="{yy}" x2="{lx0+18}" y2="{yy}" stroke="{col}" stroke-width="2"{da}/>')
        s.append(f'  <text x="{lx0+23}" y="{yy+3}" class="small">{esc(lbl)}</text>')

    # --- Panel B: 5-protein BLOSUM62 dot plot ---
    bx0, bx1, by0, by1 = 508, 864, 74, 326
    s.append(f'  <text x="468" y="22" class="title">(b) Real protein benchmark (BLOSUM62, Table 5)</text>')
    s.append(f'  <text x="468" y="36" class="sub">Median wall time (ms, log scale); 5 NCBI RefSeq orthologues</text>')

    prots = [
        ("TP53\n393/390 aa",     5.6,  2.0,  2.1,  2.0,  2.1),
        ("POLR2A\n1,970/1,733", 85.9,  6.2,  6.0,  5.9,  5.9),
        ("BRCA1\n1,863/1,812",  85.8,  6.9,  6.2,  5.8,  5.9),
        ("MTOR\n2,549/2,549",  152.1,  9.9,  9.5,  8.3,  6.5),
        ("ATM\n3,056/3,066",   226.5, 13.1, 12.4, 10.4,  8.1),
    ]
    pmin, pmax = math.log10(1.2), math.log10(350.0)
    def bpy(v): return by1 - (math.log10(v) - pmin) / (pmax - pmin) * (by1 - by0)

    for yv in [2, 5, 10, 30, 100, 300]:
        y = bpy(yv)
        s.append(f'  <line x1="{bx0}" y1="{y:.1f}" x2="{bx1}" y2="{y:.1f}" stroke="{C_LIGHT}" stroke-width="1"/>')
        s.append(f'  <text x="{bx0-6}" y="{y+3:.1f}" text-anchor="end" class="tick">{yv}</text>')
    s.append(f'  <rect x="{bx0}" y="{by0}" width="{bx1-bx0}" height="{by1-by0}" fill="none" stroke="#94a3b8" stroke-width="1"/>')
    s.append(f'  <text x="470" y="{(by0+by1)/2:.1f}" text-anchor="middle" class="lblb" transform="rotate(-90 470 {(by0+by1)/2:.1f})">Median wall time (ms, log)</text>')

    gstep = (bx1 - bx0) / len(prots)
    colors_p = [C_GRAY, C_AMBER, C_NAVY, C_BLUE, C_TEAL]
    names_p = ["seq-align (1t)", "parasail scan (1t)", "Nucleics full (1t)", "Nucleics --lin (1t)", "Nucleics --lin (8t)"]
    offsets = [-14, -7, 0, 7, 14]
    for gi, (pname, *vals) in enumerate(prots):
        cx = bx0 + (gi + 0.5) * gstep
        s.append(f'  <line x1="{cx:.1f}" y1="{by0}" x2="{cx:.1f}" y2="{by1}" stroke="#f1f5f9" stroke-width="1"/>')
        # Subtle vertical span connecting min to max
        s.append(f'  <line x1="{cx:.1f}" y1="{bpy(max(vals)):.1f}" x2="{cx:.1f}" y2="{bpy(min(vals)):.1f}" stroke="#cbd5e1" stroke-width="1.2"/>')
        for bi, val in enumerate(vals):
            x = cx + offsets[bi]
            y = bpy(val)
            if bi == 0:
                s.append(f'  <polygon points="{x:.1f},{y-4.2:.1f} {x+4.2:.1f},{y:.1f} {x:.1f},{y+4.2:.1f} {x-4.2:.1f},{y:.1f}" fill="{colors_p[bi]}" stroke="#fff" stroke-width="0.8"/>')
            elif bi == 1:
                s.append(f'  <polygon points="{x:.1f},{y-4.5:.1f} {x-4.2:.1f},{y+3.8:.1f} {x+4.2:.1f},{y+3.8:.1f}" fill="{colors_p[bi]}" stroke="#fff" stroke-width="0.8"/>')
            else:
                s.append(f'  <circle cx="{x:.1f}" cy="{y:.1f}" r="4.0" fill="{colors_p[bi]}" stroke="#fff" stroke-width="0.8"/>')
        for li, part in enumerate(pname.split("\n")):
            s.append(f'  <text x="{cx:.1f}" y="{by1+14+li*11}" text-anchor="middle" class="small">{esc(part)}</text>')

    # Legend B
    for bi, (nm, col) in enumerate(zip(names_p, colors_p)):
        lx = bx0 + (bi % 3) * 120
        ly = 52 if bi < 3 else 64
        if bi == 0:
            s.append(f'  <polygon points="{lx+4},{ly-7} {lx+8},{ly-3} {lx+4},{ly+1} {lx},{ly-3}" fill="{col}"/>')
        elif bi == 1:
            s.append(f'  <polygon points="{lx+4},{ly-7} {lx},{ly+1} {lx+8},{ly+1}" fill="{col}"/>')
        else:
            s.append(f'  <circle cx="{lx+4}" cy="{ly-3}" r="3.6" fill="{col}"/>')
        s.append(f'  <text x="{lx+12}" y="{ly}" class="small">{esc(nm)}</text>')

    s.append('</svg>\n')
    with open(os.path.join(OUT, "fig3_scaling_proteins.svg"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(s))

if __name__ == "__main__":
    make_fig1()
    make_fig2()
    make_fig3()
    print("Generated 3 SVG figures in", OUT)
