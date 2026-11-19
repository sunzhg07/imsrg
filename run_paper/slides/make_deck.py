#!/usr/bin/env python3
"""MR-EOM rho3 / charge-exchange fixes: old vs new results deck (python-pptx)."""
from pptx import Presentation
from pptx.chart.data import CategoryChartData
from pptx.dml.color import RGBColor
from pptx.enum.chart import XL_CHART_TYPE, XL_LEGEND_POSITION, XL_LABEL_POSITION
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.util import Inches, Pt

INK = RGBColor(0x14, 0x21, 0x3D)      # dark navy
INK2 = RGBColor(0x3A, 0x4A, 0x66)     # muted text
AMBER = RGBColor(0xE0, 0x9F, 0x1F)    # accent
NEW = RGBColor(0x1F, 0x8A, 0x70)      # new / fixed (teal-green)
OLD = RGBColor(0xC2, 0x4E, 0x3A)      # old / broken (brick)
REF = RGBColor(0x6B, 0x76, 0x88)      # KSHELL reference (slate)
RHO2 = RGBColor(0xD9, 0xA4, 0x41)     # rho2-only
TINT = RGBColor(0xEE, 0xF2, 0xF7)     # card background
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
HEAD, BODY = "Cambria", "Calibri"

prs = Presentation()
prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
BLANK = prs.slide_layouts[6]
W = 13.333


def bg(slide, color):
    f = slide.background.fill
    f.solid()
    f.fore_color.rgb = color


def text(slide, x, y, w, h, paras, size=16, color=INK, font=BODY, bold=False,
         align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP, space=6, name=None):
    """paras: list of str | list of (str, dict) runs | dict(bullet=..., runs=...)"""
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    if name:
        tb.name = name
    tf = tb.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    tf.margin_left = tf.margin_right = Inches(0.05)
    tf.margin_top = tf.margin_bottom = Inches(0.03)
    if isinstance(paras, str):
        paras = [paras]
    for i, p in enumerate(paras):
        para = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        para.alignment = align
        para.space_after = Pt(space)
        bullet = False
        if isinstance(p, dict):
            bullet, p = p.get("bullet", False), p["runs"]
        if isinstance(p, str):
            p = [(p, {})]
        if bullet:
            pPr = para._p.get_or_add_pPr()
            pPr.set("marL", str(int(Inches(0.28))))
            pPr.set("indent", str(int(-Inches(0.22))))
            bu = pPr.makeelement("{http://schemas.openxmlformats.org/drawingml/2006/main}buChar", {"char": "•"})
            pPr.append(bu)
        for s, o in p:
            r = para.add_run()
            r.text = s
            f = r.font
            f.name = o.get("font", font)
            f.size = Pt(o.get("size", size))
            f.bold = o.get("bold", bold)
            f.italic = o.get("italic", False)
            f.color.rgb = o.get("color", color)
    return tb


def title(slide, s, sub=None, dark=False):
    text(slide, 0.6, 0.35, W - 1.2, 0.8, [[(s, {"bold": True})]], size=32,
         color=WHITE if dark else INK, font=HEAD, name="Title")
    if sub:
        text(slide, 0.6, 1.1, W - 1.2, 0.5, sub, size=16,
             color=RGBColor(0xC8, 0xD2, 0xE0) if dark else INK2, name="Subtitle")


def card(slide, x, y, w, h, fill=TINT, name=None):
    s = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(h))
    s.adjustments[0] = 0.06
    s.fill.solid()
    s.fill.fore_color.rgb = fill
    s.line.fill.background()
    s.shadow.inherit = False
    if name:
        s.name = name
    return s


def badge(slide, x, y, label, fill=INK, d=0.55, size=16):
    s = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(x), Inches(y), Inches(d), Inches(d))
    s.fill.solid()
    s.fill.fore_color.rgb = fill
    s.line.fill.background()
    s.shadow.inherit = False
    tf = s.text_frame
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.CENTER
    r = p.add_run()
    r.text = label
    r.font.size, r.font.bold, r.font.name = Pt(size), True, BODY
    r.font.color.rgb = WHITE
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    return s


def stat(slide, x, y, w, big, label, color=NEW, big_size=44):
    card(slide, x, y, w, 2.0)
    text(slide, x + 0.25, y + 0.2, w - 0.5, 0.95, [[(big, {"bold": True})]], size=big_size,
         color=color, font=HEAD)
    text(slide, x + 0.25, y + 1.15, w - 0.5, 0.8, label, size=14, color=INK2)


def table(slide, x, y, w, rows, colw, size=13, head_fill=INK, row_h=0.42, colors=None):
    nr, nc = len(rows), len(rows[0])
    gt = slide.shapes.add_table(nr, nc, Inches(x), Inches(y), Inches(w), Inches(row_h * nr))
    t = gt.table
    for j, cw in enumerate(colw):
        t.columns[j].width = Inches(cw)
    for i, row in enumerate(rows):
        t.rows[i].height = Inches(row_h)
        for j, val in enumerate(row):
            c = t.cell(i, j)
            c.margin_left = c.margin_right = Inches(0.08)
            c.margin_top = c.margin_bottom = Inches(0.03)
            c.vertical_anchor = MSO_ANCHOR.MIDDLE
            c.fill.solid()
            c.fill.fore_color.rgb = head_fill if i == 0 else (WHITE if i % 2 else TINT)
            tf = c.text_frame
            tf.word_wrap = True
            p = tf.paragraphs[0]
            p.alignment = PP_ALIGN.LEFT if j == 0 else PP_ALIGN.CENTER
            r = p.add_run()
            r.text = str(val)
            r.font.size, r.font.name = Pt(size), BODY
            r.font.bold = i == 0
            col = WHITE if i == 0 else INK
            if colors and i > 0 and (i, j) in colors:
                col = colors[(i, j)]
                r.font.bold = True
            r.font.color.rgb = col
    return gt


def notes(slide, s):
    slide.notes_slide.notes_text_frame.text = s


def footer(slide, s):
    text(slide, 0.6, 6.95, W - 1.2, 0.35, s, size=11, color=INK2)


# ---------------------------------------------------------------- 1 title
s = prs.slides.add_slide(BLANK)
bg(s, INK)
text(s, 0.8, 1.9, 11.5, 1.6, [[("MR-EOM: ρ₃ norm terms and charge-exchange fixes", {"bold": True})]],
     size=40, color=WHITE, font=HEAD)
text(s, 0.8, 3.55, 11.5, 0.6, "What changed since origin/devel 4d501bfc, why, and old vs new results",
     size=20, color=RGBColor(0xC8, 0xD2, 0xE0))
text(s, 0.8, 5.6, 11.5, 0.9, ["He8 / Li8 at emax 3  ·  Ne20 at emax 4  ·  EM1.8/2.0, ħω = 16 MeV",
                              "All changes are uncommitted in the working tree (src/EOM.cc, EOM.hh, pyIMSRG.cc)"],
     size=16, color=AMBER)

# ---------------------------------------------------------------- 2 bottom line
s = prs.slides.add_slide(BLANK)
bg(s, WHITE)
title(s, "Bottom line", "The metric is now exact and positive; charge exchange is no longer overbound")
stat(s, 0.6, 1.8, 3.85, "≤ 1.4e-11", "Kernel vs commutator, full values, ρ₃ on/off, all He8 blocks incl. itz = 1")
stat(s, 4.74, 1.8, 3.85, "5 → 0", "Negative metric modes with ρ₃ in He8 → Li8 2⁺ (old run vs new run)")
stat(s, 8.88, 1.8, 3.85, "−34 → −4.9", "Lowest Li8 2⁺ with ρ₃, MeV rel. to He8 gs (KSHELL −11.3)", color=OLD,
     big_size=40)
text(s, 0.6, 4.2, 12.1, 2.6, [
    {"bullet": True, "runs": [("ρ₃ terms added to the scalar and tensor norm kernels; ", {"bold": True}),
                              ("verified against the commutator route.", {})]},
    {"bullet": True, "runs": [("Tensor ladder was dropping half the 2b configs ", {"bold": True}),
                              ("(creator pair in the stored ket). Affects itz = 0 tensor runs too.", {})]},
    {"bullet": True, "runs": [("Charge exchange now has a signed ΔTz; ", {"bold": True}),
                              ("mixing n→p and p→n caused the negative modes and the −34 MeV state.", {})]},
    {"bullet": True, "runs": [("Still open: ", {"bold": True, "color": OLD}),
                              ("Li8 now underbound by 4–6 MeV (missing vv←vv H-side ρ₃/ρ₄); "
                               "most Arnoldi solves not converged.", {})]},
], size=16, space=8)
notes(s, "Old = He8 run 214864 (ρ3 kernel added, before the B fixes). New = run 215057 (all fixes).")

# ---------------------------------------------------------------- 3 change map
s = prs.slides.add_slide(BLANK)
bg(s, WHITE)
title(s, "Eleven changes in two groups", "A: ρ₃ in the norm kernel   ·   B: charge-changing and tensor representation")
A = [("A1", "Scalar pphv vv←vc ρ₃ term", "ConstructNormMatrix"),
     ("A2", "Tensor ρ₃ terms (3 blocks, AMC)", "ConstructNormMatrix_tensor"),
     ("A3", "ρ₃ piece in the commutator gold", "NormMultiref_tensor"),
     ("A4", "J=0, itz≠0 dispatched to tts", "NormMultiref"),
     ("A5", "Pauli-forbidden triples → 0", "RdmThreeBody_J")]
B = [("B1", "One orientation rule for stored 2b", "StoredExcitationOrientation"),
     ("B2", "Ladder keeps swapped excitations", "GetVSEOM_ladder_multiref"),
     ("B3", "Signed ΔTz for itz ≠ 0", "SetDeltaTz, ConstructConfigs_tensor"),
     ("B4", "vv←vv transfer off when invalid", "ConstructConfigs_tensor"),
     ("B5", "Kernel built on a canonical view", "ConstructNormMatrix_tensor"),
     ("B6", "Projector keeps λ > 0 only", "ConstructProjectMatrix_blocks, SqrtMat")]
for col, (items, color, x) in enumerate(((A, INK, 0.6), (B, NEW, 6.85))):
    for k, (tag, what, where) in enumerate(items):
        y = 1.8 + k * 0.82
        card(s, x, y, 5.9, 0.7)
        badge(s, x + 0.12, y + 0.1, tag, fill=color, d=0.5, size=13)
        text(s, x + 0.75, y + 0.04, 5.05, 0.36, [[(what, {"bold": True})]], size=15)
        text(s, x + 0.75, y + 0.36, 5.05, 0.3, [[(where, {"font": "Courier New"})]], size=11, color=INK2)
footer(s, "Plus: bindings SetDeltaTz / StoredExcitationOrientation; tests and pipelines under run_paper/; "
          "AMC inputs under learn/amc_tts/eom_norm/")

# ---------------------------------------------------------------- 4 A: rho3
s = prs.slides.add_slide(BLANK)
bg(s, WHITE)
title(s, "A · ρ₃ terms in the norm kernel",
      "With 2b excitations the metric needs ρ₃; only pphv vv←vc and ppvv qv←vv reach it (plus vv←vv for ΔTz ≠ 0)")
rows = [["Change", "Where", "Why", "Check"],
        ["A1 scalar pphv vv←vc", "EOM.cc ~1026", "Missing diagram; mirrors ppvv", "Ne20 e4: ≤ 7e-11"],
        ["A2 tensor pphv vv←cv", "EOM.cc ~1810", "Tensor kernel had no ρ₃", "λ = 0–4: match"],
        ["A2 tensor ppvv qv←vv", "EOM.cc ~1876", "same", "λ = 0–4: match"],
        ["A2 tensor ppvv vv←vv", "EOM.cc ~1938", "ΔTz≠0 pair transfer, p + h line", "He8 itz 1: ≤ 1.4e-11"],
        ["A3 comm223tts · ρ₃", "EOM.cc ~3927", "Gold for the tensor ρ₃ kernel", "—"],
        ["A4 dispatch to tts", "EOM.cc ~3896", "Scalar Commutator aborts for itz ≠ 0", "0⁺ itz 1 runs"],
        ["A5 try/catch → 0", "EOM.cc ~144", "GetME_pn threw (He6 1⁻ crash)", "He6 1⁻ runs"]]
table(s, 0.6, 1.8, 12.1, rows, [3.0, 2.0, 4.6, 2.5], size=14, row_h=0.47)
card(s, 0.6, 5.75, 12.1, 1.05)
text(s, 0.85, 5.82, 11.6, 0.95, [
    [("Reduced vs unreduced bookkeeping: ", {"bold": True}),
     ("γ = ρ̄₃ / ĵ₀; λ enters only through the 9j and (−1)^λ λ̂; AMC's leading −1 dropped (λ×λ→0 sign). "
      "TDTD for a J ≠ 0 state is divided by √(2J+1); every reference passes the trace check = C(N,3).", {})]],
     size=14)

# ---------------------------------------------------------------- 5 B1/B2
s = prs.slides.add_slide(BLANK)
bg(s, WHITE)
title(s, "B1–B2 · The ladder dropped swapped excitations",
      "Storage keeps only ch_bra ≤ ch_ket, channel = 6J + 2(Tz+1) + π — orientation depends on channel order")
# diagram
card(s, 0.6, 1.85, 6.3, 3.3)
text(s, 0.85, 1.95, 5.9, 0.4, [[("Stored element ⟨A | X | B⟩, ch(A) ≤ ch(B)", {"bold": True})]], size=15)
for i, (lab, sub, col, y) in enumerate((("main  (+1)", "creator pair = stored bra  A←B", NEW, 2.55),
                                         ("swapped  (−1)", "creator pair = stored ket  B←A", OLD, 3.55))):
    badge(s, 0.9, y + 0.05, "+" if i == 0 else "−", fill=col, d=0.55, size=20)
    text(s, 1.6, y, 5.1, 0.4, [[(lab, {"bold": True, "color": col})]], size=16)
    text(s, 1.6, y + 0.38, 5.1, 0.4, sub, size=14, color=INK2)
text(s, 0.85, 4.5, 5.9, 0.6, [[("Canonical amplitude of a swapped element = (−1)^(J_ket − J_bra) · stored",
                                  {"italic": True})]], size=13, color=INK2)
# right: the bug
text(s, 7.3, 1.85, 5.4, 0.45, [[("Old ladder: main only", {"bold": True, "color": OLD})]], size=18)
text(s, 7.3, 2.35, 5.4, 1.4, ["pp(J=2)←vv(J=0) was lost while pp(J=0)←vv(J=2) was kept. "
                              "The lost configs stayed in eom_confs as dead dimensions."], size=15)
stat(s, 7.3, 3.55, 2.6, "150/282", "pphv configs dead (He8 2⁺, emax 2)", color=OLD, big_size=30)
stat(s, 10.1, 3.55, 2.6, "153/300", "ppvv configs dead (same case)", color=OLD, big_size=30)
text(s, 0.6, 5.5, 12.1, 1.3, [
    {"bullet": True, "runs": [("Fix: ", {"bold": True, "color": NEW}),
                              ("swapped excitations written as H2·hZ·herm_phase (the rule pphh already used), "
                               "for pphv, ppvv and vv×vv; configs, ladder and kernel all use StoredExcitationOrientation.", {})]},
    {"bullet": True, "runs": [("Affects itz = 0 tensor runs too: ", {"bold": True}),
                              ("He8 t02/t20, Mg24 E1.", {})]}], size=15)

# ---------------------------------------------------------------- 6 B3/B4
s = prs.slides.add_slide(BLANK)
bg(s, WHITE)
title(s, "B3–B4 · Charge exchange needs a direction",
      "itz was |ΔTz|, so one run mixed n→p and p→n")
card(s, 0.6, 1.85, 5.9, 4.0)
text(s, 0.85, 1.95, 5.4, 3.8, [
    [("Why it broke the metric", {"bold": True, "size": 18})],
    [("N = ½⟨[S⁺, S⁻]⟩ equals ⟨χ†χ⟩ only if χ†|Ψ⟩ = 0.", {"italic": True})],
    "For vv×vv elements stored backwards (nn←pn on He8), χ|He8⟩ = 0 but χ†|He8⟩ ≠ 0, so the commutator "
    "metric gave −‖χ†Ψ‖².",
    [("Result: ", {"bold": True}), ("4–5 negative modes, the −34 MeV Li8 state, and duplicate vv←vv "
                                    "configs (which hid a factor 4 in the ρ₃ term).", {})]],
     size=15, space=10)
for k, (tag, head, body, col) in enumerate((
        ("B3", "Signed ΔTz", "SetDeltaTz(−1) = n→p (your sr_dtz sign). Configs keep one direction, so every "
                             "element has one orientation. Replaces my earlier de-duplication and the ¼ factor.", NEW),
        ("B4", "vv←vv off when invalid", "If the reference has valence particles of the created species "
                                         "(Ne20, N = Z) the block is indefinite: 103 negative modes for Ne20 1⁺. "
                                         "Excluded with a warning; He8 → Li8 keeps it.", NEW))):
    y = 1.85 + k * 2.05
    badge(s, 6.9, y + 0.05, tag, fill=col, d=0.6, size=14)
    text(s, 7.7, y, 5.0, 0.45, [[(head, {"bold": True})]], size=18)
    text(s, 7.7, y + 0.45, 5.0, 1.5, body, size=15)
text(s, 0.6, 6.1, 12.1, 0.7, [[("Check: ", {"bold": True}),
                                ("[Tz, χ] = ΔTz·χ for every non-null config — 0 failures in 1274 (1⁺) and 1598 (2⁺).", {})]],
     size=15)

# ---------------------------------------------------------------- 7 B5/B6
s = prs.slides.add_slide(BLANK)
bg(s, WHITE)
title(s, "B5–B6 · Canonical kernel and a positive projector")
for k, (tag, head, why, what, x) in enumerate((
        ("B5", "Kernel on a canonical view",
         "All Wick strings assume the creator pair is the bra; swapped configs (now present after B2) "
         "were evaluated wrongly. vv×vv ρ₂ block: kernel +5.60 vs commutator −15.95 (He8 1⁺ itz 1).",
         "Relabel swapped configs creator←annihilator, build, map back N = D·N_canon·D, "
         "D = (−1)^(J_ket − J_bra). eom_confs restored by RAII. Formulas unchanged.", 0.6),
        ("B6", "Projector keeps λ > 0 only",
         "Both projectors kept |λ|, so negative-norm modes entered the Krylov space and ArnoldiSolveH2 "
         "(positive metric assumed) turned them into spurious low or complex states.",
         "Symmetric blocks keep λ ≥ 1e-6·max|λ|; number of dropped negative modes printed. "
         "Only change to scalar 0⁺ runs.", 6.85))):
    card(s, x, 1.55, 5.9, 5.2)
    badge(s, x + 0.25, 1.75, tag, fill=NEW, d=0.65, size=15)
    text(s, x + 1.05, 1.8, 4.7, 0.6, [[(head, {"bold": True})]], size=20)
    text(s, x + 0.3, 2.65, 5.3, 4.0, [
        [("Why", {"bold": True, "color": OLD})], why,
        [("What", {"bold": True, "color": NEW})], what], size=15, space=8)

# ---------------------------------------------------------------- 8 validation
s = prs.slides.add_slide(BLANK)
bg(s, WHITE)
title(s, "Validation after all changes",
      "Full values now; the old ρ₃-only comparison "
      "missed the vv×vv ρ₂ bug")
rows = [["Check", "Case", "Result"],
        ["Kernel = commutator (ρ₃ off/on, a≠b too)", "He8 0⁺, emax 3", "1.5e-14"],
        ["", "He8 2⁺, emax 3", "1.9e-12"],
        ["", "He8 1⁺/2⁺/3⁺ itz 1", "≤ 1.4e-11"],
        ["Scalar ρ₃ with real kshell_threebody ρ₃", "Ne20 0⁺, emax 4", "≤ 7e-11"],
        ["Metric negative modes with ρ₃", "He8 all 13 jobs, Ne20 0⁺", "0  (min ≈ −1e-10)"],
        ["Metric negative modes, ρ₂ only", "He8 / Ne20", "2–6 / 18  (expected)"],
        ["[Tz, χ] = ΔTz·χ", "Li8 1⁺ / 2⁺ configs", "0 / 1274,  0 / 1598 fail"],
        ["ρ₃ trace = C(N,3)", "all KSHELL refs", "exact"]]
table(s, 0.6, 1.85, 12.1, rows, [5.4, 3.7, 3.0], size=15, row_h=0.52,
      colors={(i, 2): NEW for i in range(1, 9) if i != 6})

# ---------------------------------------------------------------- 9 old vs new metric
s = prs.slides.add_slide(BLANK)
bg(s, WHITE)
title(s, "Old vs new: charge-exchange configuration space",
      "Old = run 214864 (ρ₃ kernel, before the B fixes)   ·   New = run 215057 (all fixes)")
rows = [["He8 → Li8", "dim old", "dim new", "neg modes ρ₃ old", "neg modes ρ₃ new", "Arnoldi ρ₃ old", "Arnoldi ρ₃ new"],
        ["1⁺  (cx1)", "2482", "1274", "4", "0", "33 steps, stalled", "71, converged"],
        ["2⁺  (cx2)", "3128", "1598", "5", "0", "20 steps, stalled", "63, converged"],
        ["3⁺  (cx3)", "2860", "1452", "3", "0", "16 steps, stalled", "99, not conv."]]
table(s, 0.6, 1.85, 12.1, rows, [1.9, 1.4, 1.4, 1.95, 1.95, 1.75, 1.75], size=15, row_h=0.6,
      colors={**{(i, 3): OLD for i in (1, 2, 3)}, **{(i, 4): NEW for i in (1, 2, 3)}})
text(s, 0.6, 4.6, 12.1, 2.2, [
    {"bullet": True, "runs": [("Dimension halves ", {"bold": True}),
                              ("because only the n→p direction is kept (B3).", {})]},
    {"bullet": True, "runs": [("Old ρ₃ metric was indefinite, ", {"bold": True}),
                              ("so the Arnoldi broke down after 16–33 steps.", {})]},
    {"bullet": True, "runs": [("Same-nucleus jobs (itz = 0): ", {"bold": True}),
                              ("tensor dim stays 2371, but null modes drop 1955 → 1660 (ρ₃): the dead "
                               "swapped configs became live (B2).", {})]}],
     size=16, space=8)

# ---------------------------------------------------------------- 10 Li8 chart
s = prs.slides.add_slide(BLANK)
bg(s, WHITE)
title(s, "Li8: overbinding gone, now underbound by 4–6 MeV",
      "Lowest Li8 state per J, E − E(He8 gs) in MeV, emax 3")
cd = CategoryChartData()
cd.categories = ["1⁺", "2⁺", "3⁺"]
series = [("Old, ρ₃", (-17.56, -34.00, -2.47), OLD),
          ("New, ρ₂ only", (-14.63, -12.37, -12.37), RHO2),
          ("New, ρ₃", (-5.04, -4.89, -3.30), NEW),
          ("KSHELL", (-10.55, -11.33, -7.56), REF)]
for name, vals, _ in series:
    cd.add_series(name, vals)
gf = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(0.6), Inches(1.6), Inches(7.9), Inches(5.3), cd)
ch = gf.chart
ch.has_legend = True
ch.legend.position = XL_LEGEND_POSITION.BOTTOM
ch.legend.include_in_layout = False
ch.legend.font.size, ch.legend.font.name = Pt(13), BODY
ch.font.name = BODY
for (name, vals, col), ser in zip(series, ch.series):
    ser.format.fill.solid()
    ser.format.fill.fore_color.rgb = col
    ser.invert_if_negative = False
    ser.data_labels.show_value = True
    ser.data_labels.number_format = '0.0'
    ser.data_labels.number_format_is_linked = False
    ser.data_labels.position = XL_LABEL_POSITION.OUTSIDE_END
    ser.data_labels.font.size = Pt(11)
    ser.data_labels.font.color.rgb = INK2
va = ch.value_axis
va.maximum_scale, va.minimum_scale = 0, -36
va.major_gridlines.format.line.color.rgb = RGBColor(0xDD, 0xE2, 0xEA)
va.tick_labels.font.size = Pt(12)
va.tick_labels.font.color.rgb = INK2
va.format.line.fill.background()
ca = ch.category_axis
ca.tick_labels.font.size = Pt(14)
ca.tick_labels.font.color.rgb = INK2
from pptx.enum.chart import XL_TICK_LABEL_POSITION
ca.tick_label_position = XL_TICK_LABEL_POSITION.HIGH
text(s, 8.85, 1.75, 3.9, 5.1, [
    [("Old ρ₃", {"bold": True, "color": OLD})],
    "Indefinite metric: −34 MeV 2⁺ state.",
    [("ρ₂ only", {"bold": True, "color": RHO2})],
    "Same −12.37 / −12.00 roots for every J — artefacts of negative modes.",
    [("New ρ₃", {"bold": True, "color": NEW})],
    "Artefacts gone; 4–6 MeV above KSHELL; 1⁺ below 2⁺ by 0.16 (KSHELL: 2⁺ lowest by 0.78).",
    [("Next roots, new ρ₃: ", {"bold": True}), ("1⁺ −3.15, −1.81 (KS −7.18, −5.84); 2⁺ −1.56 (KS −6.48)", {})]],
     size=14, space=5)

# ---------------------------------------------------------------- 11 He8 same nucleus
s = prs.slides.add_slide(BLANK)
bg(s, WHITE)
title(s, "He8 same-nucleus excitations: old vs new",
      "Lowest root, E − E(He8 gs) in MeV. Old: 59 Arnoldi steps; new: 99. None converged.")
cd = CategoryChartData()
cd.categories = ["0⁺ → 0⁺₂  (mr0gs)", "0⁺ → 2⁺  (t02)", "2⁺₁ → 2⁺  (mr2gs)"]
series = [("Old, ρ₂", (5.92, 10.34, 11.25), RGBColor(0xE8, 0xA8, 0x9A)),
          ("Old, ρ₃", (6.40, 13.00, 10.34), OLD),
          ("New, ρ₂", (6.11, 8.54, 9.20), RHO2),
          ("New, ρ₃", (5.76, 7.16, 9.20), NEW),
          ("KSHELL", (9.71, 5.27, 11.28), REF)]
for name, vals, _ in series:
    cd.add_series(name, vals)
gf = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(0.6), Inches(1.6), Inches(7.9), Inches(5.3), cd)
ch = gf.chart
ch.has_legend = True
ch.legend.position = XL_LEGEND_POSITION.BOTTOM
ch.legend.include_in_layout = False
ch.legend.font.size, ch.legend.font.name = Pt(13), BODY
ch.font.name = BODY
for (name, vals, col), ser in zip(series, ch.series):
    ser.format.fill.solid()
    ser.format.fill.fore_color.rgb = col
    ser.data_labels.show_value = True
    ser.data_labels.number_format = '0.0'
    ser.data_labels.number_format_is_linked = False
    ser.data_labels.position = XL_LABEL_POSITION.OUTSIDE_END
    ser.data_labels.font.size = Pt(11)
    ser.data_labels.font.color.rgb = INK2
va = ch.value_axis
va.minimum_scale, va.maximum_scale = 0, 14
va.major_gridlines.format.line.color.rgb = RGBColor(0xDD, 0xE2, 0xEA)
va.tick_labels.font.size = Pt(12)
va.tick_labels.font.color.rgb = INK2
va.format.line.fill.background()
ch.category_axis.tick_labels.font.size = Pt(13)
ch.category_axis.tick_labels.font.color.rgb = INK2
text(s, 8.85, 1.75, 3.9, 5.1, [
    [("t02 (tensor): 13.0 → 7.2", {"bold": True, "color": NEW})],
    "Biggest change: the B2 ladder fix made the dropped pp(2)←hh(0) configs live. KSHELL 2⁺₁ is 5.27.",
    [("Scalar jobs", {"bold": True})],
    "Only the projector changed (B6); differences of ≤ 1 MeV are mostly the longer Arnoldi run.",
    [("Caveat", {"bold": True, "color": OLD})],
    "No itz = 0 solve converged; these lowest roots can still move. A full-rank run "
    "(≈235 scalar, ≈700–770 tensor steps) settles it."], size=14, space=5)

# ---------------------------------------------------------------- 12 Ne20
s = prs.slides.add_slide(BLANK)
bg(s, WHITE)
title(s, "Ne20 0⁺ → 0⁺: ρ₃ removes spurious low states",
      "Scalar MR-EOM on the KSHELL ground state, 59 Arnoldi steps (not converged), dim 5954")
rows = [["", "ρ₂ only", "ρ₃", "KSHELL (in-valence)"],
        ["Negative metric modes", "18", "0", "—"],
        ["E_x(0)  [MeV]", "8.05", "8.60", "0⁺₂  9.66"],
        ["E_x(1)", "11.22", "23.51", ""],
        ["E_x(2)", "11.69", "25.93", ""],
        ["E_x(3)", "13.06", "26.65", ""],
        ["E_x(4)", "14.14", "27.74", ""]]
table(s, 0.6, 1.85, 7.6, rows, [2.6, 1.5, 1.5, 2.0], size=15, row_h=0.55,
      colors={(1, 1): OLD, (1, 2): NEW, (2, 2): NEW})
card(s, 8.6, 1.85, 4.1, 3.85)
text(s, 8.85, 2.0, 3.6, 3.6, [
    [("Reading", {"bold": True, "size": 18})],
    "The ρ₂-only cluster at 11–14 MeV sits on the 18 negative-norm modes; with ρ₃ it disappears.",
    "The lowest 0⁺ moves 0.55 MeV toward KSHELL.",
    [("Ran before the B fixes; for scalar runs those change only the projector.", {"italic": True})]],
     size=15, space=8)

# ---------------------------------------------------------------- 13 open
s = prs.slides.add_slide(BLANK)
bg(s, INK)
title(s, "Open items", dark=True)
items = [("1", "Converge the Arnoldi",
          "Krylov built with H2, Ritz values with H2 + H3 → slow. Run to full rank: ≈235 steps scalar, "
          "≈700–770 tensor, ~1 s each."),
         ("2", "vv←vv H side (your call)",
          "Li8 is 4–6 MeV underbound; the audit points to missing ρ₃/ρ₄-level H terms in the vv←vv transfer."),
         ("3", "Rerun Mg24 E1",
          "The old run used the ladder that dropped swapped configs (B2)."),
         ("4", "H = T² test",
          "Gives −3, −1.5 even for scalar 0⁺; not investigated."),
         ("5", "Paper / docs",
          "Deferred, as agreed. Nothing committed yet.")]
for k, (n, head, body) in enumerate(items):
    y = 1.45 + k * 1.1
    badge(s, 0.7, y + 0.05, n, fill=AMBER, d=0.6, size=18)
    text(s, 1.55, y, 11.0, 0.45, [[(head, {"bold": True})]], size=19, color=WHITE)
    text(s, 1.55, y + 0.45, 11.0, 0.6, body, size=15, color=RGBColor(0xC8, 0xD2, 0xE0))

out = "mreom_rho3_fixes.pptx"
prs.save(out)
print("wrote", out)
