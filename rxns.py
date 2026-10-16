#!/usr/bin/env python3
"""
List reactions (and daughter nuclides) for a starting isotope, using locally
downloaded ENDF-6 sublibraries (e.g. proton and alpha sublibraries from
ENDF/B-VIII.0 or TENDL-2023).

Usage examples:
    python rxns.py Fe56 --lib protons.zip alphas.zip
    python rxns.py U235 --lib C:\\data\\TENDL\\proton C:\\data\\TENDL\\alpha
    python rxns.py Al27 --lib p.zip a.zip --no-levels --csv al27.csv

--lib accepts any number of zip files or folders (searched recursively).
The projectile of each file (proton, alpha, neutron, ...) is detected from the
NSUB number in the file's own header, so the folder/zip names do not matter.

Two kinds of daughter information are reported:
  * MF=3 reactions: the daughter is DERIVED from the MT number using
    conservation of charge and mass (target + projectile - emitted particles).
  * MF=10 sections: residual nuclides listed EXPLICITLY in the file (if the
    evaluation provides them), with the excited/isomeric state number (LFS).

ENDF MT naming note: MT 103 is "(z,p)" = a proton is emitted. For an incident
proton that is inelastic proton scattering, so the daughter equals the target.
Only the standard library is needed.
"""
import argparse
import csv
import os
import re
import sys
import zipfile
from collections import defaultdict
from math import ceil

ELEMENTS = (
    "n H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni "
    "Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe "
    "Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg "
    "Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg "
    "Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()

# NSUB (sublibrary id in the file header): (symbol, Z, A, name)
PROJECTILES = {
    10: ("n", 0, 1, "neutron"),
    10010: ("p", 1, 1, "proton"),
    10020: ("d", 1, 2, "deuteron"),
    10030: ("t", 1, 3, "triton"),
    20030: ("3He", 2, 3, "helion"),
    20040: ("a", 2, 4, "alpha"),
}

MT_NAMES = {
    1: "(n,total)", 2: "(n,elastic)", 3: "(n,nonelastic)", 4: "(n,inelastic)",
    5: "(n,anything)", 16: "(n,2n)", 17: "(n,3n)", 18: "(n,fission)",
    19: "(n,f) first-chance", 20: "(n,nf) second-chance", 21: "(n,2nf) third-chance",
    22: "(n,na)", 23: "(n,n3a)", 24: "(n,2na)", 25: "(n,3na)", 28: "(n,np)",
    29: "(n,n2a)", 30: "(n,2n2a)", 32: "(n,nd)", 33: "(n,nt)", 34: "(n,n3He)",
    35: "(n,nd2a)", 36: "(n,nt2a)", 37: "(n,4n)", 38: "(n,3nf)", 41: "(n,2np)",
    42: "(n,3np)", 44: "(n,n2p)", 45: "(n,npa)", 91: "(n,n') continuum",
    101: "(n,disappearance)", 102: "(n,gamma)", 103: "(n,p)", 104: "(n,d)",
    105: "(n,t)", 106: "(n,3He)", 107: "(n,alpha)", 108: "(n,2a)", 109: "(n,3a)",
    111: "(n,2p)", 112: "(n,pa)", 113: "(n,t2a)", 114: "(n,d2a)", 115: "(n,pd)",
    116: "(n,pt)", 117: "(n,da)",
    649: "(n,p) continuum", 699: "(n,d) continuum", 749: "(n,t) continuum",
    799: "(n,3He) continuum", 849: "(n,a) continuum",
}
MT_RANGES = [  # (first, last, template, offset)
    (51, 90, "(n,n'{i}) to level {i}", 50),
    (600, 648, "(n,p{i}) to level {i}", 600),
    (650, 698, "(n,d{i}) to level {i}", 650),
    (700, 748, "(n,t{i}) to level {i}", 700),
    (750, 798, "(n,3He{i}) to level {i}", 750),
    (800, 848, "(n,a{i}) to level {i}", 800),
]
FISSION_MTS = {18, 19, 20, 21, 38}

# Emitted particles per MT: (neutrons, protons, deuterons, tritons, 3He, alphas)
EMIT = {
    102: (0, 0, 0, 0, 0, 0), 16: (2, 0, 0, 0, 0, 0), 17: (3, 0, 0, 0, 0, 0),
    37: (4, 0, 0, 0, 0, 0), 22: (1, 0, 0, 0, 0, 1), 23: (1, 0, 0, 0, 0, 3),
    24: (2, 0, 0, 0, 0, 1), 25: (3, 0, 0, 0, 0, 1), 28: (1, 1, 0, 0, 0, 0),
    29: (1, 0, 0, 0, 0, 2), 30: (2, 0, 0, 0, 0, 2), 32: (1, 0, 1, 0, 0, 0),
    33: (1, 0, 0, 1, 0, 0), 34: (1, 0, 0, 0, 1, 0), 35: (1, 0, 1, 0, 0, 2),
    36: (1, 0, 0, 1, 0, 2), 41: (2, 1, 0, 0, 0, 0), 42: (3, 1, 0, 0, 0, 0),
    44: (1, 2, 0, 0, 0, 0), 45: (1, 1, 0, 0, 0, 1), 103: (0, 1, 0, 0, 0, 0),
    104: (0, 0, 1, 0, 0, 0), 105: (0, 0, 0, 1, 0, 0), 106: (0, 0, 0, 0, 1, 0),
    107: (0, 0, 0, 0, 0, 1), 108: (0, 0, 0, 0, 0, 2), 109: (0, 0, 0, 0, 0, 3),
    111: (0, 2, 0, 0, 0, 0), 112: (0, 1, 0, 0, 0, 1), 113: (0, 0, 0, 1, 0, 2),
    114: (0, 0, 1, 0, 0, 2), 115: (0, 1, 1, 0, 0, 0), 116: (0, 1, 0, 1, 0, 0),
    117: (0, 0, 1, 0, 0, 1),
}


# ----------------------------------------------------------------- nuclides
def nuclide_name(z, a):
    if 0 <= z < len(ELEMENTS):
        return f"{ELEMENTS[z]}-{a}"
    return f"Z={z},A={a}"


def parse_target(name):
    m = re.fullmatch(r"([A-Z][a-z]?)-?(\d{1,3})", name.strip())
    if not m or m.group(1) not in ELEMENTS or m.group(1) == "n":
        sys.exit(f"Could not parse target '{name}'. Use e.g. Fe56, U235, Al27.")
    sym, a = m.group(1), int(m.group(2))
    return ELEMENTS.index(sym), sym, a


# ----------------------------------------------------------------- reactions
def emission(mt):
    if mt in EMIT:
        return EMIT[mt]
    if mt == 4 or 51 <= mt <= 91:
        return (1, 0, 0, 0, 0, 0)  # neutron emission
    if 600 <= mt <= 649: return (0, 1, 0, 0, 0, 0)
    if 650 <= mt <= 699: return (0, 0, 1, 0, 0, 0)
    if 700 <= mt <= 749: return (0, 0, 0, 1, 0, 0)
    if 750 <= mt <= 799: return (0, 0, 0, 0, 1, 0)
    if 800 <= mt <= 849: return (0, 0, 0, 0, 0, 1)
    return None


def residual(tz, ta, proj, mt):
    """Daughter nucleus = target + projectile - emitted particles."""
    if mt in FISSION_MTS:
        return "fission fragments"
    if mt == 2:
        return nuclide_name(tz, ta)  # elastic scattering
    e = emission(mt)
    if e is None:
        return "-"
    n, p, d, t, he3, al = e
    _, pz, pa, _ = proj
    a_res = ta + pa - (n + p + 2 * d + 3 * t + 3 * he3 + 4 * al)
    z_res = tz + pz - (p + d + t + 2 * he3 + 2 * al)
    return nuclide_name(z_res, a_res)


def mt_label_neutron_style(mt):
    if mt in MT_NAMES:
        return MT_NAMES[mt]
    for lo, hi, tmpl, off in MT_RANGES:
        if lo <= mt <= hi:
            return tmpl.format(i=mt - off)
    return f"MT{mt}"


def reaction_label(mt, psym):
    if psym != "n":
        if mt == 4:
            return f"({psym},n) total"
        if 51 <= mt <= 90:
            return f"({psym},n{mt - 50}) to level {mt - 50}"
        if mt == 91:
            return f"({psym},n) continuum"
    lab = mt_label_neutron_style(mt)
    if psym != "n" and lab.startswith("(n,"):
        lab = f"({psym}," + lab[3:]
    return lab


def is_level_channel(mt):
    return (51 <= mt <= 90 or 600 <= mt <= 648 or 650 <= mt <= 698
            or 700 <= mt <= 748 or 750 <= mt <= 798 or 800 <= mt <= 848)


# ----------------------------------------------------------------- ENDF I/O
def endf_float(txt):
    txt = txt.strip()
    m = re.fullmatch(r"([+-]?\d*\.?\d*)([+-]\d+)", txt)
    return float(m.group(1) + "e" + m.group(2)) if m else float(txt)


def endf_int(txt):
    return int(round(endf_float(txt)))


def parse_header(text):
    """Return (Z, A, NSUB) from the first three MF1/MT451 records, or None."""
    recs = []
    for line in text.splitlines():
        if len(line) >= 75 and line[70:72].strip() == "1" and line[72:75].strip() == "451":
            recs.append(line)
            if len(recs) == 3:
                break
    if len(recs) < 3:
        return None
    try:
        za = endf_int(recs[0][0:11])
        nsub = endf_int(recs[2][44:55])
    except ValueError:
        return None
    return za // 1000, za % 1000, nsub


def scan(text):
    """Return (sorted MF3 MT list, {MT: MF10 lines})."""
    mf3, mf10 = set(), defaultdict(list)
    for line in text.splitlines():
        if len(line) < 75:
            continue
        try:
            mf, mt = int(line[70:72]), int(line[72:75])
        except ValueError:
            continue
        if mt <= 0:
            continue
        if mf == 3:
            mf3.add(mt)
        elif mf == 10:
            mf10[mt].append(line)
    return sorted(mf3), mf10


def parse_mf10(lines):
    """Explicit residual products in one MF=10 section: [(IZAP, LFS), ...]."""
    out = []
    try:
        ns = endf_int(lines[0][44:55])
        i = 1
        for _ in range(ns):
            rec = lines[i]
            izap, lfs = endf_int(rec[22:33]), endf_int(rec[33:44])
            nr, npts = endf_int(rec[44:55]), endf_int(rec[55:66])
            out.append((izap, lfs))
            i += 1 + ceil(nr / 3) + ceil(npts / 3)
    except (ValueError, IndexError):
        pass  # return what was parsed
    return out


class Library:
    """A zip file or folder of ENDF files."""

    def __init__(self, path):
        self.path = path
        self.label = os.path.basename(os.path.normpath(path))
        if os.path.isfile(path) and zipfile.is_zipfile(path):
            self.zf = zipfile.ZipFile(path)
            self.names = [n for n in self.zf.namelist() if not n.endswith("/")]
        elif os.path.isdir(path):
            self.zf = None
            self.names = [os.path.join(r, f) for r, _, fs in os.walk(path) for f in fs]
        else:
            sys.exit(f"--lib path not found or not a zip/folder: {path}")

    def read(self, name, nbytes=None):
        if self.zf:
            with self.zf.open(name) as f:
                data = f.read(nbytes) if nbytes else f.read()
        else:
            with open(name, "rb") as f:
                data = f.read(nbytes) if nbytes else f.read()
        return data.decode("ascii", errors="replace")

    def header(self, name):
        try:
            return parse_header(self.read(name, 1024))
        except Exception:
            return None

    def find(self, z, sym, a):
        """Return [(name, nsub)] of files whose header matches target (z, a)."""
        base = lambda n: n.replace("\\", "/").rsplit("/", 1)[-1]
        pat = re.compile(rf"(^|[^A-Za-z]){sym}[-_]?0*{a}(?!\d)(?![A-Za-z])", re.I)
        matches = []
        for n in (n for n in self.names if pat.search(base(n))):
            h = self.header(n)
            if h and (h[0], h[1]) == (z, a):
                matches.append((n, h[2]))
        if matches:
            return matches
        # Fallback: file names did not match; check every file's header.
        print(f"[{self.label}] no file-name match; scanning headers of "
              f"{len(self.names)} files (slow)...", file=sys.stderr)
        for n in self.names:
            h = self.header(n)
            if h and (h[0], h[1]) == (z, a):
                matches.append((n, h[2]))
        return matches


# ----------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("target", help="Starting isotope, e.g. Fe56")
    ap.add_argument("--lib", nargs="+", required=True,
                    help="Zip files or folders holding the sublibraries")
    ap.add_argument("--all", action="store_true",
                    help="Also show rows with no defined daughter (total, nonelastic, sums)")
    ap.add_argument("--no-levels", action="store_true",
                    help="Hide discrete-level channels (MT 51-90, 600-848)")
    ap.add_argument("--csv", help="Write all results to this CSV file")
    args = ap.parse_args()

    tz, sym, ta = parse_target(args.target)
    tname = nuclide_name(tz, ta)
    csv_rows = []
    found_any = False

    for path in args.lib:
        lib = Library(path)
        matches = lib.find(tz, sym, ta)
        if not matches:
            print(f"[{lib.label}] no evaluation found for {tname}\n")
            continue
        for name, nsub in matches:
            proj = PROJECTILES.get(nsub)
            if proj is None:
                print(f"[{lib.label}] {name}: unsupported sublibrary NSUB={nsub}, skipped\n")
                continue
            found_any = True
            psym, pname = proj[0], proj[3]
            fname = name.replace("\\", "/").rsplit("/", 1)[-1]
            mf3, mf10 = scan(lib.read(name))

            print("=" * 78)
            print(f"{pname} ({psym}) + {tname}    [{lib.label} :: {fname}]")
            print("=" * 78)

            print("Reactions (MF=3), daughter derived from MT:")
            print(f"  {'MT':>4}  {'Reaction':<32} Daughter")
            hidden = 0
            for mt in mf3:
                res = residual(tz, ta, proj, mt)
                if (res == "-" and not args.all) or (args.no_levels and is_level_channel(mt)):
                    hidden += 1
                    continue
                label = reaction_label(mt, psym)
                print(f"  {mt:>4}  {label:<32} {res}")
                csv_rows.append([lib.label, fname, pname, tname, "MF3-derived",
                                 mt, label, res, ""])
            if hidden:
                print(f"  ({hidden} rows hidden: use --all / omit --no-levels to show)")

            if mf10:
                print("\nResidual nuclides listed explicitly in file (MF=10):")
                print(f"  {'MT':>4}  {'Reaction':<32} {'Residual':<12} State")
                for mt in sorted(mf10):
                    if args.no_levels and is_level_channel(mt):
                        continue
                    label = reaction_label(mt, psym)
                    for izap, lfs in parse_mf10(mf10[mt]):
                        rname = nuclide_name(izap // 1000, izap % 1000)
                        state = "ground state" if lfs == 0 else f"level/isomer {lfs}"
                        print(f"  {mt:>4}  {label:<32} {rname:<12} {state}")
                        csv_rows.append([lib.label, fname, pname, tname, "MF10-explicit",
                                         mt, label, rname, state])
            else:
                print("\n(No MF=10 residual-production sections in this file.)")
            print()

    if not found_any:
        print(f"No evaluations found for {tname} in the given libraries.")
    if args.csv and csv_rows:
        with open(args.csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["library", "file", "projectile", "target", "source",
                        "MT", "reaction", "daughter", "state"])
            w.writerows(csv_rows)
        print(f"Wrote {len(csv_rows)} rows to {args.csv}")


if __name__ == "__main__":
    main()
