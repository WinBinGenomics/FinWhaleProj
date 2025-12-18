
import dadi

# hardcoded monomorphic total for the triple
MONO_3D = {
    "ESP-ENP-GOC": 380774224,   # from your 3D mono-sum script
}

# >>> set this when working with the 3D SFS <<<
triple = "ESP-ENP-GOC"

infile  = f"{triple}_combined.sfs"          # e.g. ENP-ESP-GOC_combined.sfs
outfile = f"{triple}_combined.withMono.sfs"

fs = dadi.Spectrum.from_file(infile)

assert fs.ndim == 3, f"Expected 3D spectrum, got {fs.ndim}D"

# Totals BEFORE
seg_before = float(fs.sum())        # dadi's segregating sites (excludes masked bins)
L_before   = float(fs.data.sum())   # total sites, including mono

# Ensure monomorphic bin is unmasked so it is written to file
fs.mask[0, 0, 0] = False

# ADD mono to the (0,0,0) cell
fs.data[0, 0, 0] = float(fs.data[0, 0, 0]) + MONO_3D[triple]

# Totals AFTER
seg_after = float(fs.sum())        # should be the same (dadi still excludes monomorphic bin)
L_after   = float(fs.data.sum())   # updated total L with mono

# Save updated spectrum
fs.to_file(outfile)

# Report
print(f"=== Processing {triple} ===")
print(f"Added {MONO_3D[triple]} to bin (0,0,0). New (0,0,0) = {fs.data[0,0,0]:.6f}")
print(f"Segregating sites (exclude mono): {int(round(seg_before))} → {int(round(seg_after))}")
print(f"Total sites L (include mono):     {int(round(L_before))} → {int(round(L_after))}")
print(f"Wrote {outfile}")
