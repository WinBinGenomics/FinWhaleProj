import dadi

# hardcoded monomorphic totals
MONO = {
    "ENP-GOC": 383581853,
    "ENP-ESP": 381921730,
    "ESP-GOC": 383786069,
}

# >>> set this each time you run <<<
pair = "ENP-ESP"   # change to "ENP-ESP" or "ESP-GOC"

infile  = f"{pair}_combined.sfs"
outfile = f"{pair}_combined.withMono.sfs"

fs = dadi.Spectrum.from_file(infile)

# Totals BEFORE
seg_before = float(fs.sum())       # excludes masked cells (like (0,0))
L_before   = float(fs.data.sum())  # includes everything (our L baseline)

# Ensure (0,0) is WRITTEN to file so it "shows"
fs.mask[0, 0] = False

# ADD mono to the (0,0) cell (first bin you care about)
fs.data[0, 0] = float(fs.data[0, 0]) + MONO[pair]

# Totals AFTER
seg_after = float(fs.sum())        # unchanged (still excludes (0,0) in dadi inference)
L_after   = float(fs.data.sum())   # this is the updated total counts L

# Save updated spectrum
fs.to_file(outfile)

# Report
print(f"=== Processing {pair} ===")
print(f"Added {MONO[pair]} to bin (0,0). New (0,0) = {fs.data[0,0]:.6f}")
print(f"Segregating sites (exclude mono): {int(round(seg_before))} → {int(round(seg_after))}")
print(f"Total sites L (include mono):     {int(round(L_before))} → {int(round(L_after))}")
print(f"Wrote {outfile}")
