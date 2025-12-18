


import dadi

pair    = "ESP-ENP-GOC"           # adjust as needed., "ENP-ESP" or "ENP-GOC"
nchrom  = 21
prefix  = f"{pair}_chr"
outfile = f"{pair}_combined.sfs"

combined = None

for i in range(1, nchrom + 1):
    fname = f"{prefix}{i:02d}.sfs"
    try:
        S = dadi.Spectrum.from_file(fname)
    except FileNotFoundError:
        print(f"Skipping missing: {fname}")
        continue

    combined = S if combined is None else (combined + S)
    print(f"Added {fname}")

if combined is None:
    raise RuntimeError("No spectra were loaded. Check your prefix/path.")

combined.to_file(outfile)
print(f"\nFinal combined SFS written to {outfile}")
print("Projection sizes:", combined.sample_sizes)
print("Folded?:", combined.folded)


