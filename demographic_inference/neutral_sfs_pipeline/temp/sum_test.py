import dadi

pairs   = ["ENP-GOC", "ENP-ESP", "ESP-GOC"]
nchrom  = 21

# hardcoded monomorphic counts
mono_map = {
    "ENP-GOC": 383581853,
    "ENP-ESP": 381921730,
    "ESP-GOC": 383786069
}

for pair in pairs:
    prefix  = f"{pair}_chr"
    outfile = f"{pair}_combined.sfs"
    combined = None

    # ----- combine per-chrom spectra -----
    for i in range(1, nchrom + 1):
        fname = f"{prefix}{i:02d}.sfs"
        try:
            S = dadi.Spectrum.from_file(fname)
        except FileNotFoundError:
            print(f"{pair}: Skipping missing {fname}")
            continue
        combined = S if combined is None else (combined + S)
        print(f"{pair}: Added {fname}")

    if combined is None:
        print(f"{pair}: No spectra loaded, skipping this pair.")
        continue

    combined.to_file(outfile)
    print(f"{pair}: Final combined SFS written to {outfile}")
    print(f"{pair}: Projection sizes = {combined.sample_sizes}, Folded = {combined.folded}")

    # ----- add monomorphic to (0,0) -----
    mono_count = mono_map.get(pair, 0)
    before = int(combined.data.sum())
    combined[0, 0] += mono_count
    after = int(combined.data.sum())

    out_with_mono = f"{pair}_combined.withMono.sfs"
    combined.to_file(out_with_mono)

    print(f"{pair}: Added {mono_count} monomorphic sites to (0,0)")
    print(f"{pair}: Total sites before {before}, after {after}")
    print(f"{pair}: Monomorphic-added SFS written to {out_with_mono}\n")
