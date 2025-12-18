import dadi

def dadi_to_fsc_obs(infile: str, outfile: str, *, as_int: bool = False, unmask00: bool = True):
    """
    Convert a 2D dadi Spectrum file to fastsimcoal2 .obs format.
    - as_int=False -> keep decimals (no rounding)
    - unmask00=True -> ensure (0,0) is written (mono bin shows up)
    """
    fs = dadi.Spectrum.from_file(infile)
    if fs.ndim != 2:
        raise ValueError("This simple converter expects a 2D Spectrum.")
    if unmask00:
        fs.mask[0, 0] = False

    nrow, ncol = fs.shape
    with open(outfile, "w") as f:
        f.write("1 observation\n")
        f.write("\t" + "\t".join(f"d0_{j}" for j in range(ncol)) + "\n")
        for i in range(nrow):
            row_vals = []
            for j in range(ncol):
                v = fs.data[i, j]
                row_vals.append(str(int(round(v))) if as_int else repr(float(v)))
            f.write(f"d1_{i}\t" + "\t".join(row_vals) + "\n")

    print(f"Wrote fsc2 .obs to {outfile} (shape {nrow}x{ncol}, as_int={as_int})")

# --- usage ---
pair = "ENP-ESP"
dadi_in  = f"{pair}_combined.withMono.sfs"   # dadi file (after adding mono to (0,0))
fsc_out  = f"{pair}.fsc.obs"
dadi_to_fsc_obs(dadi_in, fsc_out, as_int=False, unmask00=True)
