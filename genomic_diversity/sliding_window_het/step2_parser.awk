#!/usr/bin/awk -f

# Project: Eastern pacific fin whale project @ CSUSM. Student -> Kaden Winspear
# April 2025
# the purpose of this script is to parse the final file of SlidingWindowhet.py. This will form files that are sample specific while keeping the correct header. This step is only necessary because the
# slidingwindowhet.py script that this project uses was originally made to fire on one VCF containing all the necessary scaffolds. However, this projects VCFs, are split by chrom / scaffold.

BEGIN {
    FS = OFS = "\t"
}

FNR == 1 {
    # Extracts the file number from FILENAME if not already done.
    # This regex looks for _VariantFiltration_([0-9]+)_ in the file name.
    if (file_num == "" && match(FILENAME, /_VariantFiltration_([0-9]+)_/, arr)) {
        file_num = arr[1]
    } else if (file_num == "") {
        file_num = "NA"  # in case the pattern isn’t found.
    }

    # Process header: get column indices for each sample.
    for (i = 1; i <= NF; i++) {
        if ($i ~ /^calls_/) {
            sample = substr($i, 7)  # Remove "calls_" prefix to get sample ID.
            calls_idx[sample] = i
        } else if ($i ~ /^hets_/) {
            sample = substr($i, 6)  # Remove "hets_" prefix to get sample ID.
            hets_idx[sample] = i
        }
    }
    next
}

# Processes data lines -> swag.
{
    for (sample in calls_idx) {
        # Creates output file name: sampleID_fileNum.txt (e.g., ENPAK19_01_temp.txt)
        file = sample "_" file_num "_temp.txt"

        # Write header to the file only once so there is no overwriting of the header which would lead to errors when calculating / plotting.
        if (!(file in written_headers)) {
            if (!system("[ -s " file " ]")) {
                # If the file exists and is non-empty, assume the header is already present, again to avoid error or rewriting header/
                written_headers[file] = 1
            } else {
                print "chrom\twindow_start\tsites_total\tcalls_" sample "\thets_" sample > file
                written_headers[file] = 1
            }
        }

        # Output the desired columns.
        print $1, $2, $3, $(calls_idx[sample]), $(hets_idx[sample]) >> file
    }
}

# lets go!
