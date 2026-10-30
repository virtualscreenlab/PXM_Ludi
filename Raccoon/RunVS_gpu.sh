#!/bin/bash
ADGPU=/media/benq/hdd/Shityakoff/Software/AutoDock-GPU/bin/autodock_gpu_128wi

count=0
skip=0
fail=0

for dir in pxm_*_receptor.protein ludi_*_receptor.protein; do
    [ -d "$dir" ] || continue

    gpf=$(find "$dir" -maxdepth 1 -name '*.gpf' | head -1)
    lig=$(find "$dir" -maxdepth 1 -name '*.pdbqt' ! -name '*receptor*' | head -1)

    if [ -z "$gpf" ] || [ -z "$lig" ]; then
        echo "SKIP $dir"
        skip=$((skip+1))
        continue
    fi

    echo "[$((count+1))] $dir"
    (
        cd "$dir" || exit 1

        if ! ls *.maps.fld >/dev/null 2>&1; then
            autogrid4 -p "$(basename "$gpf")" -l autogrid.log >/dev/null 2>&1
        fi

        fld=$(ls *.maps.fld 2>/dev/null | head -1)
        if [ -z "$fld" ]; then
            echo "   FAIL: no .maps.fld"
            exit 2
        fi

        "$ADGPU" --ffile "$fld" --lfile "$(basename "$lig")" --nrun 20 >/dev/null
    )
    rc=$?
    if [ $rc -eq 0 ]; then count=$((count+1)); else fail=$((fail+1)); fi
done

echo "==========================================="
echo "Done: $count docked, $skip skipped, $fail failed."
