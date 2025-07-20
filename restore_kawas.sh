#!/bin/bash

SRC="/tiara/home/mlehmann/data/FARGO3D/"
DEST="mlehmann@kawas:/scratch/data/mlehmann/FARGO3D/"

echo "⚠️  This will overwrite all differing files, even if timestamps are the same."
echo "    From: $SRC"
echo "    To:   $DEST"
echo

# Confirm
read -p "Proceed with full overwrite? (yes/no): " confirm
if [[ "$confirm" != "yes" ]]; then
    echo "Aborted."
    exit 1
fi

# Rsync with full overwrite (by checksum), plus excludes
rsync -avz --checksum --no-times \
    --exclude='sedwckIhg' \
    --exclude='outputs/' \
    --exclude='Makefile' \
    --exclude='src/makefile' \
    --exclude='*.out' \
    --exclude='*.dat' \
    --exclude='fargo3d_*' \
    "$SRC" "$DEST"
