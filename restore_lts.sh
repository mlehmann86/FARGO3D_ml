#!/bin/bash

# Define source and destination
SRC="/tiara/home/mlehmann/data/FARGO3D/"
DST="/theory/lts/mlehmann/FARGO3D/"

# Enable dry run (true/false)
dryrun=false

# Log file for errors
ERROR_LOG="error.txt"
rm -f "$ERROR_LOG"

# Compose rsync options
options="-a --no-perms --no-times --human-readable --progress"

# Check dry run flag
case "$dryrun" in
  true|yes|1)
    options="$options --dry-run"
    echo "🟡 Dry run mode enabled — no files will be copied."
    ;;
  *)
    echo "⚪ Running actual copy — files will be copied and overwritten where applicable."
    ;;
esac

# Print command summary
echo ""
echo "Command: rsync $options [...]"
echo "Logging any errors to $ERROR_LOG"
echo ""

# Run rsync and capture errors
rsync $options \
  --exclude='src/makefile' \
  --exclude='Makefile' \
  --exclude='outputs/' \
  --exclude='*.dat' \
  --exclude='*.out' \
  --exclude='sedwckIhg' \
  "$SRC" "$DST" 2> "$ERROR_LOG"

# Post-run summary
exit_code=$?
if [ $exit_code -eq 0 ]; then
  echo "✅ rsync completed successfully. No errors."
elif [ $exit_code -eq 23 ]; then
  echo "⚠️  rsync completed with some skipped or failed files (code 23)."
  if [ -s "$ERROR_LOG" ]; then
    echo "📝 See the last 20 lines of $ERROR_LOG for details:"
    tail -n 20 "$ERROR_LOG"
  else
    echo "ℹ️  No specific error output captured."
  fi
else
  echo "❌ rsync failed with exit code $exit_code."
  echo "📝 Check full log in $ERROR_LOG"
fi
