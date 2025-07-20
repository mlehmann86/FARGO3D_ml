


cd /tiara/home/mlehmann/data/FARGO3D/outputs/

process_id=$!
wait $process_id

module add "idl/8.7.3"

process_id=$!
wait $process_id

idl


