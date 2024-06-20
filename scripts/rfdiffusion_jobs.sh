# 1. tr_r366_381_fixed
python run_rfdiffusion.py \
--job_name=tr_r366_381_fixed \
--contigs="[A6-39/A92-178/A313-342/A429-439/A543-564/A986-1011/A1102-1113/A1207-1231/0 B33-365/15-15/B381-600]" \
--hotspot_residues="[A318,A320,A322,A324]" \
--logs_dir=/home/path/to/logs

# 2. tr_r366_381_flexible
python run_rfdiffusion.py \
--job_name=tr_r366_381_flexible \
--contigs="[A6-39/A92-178/A313-342/A429-439/A543-564/A986-1011/A1102-1113/A1207-1231/0 B33-365/5-25/B381-600]" \
--hotspot_residues="[A318,A320,A322,A324]" \
--logs_dir=/home/path/to/logs

# 3. tr_r366_448_fixed
python run_rfdiffusion.py \
--job_name=tr_r366_448_fixed \
--contigs="[A6-39/A92-178/A313-342/A429-439/A543-564/A986-1011/A1102-1113/A1207-1231/0 B33-366/83-83/B449-600]" \
--hotspot_residues="[A318,A320,A322,A324]" \
--logs_dir=/home/path/to/logs

# 4. tr_r366_448_flexible
python run_rfdiffusion.py \
--job_name=tr_r366_448_flexible \
--contigs="[A6-39/A92-178/A313-342/A429-439/A543-564/A986-1011/A1102-1113/A1207-1231/0 B33-366/40-100/B449-600]" \
--hotspot_residues="[A318,A320,A322,A324]" \
--logs_dir=/home/path/to/logs

# 5. tr_r437_448_fixed
python run_rfdiffusion.py \
--job_name=tr_r437_448_fixed \
--contigs="[A6-39/A92-178/A313-342/A429-439/A543-564/A986-1011/A1102-1113/A1207-1231/0 B33-436/12-12/B449-600]" \
--hotspot_residues="[A318,A320,A322,A324]" \
--logs_dir=/home/path/to/logs

# 6. tr_r437_448_flexible
python run_rfdiffusion.py \
--job_name=tr_r437_448_flexible \
--contigs="[A6-39/A92-178/A313-342/A429-439/A543-564/A986-1011/A1102-1113/A1207-1231/0 B33-436/5-18/B449-600]" \
--hotspot_residues="[A318,A320,A322,A324]" \
--logs_dir=/home/path/to/logs
