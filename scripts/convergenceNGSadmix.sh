#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  scripts/convergenceNGSadmix.sh <beagle_file> <max_seed> <threads> <output_dir> <K> <conv_times> <ngsadmix_root> [start_seed] [ll_window] [q_threshold]

Arguments:
  beagle_file     Input Beagle likelihood file, plain text or .gz
  max_seed        Last seed to test, inclusive
  threads         Number of NGSadmix threads per run
  output_dir      Directory for all outputs
  K               Number of ancestral populations
  conv_times      Required number of converged runs
  ngsadmix_root   Repository root containing ./NGSadmix

Optional:
  start_seed      First seed to test, default: 1
  ll_window       Log-likelihood window for LL-based convergence, default: 3
  q_threshold     Max element-wise Q difference for Q-based convergence, default: 0.01

Exit codes:
  2   Invalid command-line usage
  3   Missing dependency
  4   Missing input file or binary
  5   Invalid numeric argument
  6   Invalid Beagle input format
  7   Failed to create output directory or files
  20  NGSadmix run failed
  21  Could not parse best likelihood from a run log
  22  Expected NGSadmix output file missing
  30  Q convergence evaluation failed
  40  No convergence reached within the requested seeds
EOF
}

log() {
    printf '[INFO] %s\n' "$*"
}

warn() {
    printf '[WARN] %s\n' "$*" >&2
}

die() {
    local code=$1
    shift
    printf '[ERROR] %s\n' "$*" >&2
    exit "$code"
}

is_pos_int() {
    [[ ${1:-} =~ ^[1-9][0-9]*$ ]]
}

is_nonneg_number() {
    [[ ${1:-} =~ ^[0-9]+([.][0-9]+)?$ ]]
}

read_beagle() {
    local path=$1
    if [[ $path == *.gz ]]; then
        gzip -cd -- "$path"
    else
        cat -- "$path"
    fi
}

validate_beagle() {
    local path=$1
    local awk_status=0

    set +e
    read_beagle "$path" | awk '
        NR == 1 {
            if (NF < 6) {
                exit 11
            }
            if (((NF - 3) % 3) != 0) {
                exit 12
            }
            next
        }
        NR == 2 {
            found_data = 1
        }
        END {
            if (NR == 0) {
                exit 13
            }
            if (!found_data) {
                exit 14
            }
        }
    '
    awk_status=$?
    set -e

    case $awk_status in
        0) ;;
        11) die 6 "Beagle header has too few columns: $path" ;;
        12) die 6 "Beagle header does not match 3 columns per individual: $path" ;;
        13) die 6 "Beagle input is empty: $path" ;;
        14) die 6 "Beagle input has no data rows: $path" ;;
        *) die 6 "Could not validate Beagle input: $path" ;;
    esac
}

extract_best_like() {
    local log_file=$1
    awk '
        /^best like=/ {
            sub("like=", "", $2)
            print $2
            found = 1
            exit
        }
        END {
            if (!found) {
                exit 1
            }
        }
    ' "$log_file"
}

ll_conv_count() {
    local likes_file=$1
    local window=$2
    awk -v window="$window" '
        NR == 1 {
            max = $2
        }
        {
            values[NR] = $2
            if ($2 > max) {
                max = $2
            }
        }
        END {
            count = 0
            for (i = 1; i <= NR; i++) {
                if ((max - values[i]) < window) {
                    count++
                }
            }
            print count + 0
        }
    ' "$likes_file"
}

if [[ $# -lt 7 || $# -gt 10 ]]; then
    usage >&2
    exit 2
fi

file=$1
max_seed=$2
threads=$3
out_dir=$4
K=$5
conv_times=$6
ngsadmix_root=$7
start_seed=${8:-1}
ll_window=${9:-3}
q_threshold=${10:-0.01}

is_pos_int "$max_seed" || die 5 "max_seed must be a positive integer"
is_pos_int "$threads" || die 5 "threads must be a positive integer"
is_pos_int "$K" || die 5 "K must be a positive integer"
is_pos_int "$conv_times" || die 5 "conv_times must be a positive integer"
is_pos_int "$start_seed" || die 5 "start_seed must be a positive integer"
is_nonneg_number "$ll_window" || die 5 "ll_window must be a non-negative number"
is_nonneg_number "$q_threshold" || die 5 "q_threshold must be a non-negative number"

(( start_seed <= max_seed )) || die 5 "start_seed must be less than or equal to max_seed"

[[ -r $file ]] || die 4 "Input file not found or not readable: $file"
[[ -d $ngsadmix_root ]] || die 4 "NGSadmix root directory not found: $ngsadmix_root"

adm_bin=$ngsadmix_root/NGSadmix
qconv_script=$ngsadmix_root/scripts/Qconv.R

[[ -x $adm_bin ]] || die 4 "NGSadmix binary not found or not executable: $adm_bin"
[[ -r $qconv_script ]] || die 4 "Q convergence script not found: $qconv_script"
command -v Rscript >/dev/null 2>&1 || die 3 "Rscript not found in PATH"
command -v awk >/dev/null 2>&1 || die 3 "awk not found in PATH"
command -v sort >/dev/null 2>&1 || die 3 "sort not found in PATH"

validate_beagle "$file"

mkdir -p -- "$out_dir" || die 7 "Could not create output directory: $out_dir"

bfile=$(basename -- "$file")
prefix=$out_dir/$bfile.$K
likes_tmp=$prefix.likes.tmp
likes_sorted=$prefix.likes
qlist=$out_dir/$K.Qlist
status_tsv=$prefix.run_status.tsv
summary_file=$prefix.summary.txt

: > "$likes_tmp" || die 7 "Could not write $likes_tmp"
: > "$qlist" || die 7 "Could not write $qlist"
printf 'seed\tngsadmix_exit\tbest_like\tll_conv_count\tq_conv_count\tresult\n' > "$status_tsv" || die 7 "Could not write $status_tsv"

log "file = $file"
log "max_seed = $max_seed"
log "threads = $threads"
log "output_dir = $out_dir"
log "K = $K"
log "conv_times = $conv_times"
log "start_seed = $start_seed"
log "ll_window = $ll_window"
log "q_threshold = $q_threshold"

converged=0
best_seed=
best_like=
best_q=
best_f=
best_log=

for seed in $(seq "$start_seed" "$max_seed"); do
    run_log=$prefix.log_$seed
    log "Running seed $seed"

    if ! "$adm_bin" -likes "$file" -K "$K" -seed "$seed" -P "$threads" -outfiles "$prefix" 2> "$run_log"; then
        run_code=$?
        printf '%s\t%s\tNA\tNA\tNA\tngsadmix_failed\n' "$seed" "$run_code" >> "$status_tsv"
        die 20 "NGSadmix failed for seed $seed with exit code $run_code. See $run_log"
    fi

    run_like=$(extract_best_like "$run_log") || {
        printf '%s\t0\tNA\tNA\tNA\tmissing_best_like\n' "$seed" >> "$status_tsv"
        die 21 "Could not parse best likelihood from $run_log"
    }

    qopt_file=$prefix.qopt
    fopt_file=$prefix.fopt.gz

    [[ -f $qopt_file ]] || die 22 "Missing expected Q output: $qopt_file"
    [[ -f $fopt_file ]] || die 22 "Missing expected F output: $fopt_file"

    qopt_seed=$prefix.qopt.$seed
    fopt_seed=$prefix.fopt.gz.$seed
    cp -- "$qopt_file" "$qopt_seed" || die 7 "Could not save $qopt_seed"
    cp -- "$fopt_file" "$fopt_seed" || die 7 "Could not save $fopt_seed"

    printf '%s\t%s\n' "$seed" "$run_like" >> "$likes_tmp"
    printf '%s\n' "$qopt_seed" >> "$qlist"

    if [[ -z ${best_like:-} ]] || awk -v a="$run_like" -v b="$best_like" 'BEGIN { exit !(a > b) }'; then
        best_seed=$seed
        best_like=$run_like
        best_q=$qopt_seed
        best_f=$fopt_seed
        best_log=$run_log
    fi

    ll_count=$(ll_conv_count "$likes_tmp" "$ll_window")
    q_count=$(Rscript "$qconv_script" "$likes_tmp" "$qlist" "$q_threshold" --count-only) || {
        printf '%s\t0\t%s\t%s\tNA\tqconv_failed\n' "$seed" "$run_like" "$ll_count" >> "$status_tsv"
        die 30 "Q convergence evaluation failed for seed $seed"
    }

    result=continue
    if (( ll_count >= conv_times )) && (( q_count >= conv_times )); then
        result=converged_both
    elif (( ll_count >= conv_times )); then
        result=converged_ll
    elif (( q_count >= conv_times )); then
        result=converged_q
    fi

    printf '%s\t0\t%s\t%s\t%s\t%s\n' "$seed" "$run_like" "$ll_count" "$q_count" "$result" >> "$status_tsv"
    log "seed $seed: best_like=$run_like ll_conv_count=$ll_count q_conv_count=$q_count result=$result"

    if [[ $result != continue ]]; then
        cp -- "$qopt_seed" "$prefix.qopt_conv" || die 7 "Could not write $prefix.qopt_conv"
        cp -- "$fopt_seed" "$prefix.fopt_conv.gz" || die 7 "Could not write $prefix.fopt_conv.gz"
        cp -- "$run_log" "$prefix.log_conv" || die 7 "Could not write $prefix.log_conv"
        converged=1
        break
    fi
done

sort -k2,2nr "$likes_tmp" > "$likes_sorted" || die 7 "Could not write $likes_sorted"

if [[ -n ${best_q:-} ]]; then
    cp -- "$best_q" "$prefix.qopt_best" || die 7 "Could not write $prefix.qopt_best"
    cp -- "$best_f" "$prefix.fopt_best.gz" || die 7 "Could not write $prefix.fopt_best.gz"
    cp -- "$best_log" "$prefix.log_best" || die 7 "Could not write $prefix.log_best"
fi

final_seed_tested=$((start_seed - 1))
if [[ -n ${seed:-} ]]; then
    final_seed_tested=$seed
fi

{
    printf 'file\t%s\n' "$file"
    printf 'K\t%s\n' "$K"
    printf 'seeds_tested\t%s-%s\n' "$start_seed" "$final_seed_tested"
    printf 'conv_times\t%s\n' "$conv_times"
    printf 'll_window\t%s\n' "$ll_window"
    printf 'q_threshold\t%s\n' "$q_threshold"
    printf 'best_seed\t%s\n' "${best_seed:-NA}"
    printf 'best_like\t%s\n' "${best_like:-NA}"
    if (( converged == 1 )); then
        printf 'converged\tyes\n'
        printf 'converged_seed\t%s\n' "$seed"
    else
        printf 'converged\tno\n'
        printf 'converged_seed\tNA\n'
    fi
} > "$summary_file" || die 7 "Could not write $summary_file"

if (( converged == 1 )); then
    log "Converged after seed $seed"
    log "Summary written to $summary_file"
    exit 0
fi

warn "No convergence reached by seed $max_seed"
warn "Best seed was ${best_seed:-NA} with likelihood ${best_like:-NA}"
exit 40
