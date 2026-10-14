#!/bin/bash
# GIFs of the MuJoCo videos for the README (ffmpeg, two-pass palette).
# Usage: docs/media/make_gifs.sh   (from anywhere; writes next to this script)
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"; X="$HERE/../../experiments/strict_xra"; PAL="$(mktemp -t pal).png"
gif() {  # input, output, speed-up, fps, width
    local f="setpts=PTS/$3,fps=$4,scale=$5:-1:flags=lanczos"
    ffmpeg -v error -y -i "$1" -vf "$f,palettegen=stats_mode=diff:max_colors=200" "$PAL"
    ffmpeg -v error -y -i "$1" -i "$PAL" -lavfi "$f [x]; [x][1:v] paletteuse=dither=bayer:bayer_scale=4:diff_mode=rectangle" "$2"
    echo "$(basename "$2"): $(du -h "$2" | cut -f1)"
}
gif "$X/scara/figs/scara_cycle.mp4"                "$HERE/scara_cycle.gif"     3 12 800
gif "$X/scara/figs/scara_red_task.mp4"             "$HERE/scara_red_task.gif"  2 15 800
gif "$X/franka/figs/red_fwd_sq/franka_red_task.mp4" "$HERE/franka_red_task.gif" 2 15 800
rm -f "$PAL"
