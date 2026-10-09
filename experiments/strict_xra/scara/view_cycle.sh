#!/bin/zsh
# Local replay window of the SCARA kitting cycle (macOS). Uses the conda otool because the system one needs the Xcode licence.
# Usage: ./view_scara.sh [speed]   (default 1.0 = real time)
HERE=${0:A:h}; ENV=/Users/jianqiang/miniforge3/envs/rab_franka
mkdir -p "$HERE/build/toolbin"; [ -e "$HERE/build/toolbin/otool" ] || ln -sf "$ENV/bin/arm64-apple-darwin20.0.0-otool" "$HERE/build/toolbin/otool"
PATH="$HERE/build/toolbin:$PATH" exec "$ENV/bin/mjpython" "$HERE/view_cycle.py" "${1:-1.0}"
