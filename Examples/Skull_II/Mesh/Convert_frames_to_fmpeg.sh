#!/bin/bash

# This script converts frames to a video file using ffmpeg.
# Usage: ./Convert_frames_to_fmpeg.sh frame_pattern

# Check if an argument is given
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 frame_pattern"
    exit 1
fi

# The pattern for frame filenames, e.g., "trui_frame_"
FRAME_PATTERN=$1

# Command to convert frames to video
ffmpeg -framerate 5 -i "${FRAME_PATTERN}%02d.png" -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -c:v libx264 -profile:v high -crf 35 -pix_fmt yuv420p -r 30 out.mp4

