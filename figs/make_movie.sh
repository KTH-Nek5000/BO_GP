#!/bin/bash
set -eu
###########
#figure numbers should be 2 digits
###########

# input check (obsolete)
#if [ $# -ne 1 ]; then
#    echo "give 1 augment: beta or gp1D" 1>&2
#    exit 1
#fi

outputFileName="opt_hist"
figs=("beta" "gp2D" "yTop")

for fig in ${figs[@]}; do
    files="$fig*.pdf" # take list of pdf files
    # pdf to png
    for file in $files; do
	filename=`basename $file .pdf` # delete extension, .pdf
	echo "convert $filename.pdf to $filename.png"
	pdftoppm -png $file png/$filename # pdf to png
	wait
	mv png/${filename}-1.png png/$filename.png # rename
    done
    echo "############### make $fig.avi from .png files###############"
    ffmpeg -y -framerate 0.5 -i png/${fig}_%02d.png png/tmp.avi # overwrite, 1 pic/sec
    wait
    echo "############# rescaling#################"
    ffmpeg -y -i png/tmp.avi -vf scale=800:600 png/$fig.avi # rescale
    rm png/tmp.avi
done

# 2 figures
#echo "################### png/${figs[0]}.avi + png/${figs[1]}.avi -> png/$outputFileName.avi ####################"
#ffmpeg -y -i png/${figs[0]}.avi -i png/${figs[1]}.avi -filter_complex "hstack" png/$outputFileName.avi
# 3 figures
echo "################### png/${figs[0]}.avi + png/${figs[1]}.avi + png/${figs[2]}.avi -> png/$outputFileName.avi ####################"
ffmpeg -y -i png/${figs[0]}.avi -i png/${figs[1]}.avi -i png/${figs[2]}.avi -filter_complex "[0:v][1:v][2:v]hstack=inputs=3[v]" -map "[v]" png/$outputFileName.avi
# obsolete
#ffmpeg -y -i png/${figs[0]}.avi -i png/${figs[1]}.avi -filter_complex "nullsrc=size=1600x600 [base]; [0:v] setpts=PTS-STARTPTS, scale=800x600 [left]; [1:v] setpts=PTS-STARTPTS, scale=800x600 [right]; [base][left] overlay=shortest=1 [tmp1]; [tmp1][right] overlay=shortest=1:x=800" -c:v libx264 -r 1 png/$outputFileName.mkv
