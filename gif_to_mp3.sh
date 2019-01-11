for i in "$@"
do
        ffmpeg -y -r 24 -i $i -c:v mpeg4 -q:v 1 ${i%".gif"}.mp4
done

