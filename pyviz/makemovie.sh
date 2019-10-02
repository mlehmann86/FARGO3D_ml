###ffmpeg -i dgrz_%04d.png -vf scale="trunc(iw/2)*2:trunc(ih/2)*2, setpts=3*PTS" -c:v libx264 -profile:v high -pix_fmt yuv420p -r 25 bump.mp4
ffmpeg -i metal_%04d.png -vf scale="trunc(iw/2)*2:trunc(ih/2)*2, setpts=3*PTS" -c:v libx264 -profile:v high -pix_fmt yuv420p -r 25 metal.mp4
ffmpeg -i dgmid_%04d.png -vf scale="trunc(iw/2)*2:trunc(ih/2)*2, setpts=3*PTS" -c:v libx264 -profile:v high -pix_fmt yuv420p -r 25 dgmid.mp4
