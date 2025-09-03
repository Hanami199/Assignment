import glob, numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio

nx = 256   # <-- set to your -x
ny = 256   # <-- set to your -y

files = sorted(glob.glob("frames/frame_*.f32"))
if not files:
    raise SystemExit("No frames found in frames/frame_*.f32")

# First pass: global min/max for consistent colors
vmin, vmax = np.inf, -np.inf
for f in files:
    a = np.fromfile(f, dtype=np.float32)
    a = a.reshape(ny, nx)
    vmin = min(vmin, float(np.nanmin(a)))
    vmax = max(vmax, float(np.nanmax(a)))

# Build frames
frames = []
cmap = plt.cm.inferno
for f in files:
    a = np.fromfile(f, dtype=np.float32).reshape(ny, nx)
    # normalize to [0,1] with global min/max
    if vmax > vmin:
        z = (a - vmin) / (vmax - vmin)
    else:
        z = np.zeros_like(a)
    img = (cmap(z)[:, :, :3] * 255).astype(np.uint8)  # RGB
    frames.append(img)

# GIF
imageio.mimsave("movie.gif", frames, duration=0.05, loop=0)  # 20 fps

# MP4 (needs imageio-ffmpeg installed)
try:
    imageio.mimsave("movie.mp4", frames, fps=20, quality=8)
except Exception as e:
    print("MP4 save skipped (install imageio-ffmpeg).", e)
