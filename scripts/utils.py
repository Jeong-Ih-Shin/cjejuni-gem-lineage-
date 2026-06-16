"""
Shared utilities: paths, color palette, plot style.

All scripts import from this module. Paths are resolved relative to the
repository root (one level above scripts/).
"""
import matplotlib
matplotlib.use('Agg')   # non-interactive backend — works in headless / WSL / SSH

from pathlib import Path
import matplotlib.pyplot as plt

# === Paths (relative to repository root) ===
REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / 'data'
GEMS_DIR = REPO_ROOT / 'gems'
INTERMEDIATE_DIR = REPO_ROOT / 'intermediate'
OUTPUT_DIR = REPO_ROOT / 'output'
FIG_DIR = OUTPUT_DIR / 'figures'

for d in (INTERMEDIATE_DIR, FIG_DIR):
    d.mkdir(parents=True, exist_ok=True)

# === Substrates measured in vitro ===
CHEMICALS = ['Lactate', 'Serine', 'Proline', 'Asparagine', 'Formate', 'Sulfite', 'Gluconate']

# === Editorial-style color palette ===
NPG = {
    'red':        '#E64B35',
    'blue':       '#4DBBD5',
    'green':      '#00A087',
    'darkblue':   '#3C5488',
    'salmon':     '#F39B7F',
    'graypurple': '#8491B4',
    'mint':       '#91D1C2',
    'brown':      '#7E6148',
}

SEM = {
    'text':    '#1a1a1a',
    'line':    '#333333',
    'grid':    '#E0E0E0',
    'neutral': '#BBBBBB',
    'lost':    NPG['salmon'],
    'gained':  NPG['green'],
    'up':      NPG['green'],
    'down':    NPG['salmon'],
}

CC_COLORS = [
    NPG['red'], NPG['blue'], NPG['green'], NPG['darkblue'],
    NPG['salmon'], NPG['graypurple'], NPG['mint'], NPG['brown'],
]


def apply_style():
    plt.rcParams.update({
        'font.family':        'sans-serif',
        'font.sans-serif':    ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
        'mathtext.fontset':   'stixsans',  # sans-serif math, matches Arial/Liberation Sans
        'mathtext.default':   'regular',  # use body font for math
        'font.size':          8,
        'axes.titlesize':     9,
        'axes.labelsize':     8,
        'xtick.labelsize':    7,
        'ytick.labelsize':    7,
        'legend.fontsize':    7,
        'xtick.labelsize':    7,
        'ytick.labelsize':    7,
        'legend.fontsize':    7,
        'axes.spines.top':    False,
        'axes.spines.right':  False,
        'axes.edgecolor':     SEM['line'],
        'axes.linewidth':     0.8,
        'axes.titleweight':   'bold',
        'axes.titlelocation': 'left',
        'axes.titlepad':      12,
        'text.color':         SEM['text'],
        'figure.facecolor':   'white',
        'savefig.facecolor':  'white',
        'savefig.dpi':        300,
        'axes.axisbelow':     True,
        'svg.fonttype':       'none',  # keep text editable in Inkscape
    })


def save_fig(fig, basename):
    """Save figure as both PNG (300 dpi) and SVG."""
    fig.savefig(FIG_DIR / f'{basename}.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIG_DIR / f'{basename}.svg', bbox_inches='tight')
    print(f'  Saved: {basename}.png + {basename}.svg')


def bounding_ellipse(points, ax, color, alpha_fill=0.18, n_std=2.0):
    """Draw covariance-based ellipse rotated to match point cluster orientation."""
    from matplotlib.patches import Ellipse
    import numpy as np
    if len(points) < 2:
        return
    cov = np.cov(points.T)
    if cov.ndim == 0:
        return
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    width, height = 2 * n_std * np.sqrt(np.maximum(vals, 0.01))
    cx, cy = points.mean(axis=0)
    ax.add_patch(Ellipse((cx, cy), width, height, angle=angle,
                         facecolor=color, edgecolor='none', alpha=alpha_fill, zorder=1))
    ax.add_patch(Ellipse((cx, cy), width, height, angle=angle,
                         facecolor='none', edgecolor=color, alpha=0.75, linewidth=1.5, zorder=2))
