"""EN: Animated demonstration with extrema of strongly differing widths.
PL: Animacja algorytmu dla ekstremów o bardzo różnych szerokościach.
Run / Uruchom: python examples/demo_python.py
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from adaptive_sampling import sample_adaptive, interpolate_samples, interval_disagreement


def multiscale_function(x, _params=None):
    """EN: Broad maximum, medium dip, very narrow resonance + mild background.
    PL: Szerokie maksimum, średnio wąskie minimum i bardzo wąski rezonans.
    """
    return (0.04 + 0.035*x
            + 1.00*np.exp(-((x-2.45)/1.10)**2)
            - 0.72*np.exp(-((x-5.30)/0.36)**2)
            + 1.18*np.exp(-((x-7.43)/0.105)**2))


def main():
    root = Path(__file__).resolve().parents[1]
    assets = root / 'assets'
    assets.mkdir(exist_ok=True)
    # EN: Uniform initial grid: NO deliberately planted samples near x=7.43.
    # PL: Jednorodne punkty startowe: BEZ punktow dodanych przy x=7.43.
    initial_x = np.linspace(0, 10, 6)
    initial_y = np.array([multiscale_function(x) for x in initial_x])
    result = sample_adaptive(multiscale_function, 0, 10, None,
                       x0=initial_x, y0=initial_y,
                       max_new=90, max_points=180,
                       abs_tol=2e-4, rel_tol=3e-4, max_step=0.22,
                       method1='linear', method2='quadratic', jitter=0,
                       priority='disagreement', point_choice='difference',
                       exploration_fraction=0.25, seed=2026)
    # EN: Prepend the actual initial configuration before the first new evaluation.
    # PL: Dodaj klatkę początkową przed pierwszym nowym wywołaniem funkcji.
    frames = list(result.history)
    if frames[0][2] is not None:
        frames.insert(0, (initial_x.copy(), initial_y.copy(), None))
    xx = np.linspace(0, 10, 1800)
    truth = np.array([multiscale_function(v) for v in xx])
    fig, (ax, ae) = plt.subplots(2, 1, figsize=(10.2, 6.5), dpi=90,
                                 sharex=True, gridspec_kw={'height_ratios':[3,1]})

    def draw(k):
        ax.clear(); ae.clear()
        x, y, added = frames[k]
        q = np.unique(np.r_[xx, x])
        lin = interpolate_samples(x, y, q, 'linear')
        qua = interpolate_samples(x, y, q, 'quadratic')
        ax.plot(xx, truth, '--', color='#777777', lw=1.5, label='Reference / Funkcja referencyjna')
        ax.plot(q, lin, color='#0072b2', lw=1.35, label='Linear / Liniowa')
        ax.plot(q, qua, color='#d55e00', lw=1.65, label='Quadratic / Kwadratowa')
        ax.scatter(x, y, s=20, color='#141414', zorder=6, label='Samples / Próbki')
        if added is not None:
            ax.scatter([added], [multiscale_function(added)], s=110, facecolors='none',
                       edgecolors='#00885c', lw=2.4, zorder=7)
        ax.set(ylabel='f(x)', ylim=(-.95, 1.45), xlim=(0, 10))
        mode = 'Initial / Start' if k == 0 else result.info['selection_modes'][k-1]
        ax.set_title(f'75% adaptive / adaptacyjne + 25% exploration / eksploracja\n{len(x)} samples  |  +{max(0, len(x)-len(initial_x))} new  |  {mode}', fontsize=11)
        ax.grid(alpha=.18)
        ax.legend(loc='upper right', fontsize=8, ncol=2)
        # EN/PL: Show interpolation disagreement, not the true approximation error.
        ae.plot(q, abs(qua-lin), color='#743aa6', lw=1.2)
        ae.set(xlabel='x', ylabel='|Q - L|')
        ae.grid(alpha=.18)
        ae.set_ylim(bottom=0)
        fig.tight_layout()

    gif = assets/'adaptive_sampling.gif'
    animation = FuncAnimation(fig, draw, frames=len(frames), interval=180, repeat=True)
    animation.save(gif, writer=PillowWriter(fps=6))
    draw(len(frames)-1)
    fig.savefig(assets/'final_state.png', dpi=150)
    plt.close(fig)
    print(f'GIF: {gif} | frames={len(frames)}, initial={len(initial_x)}, new={result.info["new_evaluations"]}, explore={result.info["exploration_count"]}, nearest_to_peak={min(abs(result.x-7.43)):.5g}, stop={result.info["reason"]}')


if __name__ == '__main__':
    main()
