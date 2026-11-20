"""Faithful Python translation of the user-supplied MATLAB/Octave v7 algorithm.

EN: Interfaces and selection logic mirror probkuj_adaptacyjnie.m.
PL: Interfejs i logika wyboru odpowiadają probkuj_adaptacyjnie.m.

Author / Autor: Rafał Kotyński
Code written by / Kod napisał: ChatGPT (OpenAI)
"""
from dataclasses import dataclass
from typing import Any, Callable
import numpy as np
from scipy.interpolate import PchipInterpolator, CubicSpline


@dataclass
class V7Result:
    x: np.ndarray
    y: np.ndarray
    info: dict
    history: tuple


def interpolate_samples(x, y, q, method='quadratic'):
    """EN: Interpolate using v7's node-selection rule.

    PL: Interpoluj według zasady wyboru węzłów z oryginalnej wersji v7.
    """
    x, y = np.asarray(x, float).reshape(-1), np.asarray(y, float).reshape(-1)
    q0 = np.asarray(q, float)
    shape = q0.shape
    z = np.atleast_1d(q0).reshape(-1)
    if len(x) < 3 or x.shape != y.shape or np.any(np.diff(x) <= 0):
        raise ValueError('At least three sorted, distinct x nodes are required / Wymagane 3 węzły')
    if np.any(z < x[0]) or np.any(z > x[-1]):
        raise ValueError('Queries must be inside the sampled domain / Punkty poza zakresem')
    if method == 'linear':
        out = np.interp(z, x, y)
    elif method == 'pchip':
        out = PchipInterpolator(x, y)(z)
    elif method == 'spline':
        # MATLAB's interp1(x,y,q,'spline') uses not-a-knot end conditions.
        # PL: Warunki brzegowe not-a-knot odpowiadają interpolacji spline MATLAB-a.
        out = CubicSpline(x, y, bc_type='not-a-knot')(z)
    elif method == 'quadratic':
        i = np.clip(np.searchsorted(x, z, side='right') - 1, 0, len(x)-2)
        k0 = np.clip(i-1, 0, len(x)-3)
        ks = k0[:, None] + np.arange(3)[None, :]
        xp = x[ks]
        yp = y[ks]
        scale = xp[:, 2] - xp[:, 0]
        u = (xp - xp[:, :1]) / scale[:, None]
        v = (z - xp[:, 0]) / scale
        out = np.zeros_like(v)
        for j in range(3):
            p, r = (j+1) % 3, (j+2) % 3
            out += yp[:, j]*(v-u[:, p])*(v-u[:, r])/((u[:, j]-u[:, p])*(u[:, j]-u[:, r]))
    else:
        raise ValueError('method must be linear, quadratic, pchip, or spline')
    out = np.asarray(out).reshape(shape)
    return float(out) if shape == () else out


def interval_disagreement(x, y, method1='linear', method2='pchip',
                          grid_mode='local', n_tests=20, n_global=400):
    """EN: Maximum discrepancy and its candidate location in EACH interval.

    PL: Maksymalna różnica interpolantów i jej położenie w każdym przedziale.
    """
    x = np.asarray(x, float)
    h = np.diff(x)
    mid = x[:-1] + h/2
    if grid_mode == 'local':
        t = np.arange(1, n_tests+1)/(n_tests+1)
        grid = x[:-1, None] + h[:, None]*t[None, :]
        dif = abs(interpolate_samples(x, y, grid, method2)-interpolate_samples(x, y, grid, method1))
        ids = np.argmax(dif, axis=1)
        return dif[np.arange(len(h)), ids], grid[np.arange(len(h)), ids]
    if grid_mode == 'global':
        grid = x[0] + (np.arange(1, n_global+1)/(n_global+1))*(x[-1]-x[0])
        q = np.r_[grid, mid]
        d = abs(interpolate_samples(x, y, q, method2)-interpolate_samples(x, y, q, method1))
        best = d[n_global:].copy()
        candidate = mid.copy()
        for i in range(len(h)):
            indices = np.flatnonzero((grid > x[i]) & (grid < x[i+1]))
            if len(indices):
                index = indices[np.argmax(d[indices])]
                if d[index] > best[i]:
                    best[i], candidate[i] = d[index], grid[index]
        return best, candidate
    raise ValueError('grid_mode must be local or global')


def _vertex(x, y, i):
    """EN/PL: Local parabola vertex / Wierzchołek lokalnej paraboli."""
    k0 = max(0, min(i-1, len(x)-3))
    xx, yy = x[k0:k0+3], y[k0:k0+3]
    span = xx[2] - xx[0]
    if span <= 0:
        return None
    z = (xx - xx[0])/span
    d01 = (yy[1]-yy[0])/(z[1]-z[0])
    d12 = (yy[2]-yy[1])/(z[2]-z[1])
    A = (d12-d01)/(z[2]-z[0])
    B = d01-A*(z[0]+z[1])
    if not (np.isfinite(A) and np.isfinite(B)) or abs(A) <= 128*np.finfo(float).eps*max(1, *abs(yy)):
        return None
    v = xx[0] - B/(2*A)*span
    margin = max(32*np.finfo(float).eps*max(1, abs(x[i]), abs(x[i+1])), 1e-6*(x[i+1]-x[i]))
    return float(v) if x[i]+margin < v < x[i+1]-margin else None


def sample_v7(func: Callable[[float, Any], float], xmin: float, xmax: float,
              a=None, *, x0=None, y0=None, max_new=100, max_points=200,
              abs_tol=1e-3, rel_tol=1e-2, max_step=np.inf,
              method1='linear', method2='pchip', jitter=0.04,
              grid_mode='local', n_tests=20, n_global=400,
              priority='disagreement', point_choice='difference', seed=None, exploration_fraction=0.25):
    """Sample as in MATLAB v7, without graphics (use examples/demo_python.py).

    EN: Existing x0/y0 values are never reevaluated; max_new limits the number
    of *new* expensive function calls, including newly required initial nodes.
    PL: Wartości x0/y0 są wykorzystywane ponownie; max_new ogranicza liczbę
    *nowych* wywołań funkcji, wraz z brakującymi punktami startowymi.

    EN: Differences between interpolants are indicators, not certified bounds.
    PL: Różnice interpolantów nie są ścisłym oszacowaniem błędu.
    """
    if not np.isfinite([xmin, xmax]).all() or xmin >= xmax:
        raise ValueError('Finite bounds xmin < xmax required')
    if not (isinstance(max_new, int) and max_new >= 0 and isinstance(max_points, int) and max_points >= 3):
        raise ValueError('Invalid sample budgets')
    if not (np.isfinite(abs_tol) and abs_tol >= 0 and np.isfinite(rel_tol) and rel_tol >= 0):
        raise ValueError('Invalid tolerances')
    if not (max_step > 0 and not np.isnan(max_step)) or not 0 <= jitter <= .25:
        raise ValueError('Invalid maximum step or jitter')
    if method1 == method2 or any(m not in ('linear','quadratic','pchip','spline') for m in (method1, method2)):
        raise ValueError('Choose two different supported interpolation methods')
    if not np.isfinite(exploration_fraction) or not 0 <= exploration_fraction <= 1:
        raise ValueError('exploration_fraction must be in [0, 1]')
    if priority not in ('disagreement', 'step') or point_choice not in ('difference','vertex'):
        raise ValueError('Invalid priority / point choice')
    if grid_mode not in ('local','global') or n_tests < 1 or n_global < 1:
        raise ValueError('Invalid diagnostic grid')
    if (x0 is None) != (y0 is None):
        raise ValueError('x0 and y0 must both be provided')
    x = np.array([], dtype=float) if x0 is None else np.asarray(x0, dtype=float).reshape(-1)
    y = np.array([], dtype=float) if y0 is None else np.asarray(y0, dtype=float).reshape(-1)
    if len(x) != len(y) or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Invalid initial samples')
    order = np.argsort(x)
    x, y = x[order], y[order]
    if np.any(np.diff(x) <= 0) or np.any(x < xmin) or np.any(x > xmax):
        raise ValueError('Invalid initial x coordinates')
    rng = np.random.default_rng(seed)
    n_old = len(x)
    history = []
    new_x_history, delta_history = [], []
    vertices = 0
    selection_modes = []  # EN: One mode per newly selected interval; PL: tryb każdego wyboru.

    def evaluate(z):
        result = np.asarray(func(float(z), a))
        if result.ndim != 0 or not np.isrealobj(result) or not np.isfinite(result):
            raise ValueError('func(x,a) must return a finite real scalar')
        return float(result)

    def insert(z, deviation):
        nonlocal x, y
        value = evaluate(z)
        i = np.searchsorted(x, z)
        x = np.insert(x, i, z)
        y = np.insert(y, i, value)
        new_x_history.append(z)
        delta_history.append(deviation)
        history.append((x.copy(), y.copy(), z))

    needed = [xmin, (xmin+xmax)/2, xmax] if len(x) < 3 else [xmin, xmax]
    missing = [z for z in needed if z not in x]
    if len(x)+len(missing) > max_points or len(missing) > max_new:
        raise ValueError('Insufficient budget for initial samples')
    for z in missing:
        if z == (xmin+xmax)/2 and jitter:
            z += rng.uniform(-jitter,jitter)*(xmax-xmin)
            z = np.clip(z, xmin+.01*(xmax-xmin), xmax-.01*(xmax-xmin))
        insert(float(z), np.nan)
    if len(x) < 3:
        raise ValueError('At least three initial samples are needed')
    if not history:
        history.append((x.copy(), y.copy(), None))
    reason = 'point_limit'
    while True:
        differences, candidates = interval_disagreement(x, y, method1, method2, grid_mode, n_tests, n_global)
        widths = np.diff(x)
        limit = abs_tol + rel_tol*(y.max()-y.min())
        max_d = differences.max()
        wide = np.flatnonzero(widths > max_step)
        if not len(wide) and max_d <= limit:
            reason = 'tolerance'
            break
        if len(new_x_history) >= max_new or len(x) >= max_points:
            reason = 'point_limit'
            break
        numerical_step = 16*np.finfo(float).eps*max(1, abs(xmin), abs(xmax))
        eligible = np.flatnonzero(widths > 2*numerical_step)
        if not len(eligible):
            reason = 'numerical_resolution'
            break
        # EN: Independent Bernoulli decision at each iteration: 25% exploration
        # by bisecting the longest eligible interval, 75% existing criterion.
        # PL: Losujemy niezależnie w każdej iteracji: z prawdopodobieństwem
        # 25% eksplorujemy środek najdłuższego przedziału, a w 75% adaptujemy.
        explore = rng.random() < exploration_fraction
        by_step = (priority == 'step' and len(wide)) or (max_d <= limit and len(wide))
        if explore:
            i = eligible[np.argmax(widths[eligible])]
            z = float(x[i] + widths[i]/2)
            selection_mode = 'exploration'
        elif by_step:
            possible = np.intersect1d(eligible, wide)
            if not len(possible):
                reason = 'numerical_resolution'
                break
            i = possible[np.argmax(widths[possible])]
            z = float(x[i] + widths[i]/2)
            selection_mode = 'step_fallback'
        else:
            i = eligible[np.argmax(differences[eligible])]
            scale = 64*np.finfo(float).eps*max(1, *abs(y))
            z = float(x[i] + widths[i]/2 if differences[i] <= scale else candidates[i])
            if point_choice == 'vertex':
                vertex = _vertex(x, y, int(i))
                if vertex is not None:
                    z = vertex
                    vertices += 1
            if jitter:
                z += rng.uniform(-jitter, jitter)*widths[i]
                margin = max(16*np.finfo(float).eps*max(1, abs(x[i]), abs(x[i+1])), .01*widths[i])
                z = float(np.clip(z, x[i]+margin, x[i+1]-margin)) if margin*2 < widths[i] else float(x[i]+widths[i]/2)
            selection_mode = 'disagreement'
        # EN: Exploration/fallback points are EXACT midpoints: no vertex/jitter.
        # PL: W eksploracji (i kroku awaryjnym) wybieramy dokładny środek.
        if not x[i] < z < x[i+1]:
            reason = 'numerical_resolution'
            break
        insert(z, float(max_d))
        selection_modes.append(selection_mode)
    final_d, _ = interval_disagreement(x, y, method1, method2, grid_mode, n_tests, n_global)
    spacing = np.diff(x)
    info = dict(method1=method1, method2=method2, grid_mode=grid_mode,
                n_tests=n_tests, n_global=n_global, priority=priority,
                point_choice=point_choice, jitter=jitter,
                exploration_fraction=float(exploration_fraction),
                selection_modes=tuple(selection_modes),
                exploration_count=selection_modes.count('exploration'),
                disagreement_count=selection_modes.count('disagreement'),
                step_fallback_count=selection_modes.count('step_fallback'),
                new_evaluations=len(new_x_history), previous_samples=n_old,
                total_points=len(x), new_x=np.asarray(new_x_history),
                max_discrepancy_history=np.asarray(delta_history),
                chosen_vertices=vertices, max_discrepancy=float(final_d.max()),
                max_spacing=float(spacing.max()), min_spacing=float(spacing.min()),
                mean_spacing=float((x[-1]-x[0])/(len(x)-1)),
                mean_to_min_spacing=float((x[-1]-x[0])/(len(x)-1)/spacing.min()),
                tolerance=float(abs_tol+rel_tol*(y.max()-y.min())), reason=reason)
    return V7Result(x=x, y=y, info=info, history=tuple(history))
