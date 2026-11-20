"""EN: Regression tests for the translated sampling behavior.
PL: Testy regresyjne algorytmu próbkowania.
"""
import numpy as np
from adaptive_sampling import sample_adaptive, interpolate_samples, interval_disagreement


def test_interpolants_reproduce_nodes():
    x=np.array([0,.2,.6,1.1,2.]); y=np.array([0,.5,-.1,.8,.2])
    for method in ('linear','quadratic','pchip','spline'):
        np.testing.assert_allclose(interpolate_samples(x,y,x,method),y,atol=1e-12)


def test_quadratic_reproduces_polynomial():
    x=np.array([0,.2,.6,1.1,2.]); q=np.linspace(0,2,301)
    np.testing.assert_allclose(interpolate_samples(x,x*x-3*x+1,q,'quadratic'),q*q-3*q+1,atol=1e-12)


def test_linear_terminates_at_three_samples():
    r=sample_adaptive(lambda x,a: 2*x-0.6,-2,3,[],jitter=0,abs_tol=1e-9,rel_tol=0)
    assert len(r.x)==3 and r.info['reason']=='tolerance'


def test_no_reevaluation_on_resume():
    calls=[]
    def f(x,a):
        calls.append(x); return np.sin(3*x)
    first=sample_adaptive(f,0,2,[],max_new=15,max_points=25,jitter=0,abs_tol=1e-5,rel_tol=0)
    calls.clear()
    second=sample_adaptive(f,0,2,[],x0=first.x,y0=first.y,max_new=10,max_points=40,jitter=0,max_step=.10)
    assert len(calls)==second.info['new_evaluations']
    assert not set(calls).intersection(set(first.x))
    assert set(first.x).issubset(set(second.x))


def test_discrepancy_priority():
    x=np.array([0.,.5,1.,1.5,5.]); y=np.array([0.,1.,0.,0.,0.])
    diff,candidate=interval_disagreement(x,y,'linear','quadratic',n_tests=20)
    chosen=int(np.argmax(diff))
    assert chosen != int(np.argmax(np.diff(x)))
    r=sample_adaptive(lambda z,a: z*z,0,5,[],x0=x,y0=y,max_new=1,max_points=6,jitter=0,
                abs_tol=0,rel_tol=0,max_step=.1,method2='quadratic',
                exploration_fraction=0)
    assert x[chosen] < r.info['new_x'][0] < x[chosen+1]


def test_global_grid_and_vertex():
    r=sample_adaptive(lambda x,a: -(x-.35)**2,0,1,[],max_new=8,max_points=12,
                grid_mode='global',point_choice='vertex',method2='quadratic',jitter=0,
                abs_tol=0,rel_tol=0)
    assert len(r.x)>=3 and r.info['grid_mode']=='global'


def test_budget_and_validation():
    r=sample_adaptive(lambda x,a:x*x,0,1,[],max_new=3,max_points=3,jitter=0,abs_tol=0,rel_tol=0)
    assert r.info['new_evaluations']==3 and r.info['reason']=='point_limit'


def test_exploration_exact_longest_interval_midpoint():
    # EN: With 100% exploration every new x bisects the longest interval.
    # PL: Przy 100% eksploracji dzielimy najdluzszy przedzial na pol.
    x0=np.array([0., .5, 1., 4.])
    y0=np.array([0., 1., 0., 0.])
    r=sample_adaptive(lambda x,a:x*x,0,4,x0=x0,y0=y0,
                max_new=4,max_points=8,exploration_fraction=1.,
                jitter=.20,point_choice='vertex',seed=4,abs_tol=0,rel_tol=0)
    assert r.info['exploration_count']==4
    state=x0.copy()
    for nx in r.info['new_x']:
        i=int(np.argmax(np.diff(state)))
        assert abs(nx-0.5*(state[i]+state[i+1]))<1e-12
        state=np.sort(np.append(state,nx))


def test_mixed_sampling_seed_and_no_hidden_initial_resonance():
    # EN: The narrowly peaked function has uniform initial samples only.
    # PL: Funkcja z waskim pikiem ma jedynie rownomierne punkty poczatkowe.
    f=lambda x,a: .04+.035*x+np.exp(-((x-2.45)/1.1)**2)-.72*np.exp(-((x-5.3)/.36)**2)+1.18*np.exp(-((x-7.43)/.105)**2)
    x0=np.linspace(0,10,6)
    options=dict(x0=x0,y0=np.array([f(x,None) for x in x0]),
                 max_new=90,max_points=180,abs_tol=2e-4,rel_tol=3e-4,
                 max_step=.22,method1='linear',method2='quadratic',jitter=0,
                 exploration_fraction=.25,seed=2026)
    r=sample_adaptive(f,0,10,**options)
    r2=sample_adaptive(f,0,10,**options)
    np.testing.assert_array_equal(r.x,r2.x)
    assert r.info['exploration_count']>0
    assert min(abs(r.x-7.43))<0.03
    assert not np.any(abs(x0-7.43)<0.3)
    assert len(r.info['selection_modes'])==90
