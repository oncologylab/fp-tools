"""Independent numerical fixtures; never reads scientific inputs."""
import unittest
import numpy as np
from scipy.special import expit, roots_legendre
from scipy.stats import binom, betabinom, norm
from scipy.optimize import brentq
from shrinkage import posterior, fit_dispersion, threshold


def oracle(c, n, p, rho, prior_sd=1.):
    # Independent fixed Gauss-Legendre integration of SciPy distribution PMFs;
    # no production likelihood, mode search, or adaptive integrator is reused.
    nodes, weights = roots_legendre(768)
    def integrals(end, power=0):
        t = (nodes+1)*(end+20*prior_sd)/2-20*prior_sd
        mass = norm.pdf(t, scale=prior_sd)
        for ci, ni, pi in zip(c, n, p):
            q = expit(np.log(pi/(1-pi))-t) if 0 < pi < 1 else np.full_like(t, pi)
            if rho and 0 < pi < 1:
                mass *= betabinom.pmf(ci, ni, q*(1/rho-1), (1-q)*(1/rho-1))
            else:
                mass *= binom.pmf(ci, ni, q)
        return float(np.dot(weights, mass*t**power)*(end+20*prior_sd)/2)
    end = 20*prior_sd
    z = integrals(end)
    return [integrals(end,1)/z, 1-integrals(0)/z] + [
        brentq(lambda x: integrals(x)/z-q, -end, end, xtol=1e-10) for q in (.025,.5,.975)]


class ModelTests(unittest.TestCase):
    def test_independent_quadrature(self):
        cases = [([0,0],[0,0],[.2,.2]), ([0,0],[1,0],[.2,.3]),
                 ([0,0],[25,30],[.2,.3]), ([15,20],[25,30],[.2,.3]),
                 ([0,0],[10,30],[1e-8,2e-8]), ([0,0],[0,25],[.2,.3]),
                 ([2,0],[100,1],[.03,.8]), ([0,4],[10,10],[0.,.3]),
                 ([10,4],[10,10],[1.,.3]), ([5,10],[200,300],[.1,.2])]
        for rho in (0., .01, .2, .95):
            for c,n,p in cases:
                with self.subTest(rho=rho,c=c,n=n,p=p):
                    actual=posterior(*c,*n,*p,rho)
                    expected=oracle(c,n,p,rho)
                    np.testing.assert_allclose([actual['mean'],actual['probability_positive']],expected[:2],atol=1e-8,rtol=0)
                    np.testing.assert_allclose([actual['low'],actual['median'],actual['high']],expected[2:],atol=1e-6,rtol=0)

    def test_prior_and_support(self):
        result=posterior(0,0,0,0,.2,.3)
        self.assertEqual(result['mean'],0.)
        self.assertTrue(result['prior_only'])
        with self.assertRaises(ValueError): posterior(1,0,1,0,0.,.2)
        with self.assertRaises(ValueError): posterior(0,0,1,0,1.,.2)
        self.assertLess(posterior(0,0,1,0,.2,.3)['probability_positive'],.95)
        self.assertLess(posterior(9,8,10,10,.2,.2)['mean'],0)

    def test_threshold(self):
        self.assertTrue(np.isnan(threshold(np.arange(199))))
        self.assertEqual(threshold(np.arange(200)),190.)
        self.assertFalse(190. > threshold(np.arange(200)))

    def test_dispersion_no_identification(self):
        result=fit_dispersion(np.zeros((200,2)),np.ones((200,2)),np.full((200,2),.2))
        self.assertFalse(result['identified'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
