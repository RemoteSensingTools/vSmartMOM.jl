"""Small analytic tests for Round7 interpolation, basis, instrument and sign."""

import unittest
import numpy as np

from build_round7_observations import (
    albedo_spectrum, bracket, finite, instrument_operator, spectral_weights,
)


class Round7ObservationTests(unittest.TestCase):
    def test_exact_and_interior_nodes(self):
        self.assertEqual(bracket([500, 750, 1000], 1000), [(2, 1.)])
        self.assertEqual(bracket([500, 750, 1000], 625), [(0, .5), (1, .5)])

    def test_forbid_extrapolation(self):
        for value in (499, 1001, float("nan")):
            with self.assertRaises(ValueError):
                bracket([500, 750, 1000], value)
        with self.assertRaises(ValueError):
            spectral_weights(np.array([0., 1.]), np.array([1.001]))

    def test_spectral_weights_preserve_linear_function(self):
        grid = np.array([10., 20., 40.])
        target = np.array([10., 17., 20., 29., 40.])
        lo, hi, w = spectral_weights(grid, target)
        np.testing.assert_allclose((1-w)*(3*grid[lo]+2) + w*(3*grid[hi]+2), 3*target+2)

    def test_canonical_legendre_coordinate(self):
        nu = np.array([10., 15., 20.])
        c = np.array([.3, .02, -.01])
        np.testing.assert_allclose(albedo_spectrum(c, nu, (10, 20)), [.27, .305, .31])
        wider = np.array([5., 10., 15., 20., 25.])
        np.testing.assert_array_equal(albedo_spectrum(c, wider, (10, 20))[1:-1],
                                      albedo_spectrum(c, nu, (10, 20)))

    def test_instrument_density_and_signed_analyzer(self):
        wavelength = np.linspace(757., 773., 3201)
        nu = 1e7/wavelength
        stokes = np.array([2., 3., 4.])[:, None] * wavelength[None, :]**2/1e7
        out = instrument_operator(nu, stokes, np.array([.5, -.2, .1, 0.]), np.array([758., 765., 772.]))
        np.testing.assert_allclose(out, 2*.5 + 3*.2 + 4*.1, atol=1e-14)

    def test_instrument_linearity_and_subtraction_closure(self):
        nu = np.linspace(1e7/773, 1e7/757, 2735)
        x = np.linspace(0, 10, len(nu))
        ray = np.array([2 + np.sin(x), .3*np.cos(x), -.1*np.sin(x)])
        cab = .97*ray
        rrs = .02*(ray + 1)
        targets = np.linspace(758, 772, 40)
        analyzer = np.array([.5, -.47, -.17, 0.])
        H = lambda st: instrument_operator(nu, st, analyzer, targets)
        delta = H(cab+rrs-ray)
        np.testing.assert_allclose(H(cab+rrs)-delta, H(ray), atol=1e-13)
        self.assertGreater(np.max(np.abs(H(cab+rrs)+delta-H(ray))), 1e-2)

    def test_ils_shoulder_guard(self):
        nu = np.linspace(1e7/772, 1e7/758, 100)
        with self.assertRaises(ValueError):
            instrument_operator(nu, np.ones((3, 100)), [.5, 0, 0], np.array([758, 772]))

    def test_mask_and_fill_rejection(self):
        for value in (np.array([1e20]), np.array([np.nan]), np.ma.array([1.], mask=[True])):
            with self.assertRaises(ValueError):
                finite(value, "test")

    def test_identical_additive_noise(self):
        rng = np.random.RandomState(3)
        y = rng.uniform(10, 20, 100)
        correction = rng.uniform(-1, 1, 100)
        noise = rng.uniform(-.1, .1, (11, 100))
        noise[-1] = 0
        new = (y-correction)[None, :] + noise
        np.testing.assert_allclose(new, (y[None, :]+noise)-correction, rtol=0, atol=5e-15)
        np.testing.assert_array_equal(new[-1], y-correction)


if __name__ == "__main__":
    unittest.main()
