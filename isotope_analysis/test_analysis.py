"""Numerical regression and limiting-case checks for analytical conversions."""
from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd
from isotope_analysis import gcms_analysis as gc
from isotope_analysis.source_fractions import gas_fraction,analyze

HERE=Path(__file__).resolve().parent


class AnalyticalTests(unittest.TestCase):
    def test_response_family(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/"gcms"
            result=gc.analyze(HERE/"inputs/ethylene_mass_spectra.csv",out).set_index("condition")
            expected={"13gas_12liquid":(644,.8675841044018286,.9149975081893805),
                      "12gas_13liquid":(374,.2039195560944763,.21545741576245442)}
            for feed,(count,low,high) in expected.items():
                row=result.loc[feed]
                self.assertEqual(row.accepted_count,count)
                np.testing.assert_allclose([row.f13_atom_min,row.f13_atom_max],[low,high],rtol=2e-10)
                self.assertGreater(row.mixed_depletion_min,0)
            fits=pd.read_csv(out/"accepted_response_fits.csv")
            self.assertEqual(len(fits),1018)
            np.testing.assert_allclose(fits.p_12_12+fits.p_12_13+fits.p_13_13,1,atol=2e-15)
            np.testing.assert_allclose(fits.random_mixed,2*fits.f13_atom*(1-fits.f13_atom),atol=2e-15)

    def test_endmembers(self):
        np.testing.assert_allclose(gas_fraction([.0107,.99],"13gas_12liquid"),[0,1])
        np.testing.assert_allclose(gas_fraction([.0107,.99],"12gas_13liquid"),[1,0])
        self.assertGreater(gas_fraction(1.,"13gas_12liquid"),1)
        for feed in ["13gas_12liquid","12gas_13liquid"]:
            for g in [0.,.1,.5,.9,1.]:
                f=.0107+.9793*g if feed=="13gas_12liquid" else .99-.9793*g
                self.assertAlmostEqual(float(gas_fraction(f,feed)),g)
        with self.assertRaises(ValueError):gas_fraction(.5,"13gas_12liquid",.99,.99)
        with self.assertRaises(ValueError):gas_fraction(-.1,"13gas_12liquid")

    def test_source_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            per,summary=analyze(HERE/"inputs/nmr_reporter_intervals.csv",
                HERE/"reference_outputs/gcms/ethylene_response_envelopes.csv",Path(tmp)/"sources")
            self.assertEqual(len(per),10)
            self.assertEqual(len(summary),5)
            eth=summary.set_index("product").loc["ethylene"]
            np.testing.assert_allclose([eth.gas_source_low,eth.gas_source_high],
                                      [.7909145146916548,.9234121394765523],rtol=1e-10)
            self.assertTrue(per.loc[per["product"].eq("ethanol") & per.feed.eq("13gas_12liquid"),
                                    "physical_boundary_clipped"].iloc[0])
            np.testing.assert_allclose(per.liquid_source_low,1-per.gas_source_high)
            np.testing.assert_allclose(per.liquid_source_high,1-per.gas_source_low)


if __name__=="__main__":unittest.main()
