from unittest import TestCase
import GPFA as gp
import numpy as np

__author__ = 'ruben'
__doc__ = ''


class TestGpLv(TestCase):
    def test(self):
        Y, _, _ = gp.GpLv.toy_example()
        model = gp.GpLv(Y=Y, p=2)
        self.assertEquals(model.q, 100)
