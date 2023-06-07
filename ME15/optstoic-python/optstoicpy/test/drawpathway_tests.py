import os
import unittest
from optstoicpy.script.utils import create_logger
from optstoicpy.core.pathway import Pathway
from optstoicpy.core.drawpathway import (
    draw_pathway)


class TestDrawPathway(unittest.TestCase):
    def setUp(self):
        self.logger = create_logger(name='Test core.drawpathway')
        self.pathway_fixture = {
            'flux': [-1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 1.0,
                     1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 2.0, 1.0,
                     1.0, 1.0, -1.0, 1.0],
            'iteration': 1,
            'reaction_id': ['R00200', 'R00300', 'R00658', 'R01059',
                            'R01063', 'R01512', 'R01518', 'R01519',
                            'R01538', 'R08570', 'EX_glc', 'EX_nad',
                            'EX_adp', 'EX_phosphate', 'EX_pyruvate',
                            'EX_nadh', 'EX_atp', 'EX_h2o', 'EX_nadp',
                            'EX_nadph']}

        self.p1 = Pathway(id=1,
                          name='OptStoic_pathway',
                          reaction_ids=self.pathway_fixture['reaction_id'],
                          fluxes=self.pathway_fixture['flux'])

    def test_draw_pathway(self):
        # Create png image
        draw_pathway(self.p1,
                     imageFileName='test_pathway',
                     imageFormat='png',
                     graphTitle=self.p1.name,
                     cleanup=True,
                     darkBackgroundMode=False)

        fname = 'test_pathway.png'
        self.assertEqual(os.path.exists(fname), True)

        # if os.path.exists(fname):
        #     os.remove(fname)
