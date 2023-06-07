from __future__ import print_function
import unittest
import os
from optstoicpy.script.utils import create_logger
from optstoicpy.core.pathway import (
    Pathway,
    generate_kegg_model
)


class TestPathway(unittest.TestCase):
    def setUp(self):
        self.logger = create_logger(name='Test core.Pathway')
        self.pathway_fixture = {'flux': [-1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 1.0,
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
                          name='OptStoic',
                          reaction_ids=self.pathway_fixture['reaction_id'],
                          fluxes=self.pathway_fixture['flux'])

    @unittest.skip("Need to update test!")
    def test_rearrange_pathway(self):
        self.logger.info("Test rearranging reaction order")
        self.p1.rearrange_reaction_order()

    def test_kegg_model_generation(self):
        self.logger.info(
            "Creating 'res' folder in the current directory if not exist...")
        # outputFilepath = 'res'
        # outputFilename = 'OptStoic'
        # try:
        #     os.makedirs(outputFilepath)
        # except OSError:
        #     if not os.path.isdir(outputFilepath):
        #         raise Exception

        self.logger.info("Test create KEGG model file")

        filename = "./test_kegg_model_generation.txt"
        f = open(filename, 'a+')
        kegg_model_text = generate_kegg_model(self.p1, filehandle=f)
        print(kegg_model_text)
        self.assertIn('R01512', kegg_model_text)
        self.assertIn('R01512', kegg_model_text)
        f.close()

        if os.path.exists(filename):
            os.remove(filename)
