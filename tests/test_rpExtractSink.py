"""
Created on June 17 2020

@author: Joan Hérisson
"""

# Generic for test process
from unittest import TestCase

# Specific for tool
from rpextractsink import genSink
from brs_libs      import rpCache

# Specific for tests themselves
from hashlib  import sha256
from pathlib  import Path
from tempfile import NamedTemporaryFile
from filecmp  import cmp
from os       import path as os_path



# Cette classe est un groupe de tests. Son nom DOIT commencer
# par 'Test' et la classe DOIT hériter de unittest.TestCase.
# 'Test_' prefix is mandatory
class Test_rpExtractSink(TestCase):

    rpcache = rpCache('file', ['cid_strc'])

    def test_genSink(self):
        outfile = NamedTemporaryFile(delete=True).name
        genSink(self.rpcache,
                input_sbml=os_path.join('data', 'e_coli_model.sbml'),
                output_sink=outfile)
        self.assertTrue(
            cmp(Path(outfile), os_path.join('data', 'output_sink.csv')))
        outfile.close()
