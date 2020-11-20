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



# Cette classe est un groupe de tests. Son nom DOIT commencer
# par 'Test' et la classe DOIT hériter de unittest.TestCase.
# 'Test_' prefix is mandatory
class Test_rpExtractSink(TestCase):

    rpcache = rpCache('file', ['cid_strc'])

    def test_genSink(self):
        outfile = NamedTemporaryFile(delete=True)
        genSink(self.rpcache,
                input_sbml=os_path.join('data', 'e_coli_model.sbml'),
                output_sink=outfile)
        self.assertEqual(
            sha256(Path(outfile.name).read_bytes()).hexdigest(),
            '37be178695f79f7fe58d50a799dbc6f1840793157ad40fd1596d54311222f0760d6ebb46e3ba4d72782b3ab8ac8228e9045868b3f2c519c9e3dfeb32375b7df8'
                        )
        outfile.close()
