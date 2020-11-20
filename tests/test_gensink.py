"""
Created on Jul 15 2020

@author: Joan Hérisson
"""

from _main         import Main
from os            import path as os_path
from rpextractsink import rpExtractSink

class Test(Main):
    __test__ = True

    mod_name  = 'rpextractsink'
    # cls_name  = 'rpExtractSink'
    obj       = rpExtractSink()
    func_name = 'genSink'
    cmd  = (os_path.join('data', 'e_coli_model.sbml') \
           + ' ' \
           + os_path.join('data', 'test_rpExtractSink.csv')).split()
    bap  = getattr(__import__(mod_name), 'build_args_parser')
    args = bap().parse_args(cmd)

    files = [(args.output_sink, '4c7c936f4bf862a8f22f7c2f7b7f0b063779dda93657be943cbb6cadd126ababaf4af66f0c412bce48a69f73acc2a926c14c416d1c5712da31e969adfd1aa305')]
