#!/usr/bin/env python

from rpextractsink import genSink, build_args_parser
from brs_libs      import rpCache


def _cli():
    parser = build_args_parser()
    args  = parser.parse_args()

    rpcache = rpCache('file', ['cid_strc'])
    genSink(rpcache,
            args.input_sbml,
            args.output_sink,
            args.remove_dead_end,
            args.compartment_id)


if __name__ == '__main__':
    _cli()
