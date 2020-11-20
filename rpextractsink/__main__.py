#!/usr/bin/env python

from logging       import error as logging_error
from rpextractsink import genSink, build_args_parser
from brs_libs      import rpCache

def _cli():
    parser = build_args_parser()
    args  = parser.parse_args()

    if args.remove_dead_end==True or args.remove_dead_end=='True' or args.remove_dead_end=='true' or args.remove_dead_end=='t':
        remove_dead_end = True
    elif args.remove_dead_end==False or args.remove_dead_end=='False' or args.remove_dead_end=='false' or args.remove_dead_end=='f':
        remove_dead_end = False
    else:
        logging_error('Cannot interpret input -remove_dead_end: '+str(args.remove_dead_end))

    rpcache = rpCache('file', ['cid_strc'])
    genSink(rpcache,
            args.input_sbml,
            args.output_sink,
            remove_dead_end,
            args.compartment_id)

if __name__ == '__main__':
    _cli()
