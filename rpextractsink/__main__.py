#!/usr/bin/env python

from logging  import error as logging_error
from rpextractsink  import rpExtractSink, build_args_parser

def main(input_sbml, output_sink, remove_dead_end, compartment_id):
    rpcache = rpCache.rpCache()
    rpgensink = rpTool.rpExtractSink()
    rpgensink.cid_strc = rpcache.getCIDstrc()
    rpgensink.genSink(input_sbml, output_sink, remove_dead_end, compartment_id)


def _cli():
    parser = build_args_parser()
    args  = parser.parse_args()

    if args.remove_dead_end==True or args.remove_dead_end=='True' or args.remove_dead_end=='true' or args.remove_dead_end=='t':
        remove_dead_end = True
    elif args.remove_dead_end==False or args.remove_dead_end=='False' or args.remove_dead_end=='false' or args.remove_dead_end=='f':
        remove_dead_end = False
    else:
        logging.error('Cannot interpret input -remove_dead_end: '+str(params.remove_dead_end))

    rpgensink = rpExtractSink()
    rpgensink.genSink(args.input_sbml, args.output_sink, remove_dead_end, args.compartment_id)


if __name__ == '__main__':
    _cli()
