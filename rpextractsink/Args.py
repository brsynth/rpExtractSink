from argparse import ArgumentParser, ArgumentTypeError


def build_args_parser():
    parser = ArgumentParser(prog='rpextractsink', description='Generate the sink from a model SBML by specifying the compartment')
    parser = _add_arguments(parser)

    return parser


def _add_arguments(parser):
    parser.add_argument('input_sbml',
                        type=str,
                        help="input SBML file")
    parser.add_argument('output_sink',
                        type=str,
                        help="output sink file")
    parser.add_argument('--compartment_id',
                        type=str,
                        default='MNXC3',
                        help='SBML compartment id from which to extract the chemical species')
    parser.add_argument('--remove_dead_end',
                        type=str2bool, nargs='?',
                        const=True,
                        default='True',
                        help='upon FVA evaluation, ignore chemical species that do not have any flux')
    return parser


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ArgumentTypeError('Boolean value expected.')
