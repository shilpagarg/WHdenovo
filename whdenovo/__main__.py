import sys
import importlib
import logging
import argparse
from .args import HelpfulArgumentParser
#from . import __version__
#from .args import HelpfulArgumentParser

__version__ = '1.0'
# List of all subcommands. A module of the given name must exist and define
# add_arguments() and main() functions. Documentation is taken from the first
# line of the module’s docstring.
COMMANDS = ['simulate', 
            'partition', 
            'validate', 
            'assemble'
            ]

logger = logging.getLogger(__name__)


class NiceFormatter(logging.Formatter):
	"""
	Do not prefix "INFO:" to info-level log messages (but do it for all other
	levels).

	Based on http://stackoverflow.com/a/9218261/715090 .
	"""
	def format(self, record):
		if record.levelno != logging.INFO:
			record.msg = '{}: {}'.format(record.levelname, record.msg)
		return super().format(record)


def setup_logging(debug):
	"""
	Set up logging. If debug is True, then DEBUG level messages are printed.
	"""
	handler = logging.StreamHandler()
	handler.setFormatter(NiceFormatter())
	root = logging.getLogger()
	root.addHandler(handler)
	root.setLevel(logging.DEBUG if debug else logging.INFO)


def ensure_pysam_version():
	from pysam import __version__ as pysam_version
	from distutils.version import LooseVersion
	if LooseVersion(pysam_version) < LooseVersion("0.8.1"):
		sys.exit("WhatsHap requires pysam >= 0.8.1")

def main(argv = sys.argv[1:]):
    all_args = sys.argv
    print('Running commands: whdenovo', ' '.join(all_args[1:]))
    #parser = argparse.ArgumentParser()
    parser = HelpfulArgumentParser(description=__doc__, prog='whdenovo')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    subparsers = parser.add_subparsers()
    for command_name in COMMANDS:
        module = importlib.import_module('.'+command_name, 'whdenovo')
        subparser = subparsers.add_parser(command_name, help = module.__doc__, description = module.__doc__)
        subparser.set_defaults(module=module, subparser=subparser)
        module.add_arguments(subparser)
    args = parser.parse_args(argv)
    if not hasattr(args, 'module'):
        parser.error('Please provide the name of a subcommand to run')
    else:
        module = args.module
        if hasattr(args.module, 'validate'):
            subparser = args.subparser
            args.module.validate(args, subparser)
        del args.subparser
        del args.module
        if all_args[1] == 'simulate':
            if len(all_args) >= 3:
                module.main(all_args[2], args)
            else:
                module.main(None, args)
        else:
            module.main(args)

if __name__ == '__main__':
    __name__ = 'whdenovo'
    main()
