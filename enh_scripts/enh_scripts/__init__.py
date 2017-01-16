import logging
import sys

from enh_scripts.utils.utils import open_or_stdin

HUMAN_CHROMS = ["chr%s" % i for i in range(1, 23)] + ["chrX", "chrY"]

LOG_DEFAULT_FORMAT = "%(asctime)s [%(levelname)s] - %(message)s"

LOG_LEVEL_MAP = {
    "debug": logging.DEBUG,
    "info": logging.INFO,
    "warn": logging.WARN,
}


def add_logging_args(parser):
    parser.add_argument("--log_level", choices=LOG_LEVEL_MAP.keys(), default='info')
    parser.add_argument("--log_file")


def init_logger(args):
    logger = logging.getLogger()
    logger.setLevel(LOG_LEVEL_MAP[args.log_level])
    formatter = logging.Formatter(LOG_DEFAULT_FORMAT)
    if not args.log_file:
        handler = logging.StreamHandler(stream=sys.stdout)
    else:
        handler = logging.FileHandler(args.log_file)
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger
