import logging
from logging.config import dictConfig

logging_config = dict(
    version=1,
    formatters={
        'f': {'format': '%(asctime)s %(name)s %(levelname)-8s %(message)s',
              'datefmt': '%Y-%m-%d %H:%M:%S'}
    },
    handlers={
        'h': {'class': 'logging.StreamHandler',
              'formatter': 'f',
              'level': logging.DEBUG}
    },
    root={
        'handlers': ['h'],
        'level': logging.DEBUG,
    },
)


def create_logger(name='', logging_config=logging_config):
    """Create a Logger object

    Args:
        name (str, optional): The name of the module.function
        logging_config (TYPE, optional): Logging config

    Returns:
        TYPE: Description
    """
    dictConfig(logging_config)

    logger = logging.getLogger(name)

    return logger
