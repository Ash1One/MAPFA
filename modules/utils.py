#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import errno
import logging
import logging.config
import os
import yaml
from pathlib import Path
from pathlib import PurePosixPath

def fileExists(filename):
    try:
        return filename and os.path.exists(filename) and os.path.getsize(filename) > 0
    except OSError:
        return False

def makesurePathExists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
        else:
            logger = logging.getLogger('mapfa')
            logger.warning("%s has existed.")
            if not os.listdir(path):
                logger.warning("%s is not empty.")

def getLogger(outdir, silent):
    # get the log_config file path
    config_path = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'config', 'log_config.yaml')
    # log config
    with open(config_path, 'r', encoding='utf-8') as f:
        log_config = yaml.safe_load(f)
        temp = log_config['handlers']['file']['filename']
        log_config['handlers']['file']['filename'] = str(Path.resolve(Path(outdir)).joinpath(temp))
        logging.config.dictConfig(log_config)
    if silent:
        logger = logging.getLogger('silentlogger')
    else:
        # here it should assign 'root' to name of getLogger, othrewise logger would failed to log.
        logger = logging.getLogger('mapfa')
    return logger