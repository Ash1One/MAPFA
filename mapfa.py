#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import logging
import logging.config
import errno
import os
import yaml

import modules


def confirmPathExists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def get_logger(output_dir, silent):
    confirmPathExists(output_dir)
    if silent:
        logger = logging.getLogger('silentlogger')
    return logger

def main(args=None):
    parser = argparse.ArgumentParser(description="MAPFA: Metagenomic Analysis Pipeline for Food Animals", prog="mapfa")

with open('log_config.yaml', 'r', encoding='utf-8') as f:
    log_config = yaml.load(f)
    logging.config.dictConfig(log_config)



##############################
# Pre-QC       
##############################

##############################
# readFiltering       
##############################

##############################
# QC-again       
##############################

##############################
# Assembly       
##############################

##############################
# Assembly QC      
##############################

##############################
# Binning      
##############################

##############################
# Binning QC     
##############################

mapfa_test = 'ok'
assert mapfa_test == 'ok', 'something was wrong when running mapfa!'