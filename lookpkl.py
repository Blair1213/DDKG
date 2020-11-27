# -*- coding: utf-8 -*-
# @Time    : 2020-10-22 15:07
# @Author  : xiaorui su
# @Email   :  suxiaorui19@mails.ucas.edu.cn
# @File    : lookpkl.py
# @Software : PyCharm




import os
import gc
import time

import numpy as np
from collections import defaultdict
from keras import backend as K
from keras import optimizers

from utils import load_data, pickle_load, format_filename, write_log
from models import KGCN
from config import ModelConfig, PROCESSED_DATA_DIR,  ENTITY_VOCAB_TEMPLATE, \
    RELATION_VOCAB_TEMPLATE, ADJ_ENTITY_TEMPLATE, ADJ_RELATION_TEMPLATE, LOG_DIR, PERFORMANCE_LOG, \
    DRUG_VOCAB_TEMPLATE

dataset = "kegg"
print(pickle_load(format_filename(PROCESSED_DATA_DIR,DRUG_VOCAB_TEMPLATE,dataset=dataset)))