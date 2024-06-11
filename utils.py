import os  # noqa: F401
import sys # noqa: F401
from datetime import datetime


def time_stamp():
    return '[QTL-seq:{}]'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')) 
    # strftime will transform the datetime object to string format

def clean_cmd(cmd):
    return ' '.join(cmd.split())

def call_log(out_dir, name, cmd):
    print(time_stamp(), 
          '!!ERROR!! {}\n'.format(cmd), 
          flush=True) 
          # flush the data buffer right away, you need to monitor the log immediately


    print('please check below:\n')

    with open('{}/log/{}.log'.format(out_dir, name)) as log:
        for line in log:
            print(line, end='')