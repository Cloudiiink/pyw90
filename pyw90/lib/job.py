import os, time, io
import pandas as pd

from lib.config import Config

import logging
logger = logging.getLogger(__name__)

class Job():
    r'''
    Monitoring job states.
    '''
    def __init__(self, config:Config):
        self.config = config
        self.win = config.win
        self.run = config.run
        self.usr_name = config.usr_name
        self.job_name = config.job_name
        self.num_print_check = config.num_print_check
        self.check_time = config.check_time

    def get_jobs(self):
        r'''
        Get job list in `pandas.DataFrame` format via `squeue` command.
        '''
        # TODO: Add PBS support, or allow flexible input
        jobs_str = os.popen("squeue -u '{0}'".format(self.usr_name)).read()
        df = pd.read_csv(io.StringIO(jobs_str), sep='\s+')
        # Not sure about job states codes is necessary or not
        # Ref:　https://slurm.schedmd.com/squeue.html
        df = df[df['NAME'] == self.job_name]
        return df

    def check_run_until_stop(self):
        r'''
        Check whether the `self.job_name` is running or not. The program will abort the loop until the given job has done. With `logging.Logger` level at `logging.DEBUG`, check messages will be reported in log file every `self.config.num_print_check` loops.
        '''
        check = 0
        logger.info(f'Job <{self.job_name}> is running!'.ljust(80, ' '))
        num_jobs = len(self.get_jobs())
        while num_jobs > 0:
            check += 1
            if check % self.num_print_check == 0:
                logger.debug(f'Job <{self.job_name}> is still running.'.ljust(80, ' '))
            time.sleep(self.check_time)
            num_jobs = len(self.get_jobs())

    def submit(self):
        r'''
        Submitting task of `self.config.run` file via the `sbatch` command.

        Noticed: Allow to submit only when there are no job with the same name
        '''
        num_jobs = len(self.get_jobs())
        if num_jobs == 0:
            os.system('sbatch {0}'.format(self.run))
        else:
            raise ValueError('There is already {0} *{1}* task!'.format(num_jobs, self.job_name))
