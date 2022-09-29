import os, time, io
import pandas as pd
import subprocess

from pyw90.lib.config import Config
from pyw90.utility.utility import bc

import logging
logger = logging.getLogger(__name__)

class Job():
    r'''
    Monitoring job states.
    '''
    def __init__(self, config:Config):
        self.config = config
        self.path = config.path
        self.win = config.win
        self.local = config.local
        if self.local:
            self.localrun = config.localrun
        else:
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
        # jobs_str = os.popen("squeue -u '{0}'".format(self.usr_name)).read()
        # logger.debug(jobs_str)
        # df = pd.read_csv(io.StringIO(jobs_str), sep='\s+')

        # Test: check `cwd` usage
        # ls   = subprocess.run("ls -al",
        #                         cwd=self.path, 
        #                         capture_output=True, shell=True, text=True)
        # logger.debug(ls.stdout)
        # logger.debug(self.config.environ)
        fmt = r"%.7i %9P %10j %.8u %.2t %.12M %.5C %.4D"
        jobs = subprocess.run(["squeue", "-o", fmt, "-u", self.usr_name], 
                                cwd=self.path, 
                                env=self.config.environ,
                                capture_output=True, text=True)
        # logger.debug(jobs.stdout, jobs.stderr)
        df = pd.read_csv(io.StringIO(jobs.stdout), sep='\s+')

        # Not sure about job states codes is necessary or not
        # Ref:　https://slurm.schedmd.com/squeue.html

        # Print jobs nice
        logger.debug('SQUEUE RESULTS'.ljust(80, " "))
        for l in str(df).split('\n'):
            logger.debug(l.ljust(80, " "))
            
        df = df[df['NAME'] == self.job_name]
        # df = df[df['USER'] == self.usr_name]
        # logger.debug(df)
        return df

    def check_run_until_stop(self):
        r'''
        Check whether the `self.job_name` is running or not. 
        The program will abort the loop until the given job has done.
        With `logging.Logger` level at `logging.DEBUG`, check messages will be reported in log file 
        every `self.config.num_print_check` loops.
        '''
        check = 0
        logger.info(f'Job <{self.job_name}> is running!'.ljust(80, ' '))
        num_jobs = self.get_num_jobs()
        logger.debug(f'Job <{self.job_name}> has `num_jobs` {num_jobs}'.ljust(80, ' '))
        if self.local:
            logger.info(f'Job <{self.job_name}> run with PID {self.p.pid}'.ljust(80, ' '))
        while num_jobs > 0:
            check += 1
            if check % self.num_print_check == 0:
                logger.info(f'Job <{self.job_name}> is still running.'.ljust(80, ' '))
            time.sleep(self.check_time)
            num_jobs = self.get_num_jobs()

    def get_num_jobs(self):
        if self.local:
            if self.p.poll() is None:
                return 1
            else:
                return 0
        else:
            return len(self.get_jobs())

    def submit(self):
        r'''
        Submitting task of `self.config.run` file via the `sbatch` command.

        If `self.local` is True, run the task directly according self.localrun from `.yaml` file.

        Noticed: Allow to submit only when there are no job with the same name
        '''
        if self.local:
            p = subprocess.Popen(self.localrun.split(),
                                 cwd=self.path,
                                 env=self.config.environ,
                                 start_new_session=True)
            logger.debug(f'{self.job_name} run with PID: {p.pid}'.ljust(80, ' '))
            self.p = p
        else:
            num_jobs = len(self.get_jobs())
            if num_jobs == 0:
                # os.system('sbatch {0}'.format(self.run))
                subprocess.run(["sbatch", self.run],
                                 cwd=self.path,
                                 env=self.config.environ)
            else:
                raise ValueError('There is already {0} *{1}* task!'.format(num_jobs, self.job_name))

    def cancel(self):
        r"""
        Cancel all the job with name `self.job_name`
        """
        df = self.get_jobs()
        bc.cprint(bc.RED, f'Cancel all the jobs as listed')
        print(df)
        bc.cprint(bc.RED, 'Y(es) / N(o): ')
        if input().lower()[0] == 'y':
            for i in df['JOBID']:
                # os.popen(f"scancel {i}")
                subprocess.run(["scancel", str(i)],
                                 cwd=self.path,
                                 env=self.config.environ)
                logger.info(f"CANCEL JOB: {i}")
        else:
            return
