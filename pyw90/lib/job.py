import time, io
import subprocess
from tabulate import tabulate

from pyw90.utility.utility import bc

import logging
logger = logging.getLogger(__name__)

class Job(dict):
    '''
    Monitoring job states.

    Attributes
    ----------
    path : str
        The path of the working folder.
    win : str
        The name of the win file.
    local : bool
        Whether the job is running locally.
    localrun : str
        The command to run the job locally. Example: `wannier90.x seedname`.
    run : str
        The name of the script submitting the SLURM job. 
    usr_name : str
        Username of current user.
    job_name : str
        Job name in the SLURM queue.
    num_print_check : int
        The number of loops to print the running message.
    check_time : float
        The time interval (seconds) to check the job status.
    p : subprocess.Popen
        The process of the job when running locally.
    jobs : list
        The list of jobs in the SLURM queue in dictionary format.
    '''
    def __init__(self, kwargs={},
                 path:str=None, win:str=None, local:bool=None, localrun=False, run:str=None, environ=None,
                 usr_name:str=None, job_name:str=None, num_print_check:int=1, check_time=60.):
        
        super().__init__(**kwargs)
        self.__dict__ = self
        
        self.path = kwargs.get('path', path)
        self.win = kwargs.get('win', win)
        self.local = kwargs.get('local', local)
        self.localrun = kwargs.get('localrun', localrun)
        self.run = kwargs.get('run', run)
        self.environ = kwargs.get('environ', environ)
        self.usr_name = kwargs.get('usr_name', usr_name)
        self.job_name = kwargs.get('job_name', job_name)
        self.num_print_check = kwargs.get('num_print_check', num_print_check)
        self.check_time = kwargs.get('check_time', check_time)
        self.p = None
        self.jobs            = []


    def get_jobs(self, default_usr=True, default_job=True):
        '''
        Get job list in `pandas.DataFrame` format via `squeue` command.

        Parameters
        ----------
        default_usr : bool
            Whether to use the username in ``self.usr_name``.
        default_job : bool
            Whether to use the job name in ``self.job_name``.
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
        # logger.debug(self.environ)
        
        # Use `-o` option to get the job information with specific format
        # The meaning of each column is as follows:
        #   JOBID PARTITION NAME USER ST TIME CPUS NODES
        # You can also use `squeue -O` to use a self-explanatory format
        # `fmt = r"JobID,Partition,Name,Username,StateCompact,Timeused,NumCPUs,NumNodes"`
        
        # Ref:　https://slurm.schedmd.com/squeue.html

        fmt = r"%i %P %j %u %t %M %C %D"
        # fmt = r"%.7i %9P %10j %.8u %.2t %.12M %.5C %.4D"
        keys = ['JOBID', 'PARTITION', 'NAME', 'USER', 'ST', 'TIME', 'CPUS', 'NODES']

        cmd = ["squeue", "-o", fmt]
        if default_usr:
            cmd += ["-u", self.usr_name]
        try:
            jobs = subprocess.run(cmd, cwd=self.path, 
                                  env=self.environ,
                                  capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"Error: {e.stderr}")

        self.jobs = []
        for line in jobs.stdout.strip().split('\n')[1:]:
            job = {k : v.strip() for k, v in zip(keys, line.split())}
            if default_job and job['NAME'] != self.job_name:
                continue
            self.jobs.append(job)

        logger.debug(self.jobs)
        return self.jobs
            
        # # logger.debug(jobs.stdout, jobs.stderr)
        # df = pd.read_csv(io.StringIO(jobs.stdout), sep=r'\s+')

        # # Printing jobs nicely
        # logger.debug('FULL SQUEUE RESULTS'.ljust(80, " "))
        # for l in str(df).split('\n'):
        #     logger.debug(l.ljust(80, " "))
            
        # df = df[df['NAME'] == self.job_name]
        # return df
    
    def display_jobs(self, jobs):
        '''
        Display jobs in a table format.
        '''
        if len(jobs) == 0:
            logger.info(f'No job is found.')
        
        table = tabulate(jobs, headers="keys", tablefmt='pretty')
        logger.info(table)

    def check_run_until_stop(self):
        '''
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
                logger.debug(f'Job <{self.job_name}> has `num_jobs` {num_jobs}'.ljust(80, ' '))
            time.sleep(self.check_time)
            num_jobs = self.get_num_jobs()

    def get_num_jobs(self):
        """
        Returns the number of jobs from ``self.usr_name`` with name ``self.job_name``.
        
        If the instance is running locally, it checks if the process is still running.
        If the process is running, it returns 1, otherwise it returns 0.
        
        If the instance is running via SLURM, it returns the number of jobs retrieved
        using the `get_jobs` method with `default_job=True` and `default_usr=True`.
        
        Returns
        -------
        num_jobs : int
            The number of jobs associated with the current instance.
        """
        if self.local:
            if self.p.poll() is None:
                return 1
            else:
                return 0
        else:
            return len(self.get_jobs(default_job=True, default_usr=True))

    def submit(self):
        r'''
        Submitting task of `self.config.run` file via the `sbatch` command.

        If `self.local` is True, run the task directly according self.localrun from `.yaml` file.

        Noticed: Allow to submit only when there are no job with the same name
        '''
        if self.local:
            p = subprocess.Popen(self.localrun.split(),
                                 cwd=self.path,
                                 env=self.environ,
                                 start_new_session=True)
            logger.debug(f'{self.job_name} run with PID: {p.pid}'.ljust(80, ' '))
            self.p = p
        else:
            num_jobs = len(self.get_jobs(default_job=True, default_usr=True))
            if num_jobs == 0:
                subprocess.run(["sbatch", self.run],
                                 cwd=self.path,
                                 env=self.environ)
            else:
                raise ValueError('There is already {0} *{1}* task!'.format(num_jobs, self.job_name))

    def cancel(self, job_name:str=None):
        """
        Cancel all the job with name ``job_name``. Default is ``self.job_name``.
        """
        jobs = self.get_jobs(default_job=True, default_usr=True)
        bc.cprint(bc.RED, f'Cancel all the jobs as listed')
        self.display_jobs(jobs)
        bc.cprint(bc.RED, 'Y(es) / N(o): ')
        if input().lower()[0] == 'y':
            for job in jobs:
                if job_name is not None and job['NAME'] != job_name:
                    continue
                i = job['JOBID']
                subprocess.run(["scancel", str(i)],
                                 cwd=self.path,
                                 env=self.environ)
                logger.info(f"CANCEL JOB: {i}")
        else:
            return
