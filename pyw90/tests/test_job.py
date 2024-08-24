import pytest

from io import StringIO
from pathlib import Path
import logging

from pyw90.lib.job import Job

@pytest.fixture
def job():
    config = {
        'path': '.', 
        'win': 'wannier90.win', 
        'local': None, 
        'localrun':False, 
        'run': None,
        'usr_name': 'test_user', 
        'job_name': 'w90',
    }
    return Job(kwargs=config)

def test_get_jobs(mocker, job, caplog):
    """
    Test the get_jobs function of the Job class.

    Expected behavior
    -----------------
        The function should return a pandas DataFrame containing the output of the squeue command.
        here we are mocking the `subprocess.run` function (runs `squeue -o`) to return the queueing status.

    Parameters
    ----------
    mocker
        The mocker object from the pytest-mock library.
    job : Job
        An instance of the Job class.
    """

    mock_run = mocker.patch('subprocess.run')
    mock_run.return_value.stdout = "JOBID PARTITION NAME USER ST TIME CPUS NODES\n1234567 regular w90 test_user R 00:10:00 1 1"

    jobs = job.get_jobs()
    assert len(jobs) == 1
    assert 'w90' in str(jobs)

    caplog.set_level(logging.INFO)
    job.display_jobs(jobs)
    assert 'JOBID' in caplog.text
    assert '1234567' in caplog.text

def test_check_run_until_stop(mocker, job, caplog):
    caplog.set_level(logging.DEBUG)
    mocker.patch('time.sleep', return_value=None)
    mocker.patch.object(Job, 'get_num_jobs', side_effect=[1, 1, 1, 0, 1])
    spy_get_num_jobs = mocker.spy(job, 'get_num_jobs')

    job.check_run_until_stop()

    assert f'Job <{job.job_name}> is running!' in caplog.text
    assert f'Job <{job.job_name}> is still running.' in caplog.text
    # we expect the job to be stopped when the number of jobs is 0
    assert spy_get_num_jobs.call_count == 4

# @patch.object(Job, 'get_jobs', return_value=pd.DataFrame({'JOBID': [1234567]}))
# def test_get_num_jobs(mock_get_jobs, job):
#     assert job.get_num_jobs() == 1
#     job.local = True
#     job.p = MagicMock()
#     job.p.poll.return_value = None
#     assert job.get_num_jobs() == 1
#     job.p.poll.return_value = 0
#     assert job.get_num_jobs() == 0

# @patch('subprocess.Popen')
# @patch('subprocess.run')
# def test_submit(mock_run, mock_popen, job):
#     job.local = True
#     job.submit()
#     mock_popen.assert_called_once_with(
#         job.config.localrun.split(),
#         cwd=job.config.path,
#         env=job.config.environ,
#         start_new_session=True
#     )
#     job.local = False
#     mock_run.reset_mock()
#     mock_popen.reset_mock()
#     with patch.object(Job, 'get_jobs', return_value=pd.DataFrame()):
#         job.submit()
#         mock_run.assert_called_once_with(
#             ["sbatch", job.config.run],
#             cwd=job.config.path,
#             env=job.config.environ
#         )

# @patch('builtins.input', return_value='y')
# @patch('subprocess.run')
# @patch.object(Job, 'get_jobs', return_value=pd.DataFrame({'JOBID': [1234567]}))
# def test_cancel(mock_get_jobs, mock_run, mock_input, job):
#     with patch('pyw90.utility.utility.bc.cprint') as mock_cprint:
#         job.cancel()
#         mock_run.assert_called_once_with(
#             ["scancel", '1234567'],
#             cwd=job.config.path,
#             env=job.config.environ
#         )
#         mock_cprint.assert_called()