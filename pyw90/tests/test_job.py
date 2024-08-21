import pytest

from unittest.mock import patch, MagicMock
import pandas as pd
from io import StringIO
from pathlib import Path

from pyw90.lib.job import Job
from pyw90.lib.config import Config

@pytest.fixture
def job():
    config = Config(
        yaml_file=Path(__file__).parent / 'data/auto_w90_input.yaml',
        check=False,
    )
    return Job(config)

@patch('subprocess.run')
def test_get_jobs(mock_run, job):
    mock_run.return_value.stdout = "JOBID PARTITION NAME USER ST TIME CPUS NODES\n1234567 regular w90 test_user R 00:10:00 1 1"
    df = job.get_jobs()
    expected_df = pd.read_csv(StringIO(mock_run.return_value.stdout), sep='\s+')
    print(df)
    print(expected_df)
    pd.testing.assert_frame_equal(df, expected_df)

# @patch('time.sleep', return_value=None)
# @patch.object(Job, 'get_num_jobs', side_effect=[1, 1, 0])
# def test_check_run_until_stop(mock_get_num_jobs, mock_sleep, job, caplog):
#     job.check_run_until_stop()
#     assert 'Job <test_job> is running!' in caplog.text
#     assert 'Job <test_job> is still running.' in caplog.text

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