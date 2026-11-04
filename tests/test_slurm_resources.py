"""Check scheduler options without submitting any jobs."""
import importlib.util
from pathlib import Path
import sys
import unittest

PROFILE = Path(__file__).resolve().parents[1] / 'workflow/profiles/slurm'
HAS_SNAKEMAKE = importlib.util.find_spec('snakemake') is not None
if HAS_SNAKEMAKE:
    sys.path.insert(0, str(PROFILE))
    import slurm_utils


@unittest.skipUnless(HAS_SNAKEMAKE, 'requires Snakemake')
class ResourceTests(unittest.TestCase):
    def options(self, resources):
        return slurm_utils.convert_job_properties({'threads': 4, 'resources': resources})

    def test_gpu_is_requested(self):
        self.assertEqual(self.options({'gpu': 1})['gres'], 'gpu:1')
        self.assertEqual(self.options({'gpu': 2, 'gpu_model': 'a100'})['gres'], 'gpu:a100:2')

    def test_scratch_is_preserved(self):
        self.assertEqual(self.options({'gpu': 1, 'slurm': 'gres=lscratch:100'})['gres'], 'lscratch:100,gpu:1')
        self.assertEqual(self.options({'gres': 'lscratch:50'})['gres'], 'lscratch:50')

    def test_explicit_gpu_is_authoritative(self):
        self.assertEqual(self.options({'gpu': 1, 'slurm': 'gres=gpu:a100:2,lscratch:50'})['gres'], 'gpu:a100:2,lscratch:50')
        self.assertNotIn('gres', self.options({'gpu': 1, 'slurm': 'gpus=2'}))

    def test_cpu_job_stays_cpu(self):
        self.assertNotIn('gres', self.options({}))
        self.assertEqual(self.options({})['cpus-per-task'], 4)
        with self.assertRaises(ValueError):
            self.options({'gpu': -1})
