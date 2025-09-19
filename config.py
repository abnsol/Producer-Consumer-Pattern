import os

class Config:
    """Minimal config for testing analysis_tasks.py"""

    def __init__(self):
        self.plink_dir = "./data/mock_plink"
        self.data_dir = "./data"

    @classmethod
    def from_env(cls):
        """Mimic original API but only return what we need"""
        config = cls()
        config.plink_dir = os.getenv("PLINK_DIR", "./data/mock_plink")
        config.data_dir = os.getenv("DATA_DIR", "./data")
        return config

    @classmethod
    def from_args(cls, args):
        """Mimic original API (args parsing)"""
        config = cls()
        config.plink_dir = getattr(args, "plink_dir", "./data/mock_plink")
        config.data_dir = getattr(args, "data_dir", "./data")
        return config
