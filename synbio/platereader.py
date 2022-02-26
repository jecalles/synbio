from dataclasses import dataclass

from synbio.plates import Plate
from synbio.experiment import *


__all__ = [
    "PlateReaderExperiment", "PlateReaderData"

]

@dataclass
class PlateReaderData:
    def to_csv(self):
        pass

@dataclass
class PlateReaderExperiment(Experiment):
    """
    self.data : PlateReaderData
    self.plate: Plate
    """
    pass

    def subtract_background(self):
        pass

    def plot_results(self, print_flag: bool = False, plot_rc: dict = None):
        pass