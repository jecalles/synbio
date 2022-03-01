import re
from copy import deepcopy
from datetime import date
from itertools import islice
from math import ceil, prod
from typing import Dict, List, Tuple
from dataclasses import dataclass

import numpy as np
import pandas as pd
import pint

from synbio.plates import *
from synbio.reagents import *
from synbio.experiment import *
from synbio.platereader import *

"""
TODO:
[x] Sort out Experiment inheritance tree (branch Experiment)
[x] Sort out PlateReaderExperiment
[] Sort out EchoExperiment
"""
__all__ = [
    # dataclasses
    "EchoProtocol", "EchoExperiment"
]

@dataclass
class EchoProtocol:
    src_plate: Plate = None
    dest_plate: Plate = None


class EchoExperiment(PlateReaderExperiment):
    """
    self.src_plate: Plate = None
    self.dest_plate: Plate = None
    """
    def __init__(
        self, protocol: EchoProtocol = None,
        *plate_args, **plate_kwargs,
    ):
        super().__init__(*plate_args, **plate_kwargs)

        if protocol is None:
            protocol =  self.generate_protocol()
        self.protocol = protocol

    def generate_protocol(self):
        """TODO: this"""
        pass



"""QUARANTINE ZONE

Old code

Manifest = [
    # functions
    "assign_reagents_to_wells", "update_plates",
    "make_echo_prot", "make_experiment"
]
"""
def assign_reagents_to_wells(
        conditions: List[Condition],
        reagents: Dict[Reagent, pint.Quantity],
        source_plate: Plate,
        destination_plate: Plate
) -> Tuple[List[Well], List[Well]]:
    # destination to wells
    dvols = destination_plate.well_volumes

    dest_wells = [
        Well(
            name=f"dest_{cond.name}_{i}",
            content=cond.content,
            max_vol=dvols['max_vol'], dead_vol=dvols['dead_vol'],

            vol=cond.volume
        )
        for cond in conditions
        for i in range(cond.replicates)
    ]

    # source to wells
    svols = source_plate.well_volumes
    num_wells = {
        reg: ceil(quant / svols['working_vol'])
        # number of working volumes required
        for reg, quant in reagents.items()
    }

    well_vols = {
        reg: [
            Well(
                name=f"source_{reg.name}_{i + 1}", content=reg,
                max_vol=(mv := svols['max_vol']),
                dead_vol=(dv := svols['dead_vol']),

                vol=(mv if i != (n - 1)
                     else reagents[reg] - svols['working_vol'] * (
                            n - 1) + dv)

            ) for i in range(n)
        ] for reg, n in num_wells.items()
    }

    source_wells = [w for r, wells in well_vols.items() for w in wells]

    return (source_wells, dest_wells)


def update_plates(
        plate: Plate, wells: List[Well], offset: int = 0
) -> Dict[str, str]:
    """
    Modifies plate in place by assigning wells to plate.array. Returns a python dict
    representing a "name_map" to be used


    TODO: change "offset" to an optional well_dict parameter
    """
    if (n_assign := len(wells)) > (n_plate := prod(plate.shape)) - offset:
        raise ValueError(
            f"cannot assign {n_assign} wells to {n_plate} well plate!")

    well_names = islice(plate.dict, offset, None)

    well_dict = {
        name: well
        for name, well in zip(well_names, wells)
    }
    plate.dict = well_dict

    return {name: well.name for name, well in well_dict.items()}


def make_echo_prot(
        source: Plate, destination: Plate,
        filename: str = None
) -> pd.DataFrame:
    # first, make copies of src_plate and destinatino plate
    src_plate = deepcopy(source)
    dest_plate = deepcopy(destination)

    src_wells = []
    dest_wells = []
    xfer_vols = []

    # find src_plate wells w/ reagent
    reagents = {
        rgnt
        for well in dest_plate.full_wells.values()
        for rgnt in well.recipe.keys()
    }
    src_rgnt_map = {
        rgnt: [
            well_loc
            for well_loc, well in src_plate.full_wells.items()
            if well.content == rgnt
        ] for rgnt in reagents
    }

    for dest_loc, dest_well in dest_plate.full_wells.items():
        for rgnt, amt in dest_well.recipe.items():
            # get list of source_well locations for given rgnt
            rgnt_src_locs = src_rgnt_map[rgnt]
            rgnt_src_locs.sort(reverse=True)

            # pop off wells until dest_well is satisfied
            sat = False
            while sat is not True:
                src_loc = rgnt_src_locs.pop()
                src_well = src_plate.dict[src_loc]

                dest_wells.append(dest_loc)
                src_wells.append(src_loc)

                if (avail := src_well.available_vol) >= amt:  # there's enough volume
                    sat = True
                    # do transfer
                    xfer = amt
                    src_well.volume -= amt
                    # return well to rgnt_src_locs
                    rgnt_src_locs.append(src_loc)
                else:  # not enough volume
                    # transfer what you can and adjust amt required
                    xfer = avail
                    amt -= xfer

                # do transfer
                xfer_vols.append(xfer)

            # update source_reagent_map[rgnt]
            src_rgnt_map[rgnt] = rgnt_src_locs

    # massage xfer volumes
    xfer_vols = [
        np.round(q.to('nL').magnitude, decimals=4)
        for q in xfer_vols
    ]
    num_xfer = len(xfer_vols)
    # write protocol to pd.DataFrame
    protocol = {
        'Source plate name': [1 for _ in range(num_xfer)],
        'Source well': src_wells,
        'Destination plate name': [1 for _ in range(num_xfer)],
        'Destination well': dest_wells,
        'XferVol': xfer_vols
    }

    def src_sort_key(col: pd.Series) -> List[Tuple[str, int]]:
        r = re.compile(r"\d")
        searches = [r.search(loc) for loc in col]
        ixs = [hit.start() for hit in searches]
        return [
            (loc[:ix], int(loc[ix:]))
            for loc, ix in zip(col, ixs)
        ]

    df = (pd.DataFrame
          .from_dict(protocol)
          .sort_values(by="Source well", key=src_sort_key))
    df = df[df['XferVol'] > 0]

    # optionally write file to filename
    if filename is not None:
        df.to_csv(filename)

    return df

"""
TODO: make this a class method for EchoExperiment; make standalone function 
that calls EchoExperiment under the hood.

"""
def make_experiment(
        name: str,
        conditions: List[Condition],
        src_plate: Plate, dest_plate: Plate,
        src_plate_offset: int = 0
) -> EchoExperiment:
    reagents = calc_reagents(conditions)

    src_wells, dest_wells = assign_reagents_to_wells(conditions, reagents,
                                                     src_plate, dest_plate)

    name_map = update_plates(dest_plate, dest_wells)
    _ = update_plates(src_plate, src_wells, src_plate_offset)

    # TODO: make this return an EchoExperiment object
    return Experiment(
        name=name, source_plate=src_plate, dest_plate=dest_plate,
        conditions=conditions, name_map=name_map, data=None, date=date.today()
    )