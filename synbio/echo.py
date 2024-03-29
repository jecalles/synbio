from dataclasses import dataclass, field
from typing import Callable, Dict, Tuple, List
import re

import numpy as np
import pandas as pd
import pint

from synbio.units import unit_registry as u
from synbio.platereader import *
from synbio.plates import *
from synbio.reagents import Reagent

"""
TODO:
[x] Sort out Experiment inheritance tree (branch Experiment)
[x] Sort out PlateReaderExperiment
[x] Sort out EchoExperiment
[x] Sort out EchoProtocol.generate_protocol()
[ ] Sort out EchoProtocol.check_protocol()
[ ] Sort out EchoExperiment.simulate()
"""
__all__ = [
    # dataclasses
    "EchoProtocol", "EchoExperiment"
]


@dataclass
class EchoProtocol:
    src_plate: Plate = make_384_ldv_well()
    dest_plate: Plate = make_96_well()
    dataframe: pd.DataFrame = field(default=None, repr=False)

    def __post_init__(self):
        if (self.dataframe is None) and (len(self.dest_plate.full_wells) > 0):
                self.dataframe = self.generate_protocol()


    def generate_protocol(self) -> pd.DataFrame:
        from itertools import count
        if not self.check_protocol():
            raise ValueError("source plate doesn't have enough material to "
                             "populate destination plate")
        # first, store well volumes for source plate
        src_well_vols = {
            loc: well.volume
            for loc, well in self.src_plate.dict.items()
        }
        # initialize containers for columns of protocol
        src_wells = []
        dest_wells = []
        xfer_vols = []

        # loop over destination plate; assign transfer events
        for dest_loc, dest_well in self.dest_plate.full_wells.items():
            for rgt, rgt_vol in dest_well.reagent_volumes.items():
                # satisfy volume reqs for each reagent in dest_plate
                sat = False
                while not sat:
                    src_wells_for_rgt = [
                        well for well in self.src_plate.wells_by_content[rgt]
                        if well.available_vol > 0 * u.uL
                    ]
                    try:
                        src_well = src_wells_for_rgt.pop()
                    except:
                        raise RuntimeError("could not complete assign")
                    if (avail_vol:= src_well.available_vol) >= rgt_vol:
                        # transfer rgt_vol from src_well to dest_well
                        xfer_vol = rgt_vol

                        # mark reagent as satisfied
                        sat = True
                    else:
                        # transfer what you can and adjust rgt_vol as necessary
                        xfer_vol = avail_vol
                        rgt_vol -= xfer_vol
                        # don't mark reagent sat
                        pass

                    # adjust src_well vol
                    src_well.volume -= xfer_vol
                    # append transfer info to storage variables
                    src_wells.append(src_well.location)
                    dest_wells.append(dest_well.location)
                    xfer_vols.append(xfer_vol)

        # round transfer volumes
        xfer_vols = [
            np.round(q.to('nL').magnitude, decimals=4)
            for q in xfer_vols
        ]
        # write protocol to dict
        num_xfer = len(xfer_vols)
        dict_protocol = {
            'Source plate name': [1 for _ in range(num_xfer)],
            'Source well': src_wells,
            'Destination plate name': [1 for _ in range(num_xfer)],
            'Destination well': dest_wells,
            'XferVol': xfer_vols
        }

        # convert protocol to dataframe
        def src_sort_key(col: pd.Series) -> List[Tuple[str, int]]:
            r = re.compile(r"\d")
            searches = [r.search(loc) for loc in col]
            ixs = [hit.start() for hit in searches]
            return [
                (loc[:ix], int(loc[ix:]))
                for loc, ix in zip(col, ixs)
            ]

        df = (pd.DataFrame
              .from_dict(dict_protocol)
              .sort_values(by="Source well", key=src_sort_key))
        df = df[df['XferVol'] > 0]  # ignore tranfer instructions w/ zero vol

        # finally, reset source plate well volumes
        for loc, well in self.src_plate.dict.items():
            well.volume = src_well_vols[loc]

        return df

    def check_protocol(self) -> bool:
        def check_plates(self) -> bool:
            full_src_wells = len(self.src_plate.full_wells)
            full_dest_wells = len(self.dest_plate.full_wells)

            if full_src_wells == 0:
                if full_dest_wells == 0:
                    return True
                else:
                    return False
            else:
                avail_vols = {
                    rgt: sum(
                        well.available_vol for well in well_list
                    ) for rgt, well_list in self.src_plate.wells_by_content.items()
                }
                req_vols = self.dest_plate.reagent_volumes

                return all(
                    avail_vols[rgt] >= req_vols[rgt]
                    for rgt in self.dest_plate.reagents
                )
        def check_dataframe(self) -> bool:
            # dereference attributes for convenience
            src = self.src_plate
            dest = self.dest_plate
            df = self.dataframe

            # check dest_plate is satisfied by dataframe
            for content, well_list in dest.wells_by_content.items():
                well_xfer_vols: List[Tuple[pint.Quantity, pint.Quantity]] = [
                    (well.volume,
                     df[df["Destination well"] == well.location]["XferVol"].sum() * u.nL)
                    for well in well_list
                ]
                close_list = [np.isclose(well_vol, xfer_vol) for
                              well_vol, xfer_vol in
                              well_xfer_vols]
                if not all(close_list):
                    raise ValueError(f"dest_well: {content} has wells not "
                                     f"close to spec")

            # check src_plate is satisfied by dataframe


    def to_csv(self, filename: str) -> None:
        self.protocol.to_csv(filename)


class EchoExperiment(PlateReaderExperiment):
    def __init__(
            self, protocol: EchoProtocol = EchoProtocol(),
            *plate_args, **plate_kwargs,
    ):
        super().__init__(*plate_args, **plate_kwargs)
        self.protocol = protocol

    @classmethod
    def from_plate_map(
            cls, src_plate_filepath: str, dest_plate_filepath: str,
            cond_name_map: Callable[[str, Plate], PlateReaderCondition],
            src_plate: Plate = make_384_ldv_well()
    ) -> "OwnType":
        # load up Experiment
        conditions, dest_plate = load_plate_map(dest_plate_filepath,
                                                cond_name_map)
        src_plate.load_plate_map(src_plate_filepath)
        exp = cls(conditions=conditions)

        # get source wells by reagent
        src_wells_by_reagent = src_plate.wells_by_content

        # calculate number of wells and volume required per reagent
        well_vols = exp.calc_src_wells()

        # iterate over wells by reagent; check num_wells; update volume
        for rgt, well_list in src_wells_by_reagent.items():
            num_wells, vol = well_vols[rgt]

            # check num_wells matches actual number of wells
            if len(well_list) != num_wells:
                raise ValueError("something went wrong")

            # update volume
            for well in well_list:
                well.volume = vol

        exp.protocol = EchoProtocol(src_plate=src_plate, dest_plate=dest_plate)
        return exp

    def calc_src_wells(
            self, plate: Plate = None,
            buffer_vol: pint.Quantity = 0.1 * u.uL
    ) -> Dict[Reagent, Tuple[int, pint.Quantity]]:
        if plate is None:
            plate = self.protocol.src_plate

        working_vol = plate.well_volumes["working_vol"]
        rgt_vols = {
            rgt: vol + buffer_vol
            for rgt, vol in self.reagent_volumes.items()
        }

        num_wells = {
            rg: np.ceil((vol/ working_vol).magnitude)
            for rg, vol in rgt_vols.items()
        }
        dead_vol = plate.well_volumes["dead_vol"]
        vol_per_well = {
            rg: np.round(
                (vol / num_wells[rg]) + dead_vol,
                decimals=2
           ) for rg, vol in rgt_vols.items()
        }
        return {
            rg: (num_wells[rg], vol_per_well[rg])
            for rg in self.reagents.values()
        }

    def simulate(self) -> bool:
        """ The idea is to get some sense of what an experiment is going to
        look like"""
        return self.protocol.check_protocol()