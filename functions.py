import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

import pandda_event_types

import joblib

from biopandas.pdb import PandasPdb


def map_parallel(f, iterable):
    results = joblib.Parallel(n_jobs=20,
                              verbose=50)(joblib.delayed(f)(i)
                                          for i
                                          in iterable)

    return results


def map_seriel_dict(f, dictionary):
    keys = list(dictionary.keys())
    values = list(dictionary.values())

    results = {key: f(value) for key, value in zip(keys, values)}

    return results


def map_parallel_dict(f, dictionary):
    keys = list(dictionary.keys())
    values = list(dictionary.values())

    results_list = joblib.Parallel(n_jobs=20,
                                   verbose=50)(joblib.delayed(f)(v)
                                               for v
                                               in values
                                               )

    results = {key: result for key, result in zip(keys, results_list)}

    return results


def execute(command):
    proc = subprocess.Popen(str(command),
                            shell=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            )
    stdout, stderr = proc.communicate()

    return str(stdout), str(stderr)


def get_event_rscc(event: pandda_event_types.Event):
    command = pandda_event_types.GetPanDDAEventRSCCCommand.from_event(event)

    stdout, stderr = execute(command)

    rscc = pandda_event_types.RSCC.from_phenix_stdout(stdout)

    return rscc


def get_pandda_events(pandda_fs_model: pandda_event_types.PanDDAFSModel):
    pandda_events_table = pandda_event_types.PanDDAEventTable.from_pandda_event_table_path(
        pandda_fs_model.analyses_dir.pandda_event_table_path)

    events = {}
    for index, row in pandda_events_table.iterrows():
        dtag = pandda_event_types.PanDDADtag(row["dtag"])
        event_idx = pandda_event_types.PanDDAEventIdx(row["event_idx"])
        event_id = pandda_event_types.PanDDAEventID(dtag, event_idx)

        events[event_id] = pandda_event_types.Event.from_record(row, pandda_fs_model)

    return events


# def get_closest_event_to_lig(events):
#
#     distances = {}
#     for event_id, event in events:
#         lig_com = get_lig_com()
#
#         lig_residue, distance_to_lig_com = get_distance_to_lig_com()
#
#         if event.dtag in distances:
#             distances[event.dtag][event.event_idx] = distance_to_lig_com
#         else:
#             distances[event.dtag] = {}
#             distances[event.dtag][event.event_idx] = distance_to_lig_com
#
#     closest_events_to_lig = {}
#     for dtag, events_dict in distances.items():
#         keys = list(events_dict.keys())
#         min_distance = min(keys, key=lambda key: events_dict[key])
#
#     return closest_events_to_lig

def get_closest_lig(event):
    event_model = PandasPdb().read_pdb(str(event.event_model_path))

    hetatm_table = event_model.df["HETATM"]
    ligands_table = hetatm_table[hetatm_table["residue_name"] == "LIG"]

    unqiue_ligand_residues = ligands_table["residue_number"]

    residue_tables = {}
    distances = {}
    for unqiue_ligand_residue in unqiue_ligand_residues:
        residue_table = ligands_table[ligands_table["residue_number"] == unqiue_ligand_residue]

        mean_coords = np.mean(residue_table[["x_coord", "y_coord", "z_coord"]].values,
                              axis=0,
                              )

        event_coords = event.coords

        distance = np.linalg.norm(mean_coords - event_coords)

        residue_tables[unqiue_ligand_residue] = residue_table
        distances[unqiue_ligand_residue] = distance

    closest_residue_table_key = min(residue_tables,
                                    key=lambda x: distances[x],
                                    )
    closest_residue_table = residue_tables[closest_residue_table_key]

    ppdb = PandasPdb()
    ppdb.df["HETATM"] = closest_residue_table

    return ppdb


def make_event_models(events, pandda_fs_model):
    for event_id, event in events.items():
        closest_lig = get_closest_lig(event)

        print(closest_lig)
        print(closest_lig.df["HETATM"])
        print(pandda_fs_model)
        print(pandda_fs_model.processed_datasets_dir)
        print("processed dataset dir: {}".format(pandda_fs_model.processed_datasets_dir[event.dtag]))

        event_model_path = pandda_fs_model.processed_datasets_dir[event.dtag].modelled_structures_dir / "{}.pdb".format(
                event.event_idx)

        closest_lig.to_pdb(str(event_model_path))


def get_rscc_table_from_pandda_dir(pandda_dir: Path):
    pandda_fs_model = pandda_event_types.PanDDAFSModel.from_path(pandda_dir)

    events = get_pandda_events(pandda_fs_model)

    make_event_models(events, pandda_fs_model)

    events_to_process = {}
    for event_id, event in events.items():
        if event.event_model_path.exists():
            events_to_process[event_id] = event
        else:
            print("\tCould not find path! No such file: {}".format(event.event_model_path))

    rsccs = map_seriel_dict(get_event_rscc,
                            events_to_process,
                            )
    # print(rsccs)

    rscc_table = pandda_event_types.RSCCTable.from_rsccs(rsccs)

    return rscc_table
