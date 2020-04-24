import subprocess
from pathlib import Path

import pandda_event_types

import joblib


def map_parallel(f, iterable):
    results = joblib.Parallel(n_jobs=20,
                              verbose=50)(joblib.delayed(f)(i)
                                            for i
                                            in iterable)

    return results


def map_parallel_dict(f, dictionary):
    keys = list(dictionary.keys())
    values = list(dictionary.values())

    results_list = joblib.Parallel(n_jobs=20,
                                   verbosity=50)(joblib.delayed(f)(v)
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

    print(stdout)
    print(stderr)

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


def get_rscc_table_from_pandda_dir(pandda_dir: Path):
    pandda_fs_model = pandda_event_types.PanDDAFSModel.from_path(pandda_dir)

    events = get_pandda_events(pandda_fs_model)

    events_to_process = {}
    for event_id, event in events.items():
        if event.event_model_path.exists():
            events_to_process[event_id] = event
        else:
            print("\tCould not find path! No such file: {}".format(event.event_model_path))

    rsccs = map_parallel_dict(get_event_rscc,
                              events_to_process,
                              )

    rscc_table = pandda_event_types.RSCCTable.from_rsccs(rsccs)

    return rscc_table
