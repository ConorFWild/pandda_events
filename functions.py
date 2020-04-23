import subprocess
from pathlib import Path

import types


def execute(command):
    proc = subprocess.Popen(str(command),
                            shell=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            )
    stdout, stderr = proc.communicate()

    return stdout, stderr


def get_event_rscc(event: types.Event):
    command = types.GetPanDDAEventRSCCCommand.from_event(event)

    stdout, stderr = execute(command)

    rscc = types.RSCC.from_phenix_stdout(stdout)

    return rscc


def get_pandda_events(pandda_fs_model: types.PanDDAFSModel):
    pandda_events_table = types.PanDDAEventTable(pandda_fs_model.analyses_dir.pandda_event_table_path)

    events = {}
    for index, row in pandda_events_table.iterrows():
        dtag = types.PanDDADtag(row["dtag"])
        event_idx = types.PanDDAEventIdx(row["event_idx"])
        event_id = types.PanDDAEventID(dtag, event_idx)
        events[event_id] = types.Event.from_record(row)

    return events


def get_rscc_table_from_pandda_dir(pandda_dir: Path):
    pandda_fs_model = types.PanDDAFSModel.from_path(pandda_dir)

    events = get_pandda_events(pandda_fs_model)

    rsccs = {}
    for event_id, event in events.items():
        rscc = get_event_rscc(event)
        rsccs[event_id] = rscc

    rscc_table = types.RSCCTable.from_rsccs(rsccs)

    return rscc_table

