import subprocess
from pathlib import Path

import pandda_event_types


def execute(command):

    print(str(command))
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
    pandda_events_table = pandda_event_types.PanDDAEventTable.from_pandda_event_table_path(pandda_fs_model.analyses_dir.pandda_event_table_path)

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

    rsccs = {}
    for event_id, event in events.items():
        if event.event_model is not None:
            rscc = get_event_rscc(event)
            rsccs[event_id] = rscc

    rscc_table = pandda_event_types.RSCCTable.from_rsccs(rsccs)

    return rscc_table

