import re
from pathlib import Path

import pandas as pd


class Dir(type(Path())):
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, *args, **kwargs)

    def __init__(self):
        super().__init__()
    # pass


class File(type(Path())):
    # def __init__(self, *args, **kwargs):
    #     super().__init__(*args)
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, *args, **kwargs)



class PDBFile(File):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class CCP4File(File):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class MTZFile(File):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class PanDDADir(Dir):
    def __init__(self, pathlike):
        super().__init__(pathlike)


class PanDDAProcessedDatasetsDir(Dir):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class PanDDAProcessedDatasetDir(Dir):
    def __init__(self, pathlike):
        super().__init__(pathlike)


class PanDDAAnalysesDir(Dir):
    def __init__(self, pathlike):
        super().__init__(pathlike)

    @staticmethod
    def from_pandda_dir(pandda_dir: PanDDADir):
        return PanDDAAnalysesDir(pandda_dir / "analyses")


class PanDDAEventTablePath(File):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class PanDDAFSModel:
    def __init__(self, pandda_dir: PanDDADir):
        print("pandda dir: {}".format(pandda_dir))
        self.pandda_dir = pandda_dir
        self.analyses_dir = PanDDAAnalysesDir.from_pandda_dir(pandda_dir)
        self.processed_datasets_dir = PanDDAProcessedDatasetsDir.from_pandda_dir(pandda_dir)
        self.processed_datasets_dirs = {PanDDAProcessedDatasetDir.from_path(path)
                                        for path
                                        in self.processed_datasets_dir.glob("*")
                                        }

    @staticmethod
    def from_path(path: Path):
        pandda_dir = PanDDADir(path)
        return PanDDAFSModel(pandda_dir)

    @staticmethod
    def from_string(string: str):
        return PanDDAFSModel(PanDDADir(Path(string)))


class PanDDADtag(str):
    def __new__(cls, string):
        return str.__new__(cls, string)

    def __init__(self, string):
        str.__init__(string)


class PanDDAEventIdx(int):
    def __new__(cls, integer):
        return int.__new__(cls, integer)

    def __init__(self, integer):
        int.__init__(integer)


class PanDDAEventID:
    def __init__(self, dtag, event_idx):
        self.dtag = dtag
        self.event_idx = event_idx


class PanDDAModelPath(File):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class PanDDAEventMapPath(File):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class Event:
    def __init__(self, dtag, event_idx):
        self.dtag = PanDDADtag(dtag)
        self.event_idx = PanDDAEventIdx(event_idx)
        self.model_path = PanDDAModelPath()
        self.event_map_path = PanDDAEventMapPath()

    @staticmethod
    def from_record(record):
        dtag = record["dtag"]
        event_idx = record["event_idx"]

        return Event(dtag, event_idx)


class Command:
    pass


class RSCC(float):
    def __new__(cls, value):
        return float.__new__(cls, value)

    def __init__(self, value):
        float.__init__(value)

    @staticmethod
    def from_phenix_stdout(stdout: str):
        rscc_regex = "LIG[\s]+[^\s]+[\s]+[^\s]+[\s]+[^\s]+[\s]+([\s]+)"

        matches = re.findall(rscc_regex,
                             stdout,
                             )

        rscc = float(matches[0])

        return RSCC(rscc)


class GetPanDDAEventRSCCCommand(Command):
    def __init__(self,
                 pdb_path,
                 ccp4_path,
                 ):
        command = "{env}; phenix.real_space_correlation {pdb} {ccp4} detail=residue"
        env = "module load phenix"
        formatted_command = command.format(env=env,
                                           pdb=pdb_path,
                                           ccp4=ccp4_path,
                                           )
        self.formatted_command = formatted_command

    def __repr__(self):
        return str(self.formatted_command)

    @staticmethod
    def from_event(event: Event):
        return GetPanDDAEventRSCCCommand(event.model_path,
                                         event.event_map_path,
                                         )


class RSCCTable(pd.DataFrame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod
    def from_rsccs(rsccs):
        records = []
        for event_id, rscc in rsccs.items():
            record = {}
            record["dtag"] = event_id.dtag
            record["event_idx"] = event_id.event_idx
            record["rscc"] = rscc
            records.append(records)
        RSCCTable(records)


class PanDDAEventTable(pd.DataFrame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
