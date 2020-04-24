import re
from pathlib import Path

import pandas as pd


class Dir(type(Path())):
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, args[0])

    def __init__(self, pathlike):
        super().__init__()


class File(type(Path())):
    def __init__(self, pathlike):
        super().__init__()

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

    @staticmethod
    def from_pandda_dir(pandda_dir: PanDDADir):
        return PanDDAProcessedDatasetsDir(pandda_dir / "processed_datasets")


class PanDDAModelledStucturesDir(Dir):
    def __init__(self, pathlike, modelled_structure_path):
        super().__init__(pathlike)
        self.event_model_path = modelled_structure_path

    @staticmethod
    def from_dataset_dir(pathlike, dtag):
        modelled_structure_path = PanDDAModelPath(pathlike / "{dtag}-pandda-model.pdb".format(dtag=dtag))
        return PanDDAModelledStucturesDir(pathlike,
                                          modelled_structure_path,
                                          )


class PanDDAProcessedDatasetDir(Dir):
    def __init__(self, pathlike, modelled_structures_dir, event_maps, model_path):
        super().__init__(pathlike)
        self.modelled_structures_dir = modelled_structures_dir
        self.event_maps = event_maps
        self.model_path = model_path

    @staticmethod
    def from_path(path):
        dtag = path.name
        modelled_structures_dir = PanDDAModelledStucturesDir.from_dataset_dir(path / "modelled_structures",
                                                                              dtag,
                                                                              )
        event_map_paths = list(path.glob("{}-event*".format(dtag)))
        if len(event_map_paths) > 0:
            event_maps_dict = {re.findall("event_([0-9]+)_", str(evnet_map_path))[0]: evnet_map_path
                               for evnet_map_path
                               in event_map_paths
                               }
            event_maps = {event_idx: PanDDAEventMapPath(event_map_path)
                          for event_idx, event_map_path
                          in event_maps_dict.items()
                          }
        else:
            event_maps = {}
        model_path = PanDDAModelPath(path / "{}-pandda-input.pdb".format(dtag))

        return PanDDAProcessedDatasetDir(path,
                                         modelled_structures_dir,
                                         event_maps,
                                         model_path,
                                         )


class PanDDAAnalysesDir(Dir):
    def __init__(self, pathlike):
        super().__init__(pathlike)
        self.pandda_event_table_path = PanDDAEventTablePath(pathlike / "pandda_inspect_events.csv")

    @staticmethod
    def from_pandda_dir(pandda_dir: PanDDADir):
        return PanDDAAnalysesDir(pandda_dir / "analyses")


class PanDDAEventTablePath(File):
    def __init__(self, pathlike):
        super().__init__(pathlike)


class PanDDAFSModel:
    def __init__(self, pandda_dir: PanDDADir):
        print("pandda dir: {}".format(pandda_dir))
        self.pandda_dir = pandda_dir
        self.analyses_dir = PanDDAAnalysesDir.from_pandda_dir(pandda_dir)
        self.processed_datasets_dir = PanDDAProcessedDatasetsDir.from_pandda_dir(pandda_dir)
        self.processed_datasets_dirs = {path.name: PanDDAProcessedDatasetDir.from_path(path)
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


class PanDDAModelPath(PDBFile):
    def __init__(self, pathlike):
        super().__init__(pathlike)


class PanDDAEventMapPath(CCP4File):
    def __init__(self, pathlike):
        super().__init__(pathlike)


class Event:
    def __init__(self, dtag, event_idx, event_dir, model_path, event_map_path, event_model_path):
        self.dtag = dtag
        self.event_idx = event_idx

        self.event_dir = event_dir
        self.model_path = model_path
        self.event_map_path = event_map_path

        self.event_model_path = event_model_path

    @staticmethod
    def from_record(record, pandda_fs_model: PanDDAFSModel):
        dtag = PanDDADtag(record["dtag"])
        event_idx = PanDDAEventIdx(record["event_idx"])
        event_dir = pandda_fs_model.processed_datasets_dirs[dtag]

        model_path = event_dir.model_path
        event_map_path = event_dir.event_maps[str(event_idx)]

        event_model_path = event_dir.modelled_structures_dir.event_model_path

        return Event(dtag,
                     event_idx,
                     event_dir,
                     model_path,
                     event_map_path,
                     event_model_path,
                     )


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
        print("Matches!")
        print(matches)

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
        return GetPanDDAEventRSCCCommand(event.event_model_path,
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

class PanDDAEventModel(PDBFile):
    def __init__(self, pathlike):
        super().__init__(pathlike)

    @staticmethod
    def from_modelled_structures_dir(modelled_structures_dir, dtag):
        return PanDDAEventModel(modelled_structures_dir / "{dtag}-pandda-model.pdb".format(dtag=dtag))


class PanDDAEventTable(pd.DataFrame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod
    def from_pandda_event_table_path(path: PanDDAEventTablePath):
        table = pd.read_csv(str(path))
        return PanDDAEventTable(table)
