import argparse
from pathlib import Path

from functions import get_rscc_table_from_pandda_dir


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--pandda_dir", required=True)
    parser.add_argument("-o", "--rscc_table_path", required=True)
    args = parser.parse_args()

    return args


def main():
    args = get_args()
    rscc_table = get_rscc_table_from_pandda_dir(Path(args.pandda_dir))
    print(rscc_table)
    rscc_table.to_csv(str(Path(args.rscc_table_path)))


if __name__ == "__main__":
    main()
