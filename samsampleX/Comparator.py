from logging import info

import pysam
import tempfile
import pandas as pd

from samsampleX.Loader import Loader


class Comparator:

    def __init__(self, bam_left_path: str, bam_right_path: str, out: str) -> None:
        """ """
        info(f"Comparator - Initialize Comparator")

        self.fields_of_interest = [
            "name",
            "ref_name",
            "ref_pos",
            "next_ref_name",
            "next_ref_pos",
            "length",
        ]

        bam_left = Loader(bam_left_path)
        bam_right = Loader(bam_right_path)

        # Get all query_names from the smaller BAM (right)
        bam_right_query_names = set()
        for read in bam_right.fetch():
            bam_right_query_names.add(read.query_name)

        with tempfile.NamedTemporaryFile(mode="w+") as tmpfile:

            info(f"Comparator - Create temporary file for read names")

            for line in bam_right_query_names:
                tmpfile.write(f"{line}\n")
            tmpfile.flush()

            bam_left_str = pysam.view(str(bam_left_path), "-N", tmpfile.name)
            bam_right_str = pysam.view(str(bam_right_path), "-N", tmpfile.name)

            bam_right_df = self.get_reads(bam_right_str)
            bam_left_df = self.get_reads(bam_left_str)

        info(f"Comparator - Inner join BAMs")

        self.suffixes = ("_l", "_r")
        merged_df = pd.merge(
            left=bam_left_df,
            right=bam_right_df,
            how="right",
            on="name",
            suffixes=self.suffixes,
        )

        merged_df = self.format_df(merged_df)
        merged_df.to_csv(out, sep="\t", index=False)

    def get_reads(self, read_names_str: str) -> pd.DataFrame:

        info(f"Comparator - Get reads from BAM")

        reads = []
        for line in read_names_str.splitlines():
            s = line.split("\t")
            reads.append(
                {
                    "name": s[0],
                    "ref_name": s[2],
                    "ref_pos": s[3],
                    "next_ref_name": s[6],
                    "next_ref_pos": s[7],
                    "length": s[8],
                }
            )
        return pd.DataFrame.from_records(reads)

    def format_df(self, df: pd.DataFrame) -> pd.DataFrame:

        column_order = ["name"]
        for col in [c for c in self.fields_of_interest if c != "name"]:
            for suf in self.suffixes:
                column_order.append(f"{col}{suf}")

        df = df[column_order].sort_values(by="name")

        df.loc[df["next_ref_name_l"] == "=", "next_ref_name_l"] = df["ref_name_l"]
        df.loc[df["next_ref_name_r"] == "=", "next_ref_name_r"] = df["ref_name_r"]

        return df
