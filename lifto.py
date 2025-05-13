import pandas as pd
from liftover import get_lifter
import logging
import sys
from typing import Literal
import os
import argparse
from natsort import natsort_keygen

# logging definition -> Maybe this goes to main(?)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
# añadimos un handler para stdout
c_handler = logging.StreamHandler(sys.stdout)
c_format = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
c_handler.setFormatter(c_format)
logger.addHandler(c_handler)


ASSEMBLY = {"hg19": "GRCh37", "hg38": "GRCh38"}


def check_coords(pre_data: pd.DataFrame) -> bool:
    """
    Check if coordinate system is 1-based o 0-based
    Based on SNV intervals

    1-based
    - SNVs --> START = END
    - INS --> START = END
    - DEL --> START = END - len(REF) + 1

    0-based
    - SNVs --> START = END - 1
    - INS --> START = END
    - DEL --> START = END - len(REF)
    """
    nuc = ["A", "C", "G", "T"]
    snvs = pre_data[pre_data["REF"].isin(nuc) & pre_data["ALT"].isin(nuc)]
    if all(snvs["POS"] == snvs["END"]):
        logger.debug("1-based")
        return True
    else:
        logger.debug("0-based")
        return False


def check_conversion(start, end, ref, alt, logger):
    """
    Check new POS and END agree with REF/ALT lengths and
    coordinate system type.
    TODO
    """
    pass


def enforce_col_types(data: pd.DataFrame):
    data["CHROM"] = data["CHROM"].astype(str)
    data["POS"] = data["POS"].astype(int)


def check_columns_present(
    df: pd.DataFrame, required_cols: list[str], required_are: int = 2, column_mapping: list[str] = []
) -> list[str]:
    """
    Check columns are present in dataframe
    Args:
    - df
    - required_columns DEFAULT: CHROM, POS, REF, ALT
    - required_are: posición de las columnas esenciales
    las que hacen saltar error
    - column_mapping: map each column to the required columns

    Raises:
    - IndexError: Missing columns

    Return:
    - Columns present in dataframe
    """
    if column_mapping:
        # map columns to required columns
        # transform list to dict splitting by ":" k:v and
        column_mapping = [x.split(":") for x in column_mapping]
        column_mapping = {x[0]:x[1] for x in column_mapping}
        df.rename(columns=column_mapping, inplace=True)
        logger.debug(f"Renamed columns {column_mapping} to {required_cols}")

    essentials = ["CHROM", "POS", "REF", "ALT"]
    missing = [c for c in required_cols if c not in df.columns]
    present = [c for c in required_cols if c in df.columns]
    if any([x in missing for x in essentials]):
        raise IndexError(
            f"The following columns are missing \
{', '.join(missing)}"
        )
    elif missing:
        logger.warning(f"""The following columns are missing {missing}.
The liftover coordinates will not be checked against REF-ALT lengths""")
    enforce_col_types(df)
    logger.debug(present)
    return present


def lift_coords(
    pre_data: pd.DataFrame, lifter
) -> tuple[dict[tuple:tuple], list[tuple]]:
    """
    Iterates over a dataframe and lift the coordinates using `lifter` object

    Returns:
    - tuple: Dictionary of positions and list of errors (empty if None)
    """
    post_dict = dict()
    error = list()
    n = 0

    is_ref_in_col = True if "REF" in pre_data.columns else False
    is_alt_in_col = True if "ALT" in pre_data.columns else False
    is_end_in_cols = True if "END" in pre_data.columns else False

    for _, row in pre_data.iterrows():
        keys = tuple(row.loc[pre_data.columns].to_list())
        try:
            chrom, pos = lifter.convert_coordinate(str(row.CHROM), row.POS)[0][:2]
            if is_end_in_cols:
                end = lifter.convert_coordinate(str(row.CHROM), row.END)[0][1]
            else:
                end = 0
            if is_ref_in_col:
                if is_alt_in_col:
                    post_dict[keys] = (
                        (chrom, pos, end, row.REF, row.ALT)
                        if is_end_in_cols
                        else (chrom, pos, row.REF, row.ALT)
                    )
                else:
                    post_dict[keys] = (
                        (chrom, pos, end, row.REF)
                        if is_end_in_cols
                        else (chrom, pos, row.REF)
                    )
            elif is_alt_in_col:
                post_dict[keys] = (
                    (chrom, pos, end, row.ALT)
                    if is_end_in_cols
                    else (chrom, pos, row.ALT)
                )
            else:
                post_dict[keys] = (chrom, pos, end) if end != 0 else (chrom, pos)
            # logger.debug(f'{chrom}:{row.POS}- new:{pos}')
        except IndexError:
            error.append(tuple(keys))
            n = n + 1
        except KeyError:
            error.append(tuple(keys))
            logger.debug(f"Error processing {row.CHROM}, jumping next...")
            n = n + 1
    # logger.debug(post_dict)
    return post_dict, error


def save_errors(error: list, output_file: str):
    """Save not lifted variants.
    Args:
    - error: list of errors
    - output_file: path to save the errors"""
    # create a log folder in the current directory
    cwd = os.getcwd()
    try:
        os.mkdir(cwd + "/_log/")
        logger.info(f"Directory {cwd + '/_log/'} created")
    except FileExistsError:
        logger.info(f"Directory {cwd + '/_log/'} already exists")
    # rename error file to .err
    # extract file name from path
    err_file = os.path.basename(output_file)
    err_path = os.path.join(cwd + "/_log/", f"{err_file}.err")
    # check if file exists
    if os.path.isfile(err_path):
        logger.info(f"File {err_path} already exists, overwriting...")
        os.remove(err_path)
    logger.info(
        f"Failed a total of {len(error)}, saving file to {err_path}"
    )
    with open(err_path, "w") as f:
        try:
            # write to file the error list
            f.write(f"{','.join(error)}\n")
        except TypeError as e:
            logger.error(f"Error writing to error file: {e}")
        for value in error:
            f.write(f"{','.join(map(str, value))}\n")
    logger.info(f"Error file saved at {err_path}")


def create_dir(f):
    directory = os.path.dirname(f)
    if not os.path.exists(directory):
        os.makedirs(directory)


def parse_format(format_col: str):
    """Given a FORMAT string [GT:GQ:...]
    creates a string to parse the header"""
    format_header = ""
    format_col_list = format_col.split(":")
    format_dict = {
        "GT": '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
        "GQ": '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n',
        "DP": '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n',
        "HQ": '##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">\n',
    }
    for param in format_col_list:
        try:
            format_header = format_header + format_dict[param]
            logger.info(f"Adding {param} to header")
        except KeyError as e:
            logger.warning(
                f"{param} not found in format dictionary, including a generic..."
            )
            generic = f'##FORMAT=<ID={param},Number=.,Type=String,Description="{param} from generic">\n'
            format_header = format_header + generic
    return format_header


def reorder_lifted_coords(
    vcf_df: pd.DataFrame, new_assembly: str, old_assembly: str
) -> pd.DataFrame:
    """Helper function to undo the mess coming y
    from the lifted dataframe"""
    # fix order (as CHROM_{new} is the target col, not CHROM)
    colist = vcf_df.columns.to_list()
    ch, ch_n = colist.index(f"CHROM_{new_assembly}"), colist.index("CHROM")
    colist[ch], colist[ch_n] = colist[ch_n], colist[ch]
    # same with POS <-> POS_{new}
    # logger.debug(vcf_df['POS'].head())
    # logger.debug(vcf_df[f'POS_{new_assembly}'].head())
    pos, pos_n = colist.index(f"POS_{new_assembly}"), colist.index("POS")
    colist[pos], colist[pos_n] = colist[pos_n], colist[pos]

    # assign new col order to df
    vcf_df = vcf_df.reindex(columns=colist)

    vcf_df["INFO"] = (
        f"{ASSEMBLY[old_assembly]}=" + vcf_df["CHROM"] + ":" + vcf_df["POS"].astype(str)
    )
    vcf_df.drop(
        columns={"CHROM", "POS", f"REF_{new_assembly}", f"ALT_{new_assembly}"},
        inplace=True,
    )  # REF=REF_NEW, redundant. Should be removed in previous step
    logger.info("Old assembly coords added to INFO field")
    return vcf_df


def reorder_lifted_2(data: pd.DataFrame, new: str, old: str) -> pd.DataFrame:
    """
    Modify the DataFrame so that CHROM_hg38 and POS_hg38 become the new CHROM and POS,
    and the old CHROM and POS are renamed to CHROM_hg19 and POS_hg19 (or any specified old assembly).

    Args:
        - data: DataFrame containing 'CHROM', 'POS', 'CHROM_hg38', 'POS_hg38'.
        - new: Name of the new assembly (e.g., 'hg38').
        - old: Name of the old assembly (e.g., 'hg19').

    Returns:
        - Modified DataFrame with updated column names.
    """

    data.rename(columns={"CHROM": f"CHROM_{old}", "POS": f"POS_{old}"}, inplace=True)

    # Rename CHROM_hg38 and POS_hg38 to CHROM and POS
    data.rename(columns={f"CHROM_{new}": "CHROM", f"POS_{new}": "POS"}, inplace=True)

    data["INFO"] = (
        f"{ASSEMBLY[old]}="
        + data[f"CHROM_{old}"]
        + ":"
        + data[f"POS_{old}"].astype(str)
    )
    data.drop(
        columns={f"CHROM_{old}", f"POS_{old}", f"REF_{new}", f"ALT_{new}"}, inplace=True
    )  # REF=REF_NEW, redundant
    logger.info("Old assembly coords added to INFO field")
    logger.debug(f"Columns in df are {data.columns}")

    return data

def parse_samples(samples: list, format: list) -> dict:
    """Parse samples from list to dict
    RETURNS:
    - dict: {sample_name: sample_values (0/1:20)}
    """
    sample_dict = {}
    for i,sample in enumerate(samples):
        sample_name, sample_values_sep = sample.split(":")
        sample_values = sample_values_sep.split("|")
        if len(sample_values) != len(format.split(":")):
            logger.error(
                f"Sample {sample_name} has {len(sample_dict[sample_name])} values, \
                but format has {len(format.split(':'))} values"
            )
            raise ValueError(
                f"Sample {sample_name} has {len(sample_dict[sample_name])} values, \
                but format has {len(format.split(':'))} values"
            )
        else:
            sample_dict[sample_name] = sample_values
    logger.info(f"Parsed samples to dictionary: {sample_dict}")
    return sample_dict

def write_as_vcf(
    data: pd.DataFrame,
    output_file: str,
    override: bool,
    old_assembly: str,
    new_assembly: str,
    info_cols: str,
    format_col: str,
    samples: list = [],
):
    """Parse dataframe and saves an vcf

    Params:
    - data: dataframe
    - output_file: path to output
    - override: indicates if dataframe has old and new coords(FALSE) or only
    new coords (TRUE)
    - old_assembly: original assembly, only used if `override`=TRUE
    - new_assembly: liftover target assembly, only used if `override`=TRUE
    TODO: sort and compress (bgzip)
    """

    vcf_df = data.drop(columns={"END"}, inplace=False, errors="ignore")  # type:pd.DataFrame
    vcf_df = vcf_df.drop(columns={"END_hg38"}, inplace=False, errors="ignore")  # type:pd.DataFrame

    vcf_cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    # INSERT COLUMNS
    try:
        vcf_df.insert(2, vcf_cols[2], ".")
    except ValueError:
        logger.warning("ID column already exists, using existing one")
    try:
        vcf_df.insert(5, vcf_cols[5], ".")
    except ValueError:
        logger.warning("QUAL column already exists, using existing one")
    try:
        vcf_df.insert(6, vcf_cols[6], "PASS")
    except ValueError:
        logger.warning("FILTER column already exists, using existing one")
    try:
        vcf_df.insert(7, vcf_cols[7], ".")
    except ValueError:
        logger.warning("INFO column already exists, using existing one")

    if format_col and samples:
        i = 8
        vcf_df.insert(i, "FORMAT", format_col)  # not fine
        vcf_cols.insert(i,"FORMAT")
        logger.info(f"Added FORMAT field = {format_col}")
        # append to list in given position
        samples_dict = parse_samples(samples,format_col)
        insert_sample_values(vcf_df, samples_dict, vcf_cols)
    elif not format_col and samples:
        logger.error(
            "FORMAT column not specified")
        raise ValueError("FORMAT column not specified. Use --format_column")
    elif not samples and format_col:
        logger.error(
            "SAMPLE column not specified")
        raise ValueError("SAMPLE column not specified. Use --samples")

    len_cols = len(vcf_cols)
    # reorder column to the right position, leeve the rest as is
    vcf_df = vcf_df.reindex(columns=vcf_cols + vcf_df.columns[len_cols:].to_list())
    logger.debug(f"Columns in vcf are {vcf_df.columns}")
    # CHANGE chr TODO
    vcf_df["CHROM"] = vcf_df["CHROM"].str.replace("chr", "")
    logger.info("Removing `chr` string from CHROM column")

    if not override:
        vcf_df[f"CHROM_{new_assembly}"] = vcf_df[f"CHROM_{new_assembly}"].str.replace(
            "chr", ""
        )

        vcf_df = reorder_lifted_2(vcf_df, new=new_assembly, old=old_assembly)

    # colist = vcf_df.columns.to_list()
    # info_index = colist.index('INFO')
    # logger.debug(f'INFO is in {info_index} and POS_hg38 is in {colist.index("POS_hg38")}')
    # add info_cols to new INFO field
    nondropcols = []
    if info_cols:
        logger.info(
            f"Adding {info_cols} to INFO field.")

        vcf_df["INFO"] = ""
        nondropcols = info_cols.split(",")
        for i, column in enumerate(nondropcols):
            # convert to string TODO: review the type in the header and adjust accordingly instead of changing the type
            vcf_df.loc[:, column] = vcf_df.loc[:, column].fillna(".").astype(str)
            # replace invalid characters
            vcf_df.loc[:, column] = vcf_df.loc[:, column].astype(str).str.replace(" ", "_")
            # vcf_df.loc[:, column] = vcf_df.loc[:, column].str.replace(':',',')
            vcf_df.loc[:, column] = vcf_df.loc[:, column].astype(str).str.replace(";", ",")
            if vcf_df.columns.duplicated().any():
                # log wich columns are duplicated
                logger.debug(f"Duplicated columns {vcf_df.columns[vcf_df.columns.duplicated()]}")
                logger.warning(f"Duplicate column names detected. Removing duplicates.")
                vcf_df = vcf_df.loc[:, ~vcf_df.columns.duplicated()]
            vcf_df.loc[:, "INFO"] = (
                vcf_df.loc[:, "INFO"] + f"{column}=" + vcf_df.loc[:, column] + ";"
            )
        # remove the last ";"
        logger.debug(vcf_df["INFO"].head())
        vcf_df["INFO"] = vcf_df.loc[:,"INFO"].astype(str).str[:-1]

    logger.debug(f"Current columns in vcf are {vcf_df.columns.to_list()}")

    # ordering and dropping columns
    vcf_df = vcf_df[vcf_cols]
    logger.info(f"""Droping all columns but {vcf_df.columns}\n
{nondropcols} included in INFO field""")

    # sorting
    logger.debug(vcf_df.columns[:2])
    vcf_df.sort_values(
        by=vcf_df.columns[:2].to_list(), key=natsort_keygen(), axis=0, inplace=True
    )
    logger.info("values sorted by CHROM and POS")
    logger.info(f"Saving file to {output_file}")
    with open(output_file, "w") as f:
        f.write("""##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed (placeholder)">
##contig=<ID=1,length=249167691>
##contig=<ID=2,length=242695901>
##contig=<ID=3,length=197800245>
##contig=<ID=4,length=190915651>
##contig=<ID=5,length=180666277>
##contig=<ID=6,length=170877445>
##contig=<ID=7,length=159086805>
##contig=<ID=8,length=146293415>
##contig=<ID=9,length=141018424>
##contig=<ID=10,length=135434552>
##contig=<ID=11,length=134938471>
##contig=<ID=12,length=133763353>
##contig=<ID=13,length=115045730>
##contig=<ID=14,length=107285438>
##contig=<ID=15,length=102369712>
##contig=<ID=16,length=90141356>
##contig=<ID=17,length=81006630>
##contig=<ID=18,length=78014583>
##contig=<ID=19,length=59071322>
##contig=<ID=20,length=62906515>
##contig=<ID=21,length=48077813>
##contig=<ID=22,length=51156934>
##contig=<ID=X,length=154847490>
##contig=<ID=Y,length=62460029>
##contig=<ID=MT,length=16569>\n""")
        if not override:
            f.write(
                f'##INFO=<ID={ASSEMBLY[old_assembly]},Number=.,Type=String,Description="Coordinates in {old_assembly}">\n'
            )
        # f.write('##INFO=<ID=OTHER,Number=.,Type=String,Description="Information included in original dataframe (maybe genotype) ">\n')
        if info_cols:
            for info in info_cols.split(","):
                f.write(
                    f'##INFO=<ID={info},Number=.,Type=String,Description="Column {info} from original dataframe">\n'
                )
                logger.info(f"Adding {info} to INFO header as type String")
        if format_col:
            format_header = parse_format(format_col)
            f.write(format_header)
            logger.info("FORMAT field added to HEADER")
        f.write("#" + "\t".join(vcf_cols) + "\n")
        vcf_df.to_csv(f, sep="\t", index=False, header=False)
    logger.debug(vcf_df.shape[0])
    logger.info(f"VCF file created succesfully at {output_file}")

def insert_sample_values(vcf_df, samples_dict, vcf_cols):
    """Insert sample values in the dataframe
    Args:
    - vcf_df: dataframe to insert values
    - samples_dict: dictionary with sample names and values
    - vcf_cols: columns in the dataframe
    """
    for sample,values in samples_dict.items():
        vcf_df[sample] = ""
        try:
            vcf_df[sample] = vcf_df[values].apply(
                    lambda x: ':'.join(x.dropna().astype(str)),
                    axis=1
                )
            vcf_cols.append(sample)
            logger.info(f"Added {sample} to dataframe")
        except KeyError as e:
            logger.error(f"Error processing {sample}: {e}")
            raise KeyError(f"Error processing {sample}: {e}")
    
    return vcf_df,vcf_cols


def lifto(
    data: str | pd.DataFrame,
    old: str = "hg19",
    new: str = "hg38",
    column_mapping: list[str] = [],
    samples: list[str] = [],
    one_based: bool = False,
    override_coords: bool = True,
    lift_end: bool = True,
    show_err: bool = False,
    output_file: Literal["dict", "df", "csv"] = "df",
    remove_missing: bool = True,
    info_cols: str = "",
    format_col: str = "",
) -> pd.DataFrame | dict | None:
    """
    MAIN

    Lift the coordinates between assemblies using liftover library.
    REQUIRED COLUMNS: CHROM, POS

    Params:
    - data: path to csv file or dataframe
    - old: older assembly (def. hg19)
    - new: new assembly (def. hg38)
    - column_mapping: map columns to required columns
    - one_based: 1-based (True) o 0-based start BUG
    - append_coords: Override existing coordinates
    - show_err: Show what variants could not be converted to stdout
    - output_file: file name
    - TODO remove_missing: removes from dataframe non converted coords

    Returns:
    - pd.DataFrame | dict | None
    """
    # read File if it's not a dataframe (dataframe processing should be deprecated)
    if isinstance(data, str):
        sep = "\t" if data.endswith(".tsv") else ","
        pre_data = pd.read_csv(data, header=0, sep=sep)
    else:
        pre_data = data

    columns = ["CHROM", "POS", "END", "REF", "ALT"]
    present_cols = check_columns_present(df=pre_data,
                                         required_cols=columns,
                                         required_are=4,
                                         column_mapping=column_mapping)
    if present_cols == columns:
        check_coords(pre_data)

    lifter = get_lifter(old, new, one_based=one_based)
    result_dict, error = lift_coords(pre_data.loc[:, present_cols], lifter)
    logger.debug(list(result_dict.values())[:2])
    total = pre_data.shape[0]
    exito = len(result_dict)
    logger.info(
        f"Successfully converted: \
{len(pre_data)} - {round(exito / total, 2) * 100} %."
    )
    if show_err and error:
        save_errors(error, output_file)

    if output_file.endswith(("df", "csv", "vcf")):
        # dict to dataframe wo proper colnames
        result = pd.DataFrame(result_dict).transpose()
        result.reset_index(inplace=True)
        logger.debug(result)
        keys_as_cols = present_cols
        values_as_cols = [f"{x}_{new}" for x in present_cols]
        logger.debug(f"{keys_as_cols}\n{values_as_cols}")
        # set colnames hg19coord_names + hg38coords_names
        try:
            result.columns = keys_as_cols + values_as_cols
        except ValueError as e:
            print(
                f"{e}\ncolumns={result.columns}\n\
newcols={keys_as_cols + values_as_cols}"
            )
            raise ValueError

        enforce_col_types(result)

        # "append" data to the results
        # no need to check as non converted won't be used
        # can be save in a separate file
        result_merged = pd.merge(result, pre_data, how="left", on=present_cols)
        logger.debug(f"{result_merged.shape}")

        if override_coords:
            # new coords overwritten
            logger.debug(result_merged.head())
            result_merged[keys_as_cols] = result_merged[values_as_cols]
            # remove duplicate cols -> *_hg38 as they are copied to *
            result_merged.drop(columns=values_as_cols, inplace=True)
            cols = result_merged.columns.to_list()
            logger.debug(f"{cols}")

            # change colnames from CHROM_hgXX to CHROM
            # result_merged.columns = map(lambda x: x.replace(f'_{new}', ''),
            #                             result_merged.columns.to_list())

        if output_file.endswith("df"):
            return result_merged
        else:
            create_dir(output_file)
            if output_file.endswith("vcf"):
                write_as_vcf(
                    data=result_merged,
                    output_file=output_file,
                    old_assembly=old,
                    new_assembly=new,
                    override=override_coords,
                    info_cols=info_cols,
                    format_col=format_col,
                    samples=samples,
                )
            else:
                result_merged.to_csv(output_file, index=False)

            logger.info(f"File saved in {output_file}")

    else:
        logger.warning(f"""{len(error)} positions could not be converted
                       and will not be returned.
                       To see which positions are missing enable `--show_err`
                       Or choose other `--output_type` """)
        return result_dict


def argparsing(args: list[str] | None):
    """
    Parse args either from cli or passed as arguments"""

    parser = argparse.ArgumentParser()

    parser.add_argument("--input", type=str, help="path to data")
    parser.add_argument("--output", type=str, help="Output file name")
    parser.add_argument(
        "--old", type=str, help="assembly to convert. Default=hg19", default="hg19"
    )
    parser.add_argument(
        "--new", type=str, help="new assembly. Default=hg38", default="hg38"
    )
    parser.add_argument(
        "--column_mapping",
        type=str,
        help="Map columns in data to required columns. \
        Eg: --column_mapping CHROMOSOME:CHROM,POSITION:POS,REFERENCE:REF,ALT:ALT")
    parser.add_argument(
        "--one",
        type=bool,
        help="""Are the coordinates one based? Help: If an SNV has the same start and end is one based IF an insertion end = start+n_nucleotides is one_based Else is zero-based""",
        default=False,
    )
    parser.add_argument(
        "--override_coords",
        action="store_true",
        help="Do you want to override the original coords?",
    )
    parser.add_argument(
        "--show_err", action="store_true", default=False, help="show errors"
    )
    parser.add_argument(
        "--lift_end", action="store_true", help="do you have end coord?"
    )
    parser.add_argument(
        "--add_info_fields",
        type=str,
        help="Only if the output is vcf. \
Which columns do you want to include in the INFO field \
Eg: --add_info_fields samplegenotype,quality will include these two columns",
    )
    parser.add_argument(
        "--format_column",
        type=str,
        help="Indicate FORMAT field to use. Only for VCF output. Eg: --format_column GT:GQ:DP",
    )
    parser.add_argument(
        "--samples",
        type=str,
        help="Sample name and values to use as the header of the SAMPLE field. Eg: --sample_name SAMPLE1:GT_colname|GQ_colname,SAMPLE2:GT_colname2|GQ_colname2. This column must have the same number of values for each sample as the --format_column.",
    )
    parser.add_argument("--logging_level", type=str, default="INFO", help="Logging level")

    return parser.parse_args(args)


def main(args: list[str] | None):
    parser = argparsing(args)
    # set logging level
    logger.setLevel(getattr(logging, parser.logging_level.upper(), logging.INFO))
    lifto(
        data=parser.input,
        old=parser.old,
        new=parser.new,
        column_mapping=parser.column_mapping.split(",") if parser.column_mapping else [],
        samples=parser.samples.split(",") if parser.samples else [],
        override_coords=parser.override_coords,
        output_file=parser.output,
        show_err=parser.show_err,
        lift_end=parser.lift_end,
        info_cols=parser.add_info_fields,
        format_col=parser.format_column,
    )


if __name__ == "__main__":
    main(None)
