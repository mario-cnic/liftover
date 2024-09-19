import pandas as pd
from liftover import get_lifter
import logging
import sys
from typing import Literal
import os
import argparse


def check_coords(pre_data: pd.DataFrame, logger: logging.Logger) -> bool:
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
    nuc = ['A', 'C', 'G', 'T']
    snvs = pre_data[pre_data['REF'].isin(nuc) & pre_data['ALT'].isin(nuc)]
    if all(snvs['POS'] == snvs['END']):
        logger.debug('1-based')
        return True
    else:
        logger.debug('0-based')
        return False


def check_conversion(start, end, ref, alt, logger):
    """
    Check new POS and END agree with REF/ALT lengths and
    coordinate system type.
    TODO
    """
    pass


def enforce_col_types(data: pd.DataFrame):
    data['CHROM'] = data['CHROM'].astype(str)
    data['POS'] = data['POS'].astype(int)


def check_columns_present(df: pd.DataFrame, required_cols: list[str],
                          logger: logging.Logger,
                          required_are: int = 2) -> list[str]:
    """
    Check columns are present in dataframe
    Args:
    - df
    - required_columns
    - required_are: posición de las columnas esenciales
    las que hacen saltar error

    Raises:
    - IndexError: Missing columns

    Return:
    - Columns present in dataframe
    """
    essentials = required_cols[:required_are]
    missing = [c for c in required_cols if c not in df.columns]
    present = [c for c in required_cols if c in df.columns]
    if any([x in missing for x in essentials]):
        raise IndexError(f"The following columns are missing\
{', '.join(missing)}")
    elif missing:
        logger.warning(f'''The following columns are missing {missing}.
The liftover coordinates will not be checked against
REF-ALT lengths''')
    enforce_col_types(df)
    logger.debug(present)
    return present


def lift_coords(pre_data: pd.DataFrame, lifter,
                logger: logging.Logger) -> tuple[dict[tuple:tuple],
                                                 list[tuple]]:
    """
    Iterates over a dataframe and lift the coordinates using `lifter` object

    Returns:
    - tuple: Dictionary of positions and list of errors (empty if None)
    """
    post_dict = dict()
    error = list()
    n = 0

    is_ref_in_col = True if 'REF' in pre_data.columns else False
    is_alt_in_col = True if 'ALT' in pre_data.columns else False
    is_end_in_cols = True if 'END' in pre_data.columns else False

    for _, row in pre_data.iterrows():
        keys = tuple(row.loc[pre_data.columns].to_list())
        try:
            chrom, pos = lifter.convert_coordinate(str(row.CHROM),
                                                   row.POS)[0][:2]
            if is_end_in_cols:
                end = lifter.convert_coordinate(str(row.CHROM),
                                                row.END)[0][1]
            else:
                end = 0
            if is_ref_in_col:
                if is_alt_in_col:
                    post_dict[keys] = (chrom, pos, end, row.REF, row.ALT)\
                        if is_end_in_cols else (chrom, pos, row.REF, row.ALT)
                else:
                    post_dict[keys] = (chrom, pos, end, row.REF)\
                        if is_end_in_cols else (chrom, pos, row.REF)
            elif is_alt_in_col:
                post_dict[keys] = (chrom, pos, end, row.ALT)\
                    if is_end_in_cols else (chrom, pos, row.ALT)
            else:
                post_dict[keys] = (chrom, pos, end) if end != 0 else (chrom,
                                                                      pos)
        except IndexError:
            error.append(tuple(keys))
            n = n+1
        except KeyError:
            error.append(tuple(keys))
            logger.debug(f'Error processing {row.CHROM}, jumping next...')
            n = n+1
    # logger.debug(post_dict)
    return post_dict, error


def save_errors(data: dict, show_err: bool,
                error: list, logger):
    """Save not lifted variants.
    """
    if show_err and error:
        err_path = 'out/errors.txt'
        os.mkdir('out/') if not os.path.isdir('out/') else 0
        logger.info(f'Failed a total of {len(error)},\
                    saving file in {err_path}')
        with open(err_path, 'w') as f:
            f.write(f"{','.join(data)}\n")
            for value in error:
                f.write(f"{','.join(map(str, value))}\n")


def create_dir(f):
    directory = os.path.dirname(f)
    if not os.path.exists(directory):
        os.makedirs(directory)


def lifto(data: str | pd.DataFrame, old: str = 'hg19', new: str = 'hg38',
          one_based: bool = False, override_coords: bool = True,
          lift_end: bool = True, show_err: bool = False,
          output_file: Literal['dict', 'df', 'csv'] = 'df',
          remove_missing: bool = True) -> (pd.DataFrame | dict | None):
    """
    MAIN

    Lift the coordinates between assemblies using liftover library.
    REQUIRED COLUMNS: CHROM, POS

    Params:
    - data: path to csv file or dataframe
    - old: older assembly (def. hg19)
    - new: new assembly (def. hg38)
    - one_based: 1-based (True) o 0-based start BUG
    - append_coords: Override existing coordinates
    - show_err: Show what variants could not be converted to stdout
    - output_file: file name
    - TODO remove_missing: removes from dataframe non converted coords

    Returns:
    - pd.DataFrame | dict | None
    """
    # logging definition -> Maybe this goes to main(?)
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    # añadimos un handler para stdout
    c_handler = logging.StreamHandler(sys.stdout)
    c_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
    c_handler.setFormatter(c_format)
    logger.addHandler(c_handler)

    if isinstance(data, str):
        pre_data = pd.read_csv(data, header=0, sep=',')
    else:
        pre_data = data

    columns = ['CHROM', 'POS', 'END', 'REF', 'ALT']
    present_cols = check_columns_present(pre_data, columns, logger)
    if present_cols == columns:
        check_coords(pre_data, logger)

    lifter = get_lifter(old, new, one_based=one_based)
    result_dict, error = lift_coords(pre_data.loc[:, present_cols],
                                     lifter, logger)
    logger.debug(list(result_dict.values())[:2])
    total = pre_data.shape[0]
    exito = len(result_dict)
    logger.info(f'Successfully converted: \
{len(pre_data)} - {round(exito/total, 2)*100} %.')

    save_errors(result_dict, show_err, error, logger)

    if output_file.endswith(('df', 'csv', 'vcf')):
        # dict to dataframe wo proper colnames
        result = pd.DataFrame(result_dict).transpose()
        result.reset_index(inplace=True)
        logger.debug(result)
        keys_as_cols = present_cols
        values_as_cols = [f'{x}_{new}' for x in present_cols]
        logger.debug(f'{keys_as_cols}\n{values_as_cols}')
        # set colnames hg19coord_names + hg38coords_names
        try:
            result.columns = keys_as_cols + values_as_cols
        except ValueError as e:
            print(f'{e}\ncolumns={result.columns}\n\
newcols={keys_as_cols + values_as_cols}')
            raise ValueError

        enforce_col_types(result)

        # "append" data to the results
        # no need to check as non converted won't be used
        # can be save in a separate file
        result_merged = pd.merge(result, pre_data, how='left', on=present_cols)
        logger.debug(f'{result_merged.shape}')

        if override_coords:

            # new coords overwritten
            logger.debug(result_merged.head())
            result_merged[keys_as_cols] = result_merged[values_as_cols]
            # remove duplicate cols -> *_hg38 as they are copied to *
            result_merged.drop(columns=values_as_cols, inplace=True)
            cols = result_merged.columns.to_list()
            logger.debug(f'{cols}')

            # change colnames from CHROM_hgXX to CHROM
            # result_merged.columns = map(lambda x: x.replace(f'_{new}', ''),
            #                             result_merged.columns.to_list())

        if output_file.endswith('df'):
            return result_merged
        else:
            create_dir(output_file)
            if output_file.endswith('vcf'):
                vcf_cols = [f'CHROM_{new}', f'POS_{new}', 'ID', 'REF', 'ALT', 'QUAL',
                            'FILTER', 'INFO']
                # FIX POSITIONS <-
                result_merged.insert(2, 'ID', '.')
                result_merged.insert(6, 'QUAL', '.')
                result_merged.insert(7, 'FILTER', '.')
                result_merged.insert(8, 'INFO', '.') 
                vcf_df = result_merged.drop(columns={'END'}, inplace=False)  # type:pd.DataFrame
                if not override_coords:
                    vcf_df['INFO'] = f'GRC{old}='+vcf_df['CHROM'] + ':'+vcf_df['POS'].astype(str) +\
                          '|'.join(vcf_df.columns.tolist()[9:]) +\
                          vcf_df.loc[:, ~vcf_df.columns.isin(vcf_cols)].astype(str).agg(';'.join, axis=1)
                    vcf_df.drop(vcf_df.iloc[:, 9: vcf_df.shape[1]+1], axis=1, inplace=True)
                with open(output_file, 'w') as f:
                    f.write('##fileformat=VCFv4.1\n#')
                    f.write('\t'.join(vcf_cols)+'\n')
                    vcf_df.to_csv(f, sep='\t', index=False, header=False)
            else:
                result_merged.to_csv(output_file, index=False)

            logger.info(f'File saved in {output_file}')

    else:
        logger.warning(f'''{len(error)} positions could not be converted
                       and will not be returned.
                       To see which positions are missing enable `--show_err`
                       Or choose other `--output_type` ''')
        return result_dict


def argparsing(args: list[str] | None):
    """
    Parse args either from cli or passed as arguments"""

    parser = argparse.ArgumentParser()

    parser.add_argument('--in_file', type=str,
                        help='path to data')
    parser.add_argument('--old', type=str,
                        help='assembly to convert', default='hg19')
    parser.add_argument('--new', type=str,
                        help='new assembly', default='hg38')
    parser.add_argument('--one', type=bool,
                        help='''are the coordinates one based?
Help:
If an SNV has the same start and end is one based
IF an insertion end = start+n_nucleotides is one_based
                        Else is zero-based''',
                        default=False)
    parser.add_argument('--override_coords',
                        action='store_true',
                        help='Do you want to override\
the original coords?',
                        default=False)
    parser.add_argument('--show_err', type=bool,
                        help='print errors', default=False)
    parser.add_argument('--output', type=str,
                        help='''''')
    parser.add_argument('--lift_end', type=bool,
                        help='do you have end coord?', default=False)

    return parser.parse_args(args)


def main(args: list[str] | None):
    parser = argparsing(args)

    lifto(parser.in_file, old=parser.old, new=parser.new,
          override_coords=parser.override_coords,
          output_file=parser.output, show_err=parser.show_err,
          lift_end=parser.lift_end)


if __name__ == '__main__':
    main(None)
