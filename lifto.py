import pandas as pd
from liftover import get_lifter
import logging
import sys
from typing import Literal
import os
import argparse

# logging definition -> Maybe this goes to main(?)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
# añadimos un handler para stdout
c_handler = logging.StreamHandler(sys.stdout)
c_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
c_handler.setFormatter(c_format)
logger.addHandler(c_handler)

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


def lift_coords(pre_data: pd.DataFrame, lifter) -> tuple[dict[tuple:tuple],
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
            logger.debug(f'{chrom}:{row.POS}- new:{pos}')
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
                error: list):
    """Save not lifted variants.
    """
    if show_err and error:
        err_path = 'out/errors.txt'
        cwd = os.getcwd()
        os.mkdir(cwd+'out/') if not os.path.isdir('out/') else 0
        logger.info(f'Failed a total of {len(error)},\
                    saving file in {err_path}')
        with open(err_path, 'w') as f:
            f.write(f"{','.join(data)}\n")
            for value in error:
                f.write(f"{','.join(map(str, value))}\n")
        logger.info(f'Error file saved at {err_path}')


def create_dir(f):
    directory = os.path.dirname(f)
    if not os.path.exists(directory):
        os.makedirs(directory)


def write_as_vcf(data: pd.DataFrame, output_file: str,
                 override: bool,
                 old_assembly: str, new_assembly: str):
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
    vcf_df = data.drop(columns={'END'}, inplace=False)  # type:pd.DataFrame
    vcf_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
                'FILTER', 'INFO']
    # INSERT COLUMNS
    vcf_df.insert(2, 'ID', '.')
    vcf_df.insert(5, 'QUAL', '.')
    vcf_df.insert(6, 'FILTER', 'PASS')
    vcf_df.insert(7, 'INFO', '.')
    # CHANGE chr TODO
    vcf_df['CHROM'] = vcf_df['CHROM'].str.replace('chr', '')
    if not override:
        vcf_df[f'CHROM_{new_assembly}'] = vcf_df[f'CHROM_{new_assembly}'].str.replace('chr', '')
        vcf_df = vcf_df.drop(columns={'END_hg38'}, inplace=False)  # type:pd.DataFrame
        logger.debug(vcf_df.shape[0])
        # fix order (as CHROM_{new} is the target col, not CHROM)
        colist = vcf_df.columns.to_list()
        ch, ch_n = colist.index(f'CHROM_{new_assembly}'), colist.index('CHROM')
        colist[ch], colist[ch_n] = colist[ch_n], colist[ch]
        # same with POS <-> POS_{new}
        # logger.debug(vcf_df['POS'].head())
        # logger.debug(vcf_df[f'POS_{new_assembly}'].head())
        pos, pos_n = colist.index(f'POS_{new_assembly}'), colist.index('POS')
        colist[pos], colist[pos_n] = colist[pos_n], colist[pos]
        
        # assign new col order to df
        vcf_df = vcf_df.reindex(columns=colist)
        # logger.debug(vcf_df['POS'].head())
        # logger.debug(vcf_df[f'POS_{new_assembly}'].head())
        # insert other info in INFO field
        # logger.debug(f'pos_{old_assembly}:{vcf_df.iloc[2, 1]}')
        # logger.debug(f'pos_{new_assembly}:{vcf_df.loc[2,f"POS_{new_assembly}"]}')
        vcf_df['INFO'] = f'GRC{old_assembly}='+ vcf_df['CHROM'] + ':'+vcf_df['POS'].astype(str)
        vcf_df.drop(columns={'CHROM', 'POS',
                             f'REF_{new_assembly}', f'ALT_{new_assembly}'}, inplace=True)  # REF=REF_NEW, redundant. Should be removed in previous step
    colist = vcf_df.columns.to_list()
    info_index = colist.index('INFO')
    logger.debug(f'INFO is in {info_index} and POS_hg38 is in {colist.index("POS_hg38")}')
    vcf_df.loc[:,'INFO'] = vcf_df.loc[:, 'INFO'] + ';OTHER='+'|'.join(colist[info_index+1:]) +\
                           vcf_df.iloc[:, info_index+1:].astype(str).agg('|'.join, axis=1)
    vcf_df.drop(columns=vcf_df.columns[info_index+1:], inplace=True)
    logger.debug(vcf_df.columns.to_list())
    logger.debug(vcf_df.head())
    # drop other columns
    vcf_df.drop(vcf_df.iloc[:, len(vcf_cols): vcf_df.shape[1]+1], axis=1, inplace=True)
    with open(output_file, 'w') as f:
        f.write('''##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="placeholder, no info about the filtering">
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
##contig=<ID=MT,length=16569>\n''')
        if not override:
            f.write(f'##INFO=<ID=GRC{old_assembly},Number=.,Type=String,Description="Coordinates in {old_assembly}">\n')
        f.write('##INFO=<ID=OTHER,Number=.,Type=String,Description="Information included in original dataframe (maybe genotype) ">\n')
        f.write('#'+'\t'.join(vcf_cols)+'\n')
        vcf_df.to_csv(f, sep='\t', index=False, header=False)
    logger.debug(vcf_df.shape[0])
    logger.info('VCF file created succesfully')


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

    if isinstance(data, str):
        pre_data = pd.read_csv(data, header=0, sep=',')
    else:
        pre_data = data

    columns = ['CHROM', 'POS', 'END', 'REF', 'ALT']
    present_cols = check_columns_present(pre_data, columns)
    if present_cols == columns:
        check_coords(pre_data)

    lifter = get_lifter(old, new, one_based=one_based)
    result_dict, error = lift_coords(pre_data.loc[:, present_cols],
                                     lifter)
    logger.debug(list(result_dict.values())[:2])
    total = pre_data.shape[0]
    exito = len(result_dict)
    logger.info(f'Successfully converted: \
{len(pre_data)} - {round(exito/total, 2)*100} %.')

    save_errors(result_dict, show_err, error)

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
                write_as_vcf(data=result_merged,
                             output_file=output_file,
                             old_assembly=old,
                             new_assembly=new,
                             override=override_coords)
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
