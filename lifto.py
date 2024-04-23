import pandas as pd
from liftover import get_lifter
import logging
import sys
from typing import Literal
import os
import re
import argparse


def check_coords(pre_data: pd.DataFrame, logger: logging.Logger):
    """
    Comprueba si las coordenadas son 1-based o 0-based
    En función de como lee los SNVs

    Las 1-based
    - SNVs --> START = END
    - Inserciones --> START = END
    - Delecciones --> START = END - len(REF) + 1

    Las 0-based
    - SNVs --> START = END - 1
    - Inserciones --> START = END
    - Delecciones --> START = END - len(REF)
    """
    nuc = ['A', 'C', 'G', 'T']
    snvs = pre_data[pre_data['REF'].isin(nuc) & pre_data['ALT'].isin(nuc)]
    if all(snvs['POS'] == snvs['END']):
        logger.debug('Es 1-based')
    else:
        logger.debug('Es 0-based')


def check_columns_present(df: pd.DataFrame, required_cols: list,
                          logger: logging.Logger) -> bool:
    """
    Comprueba que `required_cols` están presentes en el dataframe
    Args:
        df
        required_columns

    Raises:
        IndexError: Missing columns

    Return:
        bool
    """

    missing = [c for c in required_cols if c not in df.columns]
    if 'CHROM' in missing or 'POS' in missing:
        raise IndexError(f"Faltan las columnas {','.join(missing)}")
    elif missing:
        logger.warning(f'''Faltan las columnas {missing}. No se hará el
                        checkeo de 1-based o 0-based''')
        return False
    else:
        return True


def lifto(data: str | pd.DataFrame, old: str = 'hg19', new: str = 'hg38',
          one_based: bool = False, append_coords: bool = True,
          lift_end: bool = True, show_err: bool = False,
          output_type: Literal['dict', 'df', 'csv'] = 'df',
          remove_missing: bool = True) -> (pd.DataFrame | dict | None):
    """
    Convierte entre ensamblajes de genomas para una lista de SNV
    usando la librería liftover y devuelve un dataframe o diccionario

    Params:
    - data: ruta del archivo a convertir o dataframe
            El archivo debe tener formato .csv (comma separated values)
            Columnas obligadas -> `CHROM`, `POS`
    - old: assembly del que proviene las variantes (def. hg19)
    - new: assembly destino (def. hg38)
    - one_based: 1-based (True) o 0-based start
    - append_coords: Si `True` crea el mismo archivo, cambiando las coordenadas.
        Si `False` añade al final del archivo las nuevas coordenadas
    - show_err: Mostrar errores en standard output
    - append_to_pd: Devolver los datos al pandas dataframe
    - remove_missing: Elimina las variantes que no han sido convertidas del
        dataframe (def. True)

    Returns:
    - pd.DataFrame | dict : dataframe o diccionario con la posición en el
        nuevo assembly.
        El diccionario tiene como keys las posiciones del assembly anterior
    """
    # logging definition
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    # añadimos un handler para stdout
    c_handler = logging.StreamHandler(sys.stdout)
    c_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
    c_handler.setFormatter(c_format)
    logger.addHandler(c_handler)

    lifter = get_lifter(old, new, one_based=one_based)

    if isinstance(data, str):
        pre_data = pd.read_csv(data, header=0, sep=',')
    else:
        pre_data = data

    post_dict = dict()
    error = list()
    n = 0
    columns = ['CHROM', 'POS', 'END', 'REF', 'ALT']
    if check_columns_present(pre_data, columns, logger):
        check_coords(pre_data, logger)

    if 'END' not in pre_data.columns and not lift_end:
        pre_data['END'] = 0

    for _, row in pre_data.iterrows():
        try:
            chrom, pos = lifter.convert_coordinate(str(row.CHROM),
                                                   row.POS)[0][:2]
            if lift_end:
                end = lifter.convert_coordinate(str(row.CHROM),
                                                row.END)[0][1]
            else:
                end = 0
            post_dict[row.CHROM, row.POS, row.END] = chrom, pos, end
        except IndexError:
            error.append((row.CHROM, row.POS, row.END))
            n = n+1
        except KeyError:
            error.append((row.CHROM, row.POS, row.END))
            logging.debug(f'Error al procesar {row.CHROM}, saltando al siguiente...')
            n = n+1  # sobra
    exito = len(post_dict)
    total = pre_data.shape[0]
    logger.info(f'Convertidas con éxito {len(post_dict)} - {round(exito/total,2)*100} %.')
    if n > 0:
        logger.warning(f'Han fallado un total de {n}')
    if show_err and n > 0:
        os.mkdir('out/') if not os.path.isdir('out/') else 0
        logger.info(f'Han fallado un total de {n} conversiones, guardando en out/')
        with open("out/errors.txt", 'w') as f:
            for value in error:
                f.write(f'{value}\n')

    if output_type in ['df', 'csv']:
        result = pd.DataFrame(post_dict).transpose()
        result.reset_index(inplace=True)
        result.columns = ['CHROM', 'POS', 'END',
                          f'CHROM{new}', f'POS{new}', f'END{new}']
        result_merged = pd.merge(pre_data, result, how='left', on=['CHROM',
                                                                   'POS',
                                                                   'END'])
        if not append_coords:
            result_merged[['CHROM', 'POS', 'END']] = result_merged[[f'CHROM{new}',
                                                                    f'POS{new}',
                                                                    f'END{new}']]
            result_merged.drop(columns=[f'CHROM{new}', 
                                        f'POS{new}', f'END{new}'], inplace=True)
        if output_type == 'df':
            return result_merged
        else:
            os.mkdir('out/') if not os.path.isdir('out/') else 0
            filename = os.path.basename(data) if isinstance(data, str) else 'dataframe'
            filename_wo_ex = re.sub(r'\.\w+$', '', filename)
            out_path = f'out/{old}_to_{new}_from_{filename_wo_ex}.csv'
            result_merged.to_csv(out_path, index=False)
            logger.info(f'Archivo guardado en {out_path}')

    else:
        logger.warning(f'Las variantes no convertidas n = {n} no se han devuelto')
        return post_dict


def parse_args():
    parser = argparse.ArgumentParser(
    )
    parser.add_argument('--in_file', type=str, help='path to data')
    parser.add_argument('--old', type=str, help='assembly to convert', default='hg19')
    parser.add_argument('--new', type=str, help='new assembly', default='hg38')
    parser.add_argument('--one', type=bool, help='''are the coordinates one based?
                        Help:
                        If an SNV has the same start and end is one based
                        IF an insertion end = start+n_nucleotides is one_based
                        Else is zero-based''', default=False)
    parser.add_argument('--append_coords', type=bool,
                        help='When returning the data, do you want the substitute the original coord?',
                        default=False)
    parser.add_argument('--show_err', type=bool, help='print errors', default=False)
    parser.add_argument('--output_type', type=str, help='print errors', default='csv')
    parser.add_argument('--lift_end', type=bool, help='do you have end coord?',default=False)

    return parser.parse_args()


if __name__ == '__main__':
    parser = parse_args()
    lifto(parser.in_file, old=parser.old, new=parser.new,
          append_coords=parser.append_coords,
          output_type=parser.output_type, show_err=parser.show_err, lift_end=parser.lift_end)
