## Descripción

Script para convertir entre coordenadas

## Requisitos

Para ejecutarlo requiere
- Instalar python (comprobar si está instalado con `python --version` o `python3 --version`)
- (Opcional) crear entorno virtual `python -m venv liftover_venv`
- (Opcional) activar entorno virtual `source liftover_venv/bin/activate`
- Instalar pip (comprobar si está instalado con `pip --version`)
- Ejecutar el comando `pip install requirements.txt`

## Comandos
- Ejecutar el script `python lifto.py -h` para ver los parámetros que acepta

  -h, --help            show this help message and exit
  --in_file IN_FILE     path to data
  --old OLD             assembly to convert
  --new NEW             new assembly
  --one ONE             are the coordinates one based? Help: If an SNV has the same start and end is one based IF an insertion end = start+n_nucleotides is
                        one_based Else is zero-based
  --append_coords APPEND_COORDS
                        When returning the data, do you want the substitute the original coord?
  --show_err SHOW_ERR   print errors
  --output_type OUTPUT_TYPE
                        print errors
  --lift_end LIFT_END   do you have end coord?

- Ej: `python lifto.py --in_file 'ruta/a/mi/archivo.csv' --output_type 'csv'`

Archivo de entrada -> archivo CSV con, al menos, los siguientes campos y con el siguiente nombre (cabecera)
- CHROM -> cromosoma, da igual el formato (acepta chr1 y 1)
- POS -> Posición o posición de inicio
- END -> posición final (si no se dispone y no se necesita se creará nueva columna)

Por defecto guarda los resultados dentro de una nueva carpeta *out/*

## Aviso

El script únicamente convierte las coordenadas de inicio `POS` y final `END` para el sistema de coordenadas basado en 1
ver ![esta explicación](https://www.biostars.org/p/84686/).

- SNVs no tienen rango (start,start+1) sino que son start = end (POS = END).
- Las inserciones se interpretan como POS = END, independientemente del nº de nucleótidos insertados
- Las delecciones, en cambio, tienen en cuenta el número de nucleótidos perdidos (start,start-len(ref)).

## TODO

- Transformar las INDELs a formato VCF **nota** las INDELs se indican con un `-`, sin embargo el standard VCF (usado por gnomAD) no emplea esta nomenclatura puesto que todos las variantes deben estar normalizadas y *left aligned*, por lo que es necesario anotar el nucleótido previo. https://www.biostars.org/p/356458/ 
- Una forma de llevar a cabo todo este proceso es hacer una anotación con VEP usando el formato propio de VEP http://www.ensembl.org/info/docs/tools/vep/vep_formats.html#default al parsear a este formato y anotar con VEP podemos obtener tanto sus coordenadas en el assembly 38 como el formato VCF


por Mario Ruiz Pérez