# Indel Mapper

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

## Set up

Set up your environment:

```shell
virtualenv indel-mapper
source indel-mapper/bin/activate
pip3 install -r requirements.txt
```

## Web application

Run the application:

```shell
python3 app.py
```

Run the tests:

```shell
python3 -m unittest
```

Deploy to Heroku:

```shell
git push heroku master
```

## Command line application

Example:

```shell
python3 run.py -a ~/Documents/bowtie2_results.sam -r ~/Documents/references.csv -o ~/Documents/results.csv
```

The arguments are:

* `-a` or `--alignment` Alignment SAM file
* `-r` or `--reference` Reference CSV file
* `-o` or `--output` Filename to save the output as

## License

Indel Mapper is licensed under Version 3 of the GNU General Public License.
