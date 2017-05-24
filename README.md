# Indel Mapper

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

## Set up

Set up your environment:

```shell
virtualenv indel-mapper
source indel-mapper/bin/activate
pip3 install -r requirements.txt
```

Install Redis

## Web application

Run the application:

```shell
python3 app.py
```

Run Redis

```shell
src/redis-server
```

Run Celery

```shell
indel-mapper/bin/celery worker -A app.celery --loglevel=info
```

## Testing

Run the tests:

```shell
python3 -m unittest
```

## Deploying

Deploy to Heroku:

```shell
git push heroku master
```

## Command line application

Example:

```shell
python3 run.py -a ~/Documents/bowtie2_results.sam -r ~/Documents/references.csv -o ~/Documents/results.csv
```
There are three required arguments:

* `-a` or `--alignment` Alignment SAM file
* `-r` or `--reference` Reference CSV file
* `-o` or `--output` Output file, in CSV

## License

Indel Mapper is licensed under Version 3 of the GNU General Public License.
